# ==============================================================================
# 04_Single-Cell-Mixed-Models.R
# ==============================================================================

# ─── LIBRARIES ───
suppressPackageStartupMessages({
  library(tidyverse)
  library(car)
  library(emmeans)
  library(effectsize)
  library(lme4)
  library(lmerTest)
  library(ggsignif)
  library(patchwork)
  library(openxlsx)
  library(broom.mixed)
})

# ─── COMMAND-LINE ARGS ───
args <- commandArgs(trailingOnly = TRUE)
SEG_DIR <- file.path("Confocal", "results_segmentation")   # outputs of the Python steps 02-03
master_cells_csv  <- file.path(SEG_DIR, "master_cells.csv")
master_coloc_csv  <- file.path(SEG_DIR, "master_coloc.csv")   # Costes p-values (step 03)
out_dir           <- file.path("Confocal", "results_confocal_stats")
exclude_csv       <- NULL

for (i in seq_along(args)) {
  if (args[i] == "--master-cells"  && i < length(args)) master_cells_csv  <- args[i + 1]
  if (args[i] == "--master-coloc"  && i < length(args)) master_coloc_csv  <- args[i + 1]
  if (args[i] == "--out-dir"       && i < length(args)) out_dir           <- args[i + 1]
  if (args[i] == "--exclude"       && i < length(args)) exclude_csv       <- args[i + 1]
}

out_dir_pdf  <- file.path(out_dir, "pdf")
out_dir_png  <- file.path(out_dir, "png")
out_dir_data <- file.path(out_dir, "data")
for (d in c(out_dir_pdf, out_dir_png, out_dir_data)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ─── DESIGN ───
virus_levels    <- c("None", "V1", "V2")
antibody_levels <- c("Syn1", "EPR20535")
group_order     <- c("Control", "V1", "V2")

my_comparisons_M1 <- list(c("Control", "V1"), c("Control", "V2"), c("V1", "V2"))

palette_grbr <- c("Control" = "#D5D8DC", "V1" = "#A8C7F0", "V2" = "#F2B4B4")
group_colors <- palette_grbr
antibody_pair_colors <- c("Syn1" = "#4285F4", "EPR20535" = "#E07B3D")
infection_colors <- c("Bystander" = "#9AA5B0", "Infected" = "#C0392B")
paracrine_colors <- c("Control" = "#D5D8DC", "V1.Bystander" = "#A8C7F0",
                      "V2.Bystander" = "#F2B4B4")

ADJUST_METHOD       <- "holm"
MIN_CELLS_FOR_LMER  <- 20
MIN_CELLS_PER_GROUP <- 5
MIN_CELLS_FOR_COLOC <- 15    
COSTES_P_THRESHOLD  <- 0.05    
OUTLIER_CLIP_K      <- 3
SHOW_ONLY_SIG       <- TRUE

# ─── BOOTSTRAP CONSTANTS & SEED ───
N_BOOTSTRAP       <- 5000
BOOTSTRAP_SEED    <- 20260524
BOOTSTRAP_MODULES <- c("M2", "M3B", "M4")
set.seed(BOOTSTRAP_SEED)

# ─── STRICT STATS: FISHER Z UTILITIES ───
fisher_z <- function(r) {
  r <- pmax(pmin(r, 0.9999), -0.9999) 
  0.5 * log((1 + r) / (1 - r))
}

inv_fisher_z <- function(z) {
  (exp(2 * z) - 1) / (exp(2 * z) + 1)
}

is_correlation_metric <- function(vcol) {
  grepl("_r$", vcol) || grepl("pearson_r", vcol, ignore.case = TRUE)
}

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  return("ns")
}

compute_clip_ylim <- function(df, vcol, gcol, k = OUTLIER_CLIP_K) {
  return(NULL)
}

get_sig_pairs <- function(pw_df, antibody_filter = NULL, only_sig = SHOW_ONLY_SIG) {
  if (is.null(pw_df) || nrow(pw_df) == 0)
    return(list(pairs = list(), labels = character()))
  d <- pw_df
  if (!is.null(antibody_filter) && "Antibody" %in% colnames(d))
    d <- d %>% filter(Antibody == antibody_filter)
  if (only_sig) d <- d %>% filter(!is.na(P_adjusted) & P_adjusted < 0.05)
  if (nrow(d) == 0) return(list(pairs = list(), labels = character()))
  
  parse_pair <- function(s) {
    if (grepl(" vs ", s, fixed = TRUE))  parts <- strsplit(s, " vs ", fixed = TRUE)[[1]]
    else if (grepl(" - ", s, fixed = TRUE)) parts <- strsplit(s, " - ", fixed = TRUE)[[1]]
    else return(NULL)
    parts <- gsub("^\\(|\\)$", "", trimws(parts))
    if (length(parts) == 2) parts else NULL
  }
  pairs <- lapply(as.character(d$Contrast), parse_pair)
  keep <- !sapply(pairs, is.null)
  list(pairs = pairs[keep], labels = as.character(d$Significance)[keep])
}

save_plot_pdf_png <- function(plt, filename_base, width = 8.5, height = 11,
                              dpi = 300) {
  pdf_path <- file.path(out_dir_pdf, paste0(filename_base, ".pdf"))
  png_path <- file.path(out_dir_png, paste0(filename_base, ".png"))
  ggsave(pdf_path, plt, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, plt, width = width, height = height, dpi = dpi)
  message("  Saved: ", basename(pdf_path), " + ", basename(png_path))
}

EXPLORATORY_CAVEAT <- paste(
  "EXPLORATORY — n=1 biological replicate.",
  "Treatment confounded with well/batch.",
  "Hypothesis-generating only.")

# ==============================================================================
# LOAD & PREPROCESS DATA
# ==============================================================================

message("Loading data...")
if (!file.exists(master_cells_csv))
  stop("master_cells.csv not found: ", master_cells_csv)
cells_raw <- read_csv(master_cells_csv, show_col_types = FALSE)
message("Loaded ", nrow(cells_raw), " cells from ",
        n_distinct(cells_raw$image_id), " images")

if (!"antibody" %in% colnames(cells_raw) && "combo" %in% colnames(cells_raw)) {
  cells_raw <- cells_raw %>%
    mutate(antibody = case_when(
      grepl("Syn1",      combo, ignore.case = TRUE) ~ "Syn1",
      grepl("EPR|Abcam", combo, ignore.case = TRUE) ~ "EPR20535",
      TRUE ~ NA_character_
    ))
}

HAS_INFECTED <- "infected" %in% colnames(cells_raw)
if (!HAS_INFECTED) {
  warning("Column 'infected' not found — Modules 2/3A/3B/3C/4 will be SKIPPED")
}

cells <- cells_raw %>%
  mutate(
    antibody = case_when(
      grepl("Syn1",              antibody, ignore.case = TRUE) ~ "Syn1",
      grepl("EPR|Abcam|C-?term", antibody, ignore.case = TRUE) ~ "EPR20535",
      TRUE ~ as.character(antibody)
    ),
    condition = factor(condition, levels = group_order),
    Antibody  = factor(antibody, levels = antibody_levels)
  ) %>%
  filter(!is.na(condition), !is.na(Antibody), !is.na(replicate))

if (any(grepl("PFF", as.character(cells_raw$condition)))) {
  cells <- cells %>% filter(!grepl("PFF", as.character(condition)))
}

if (HAS_INFECTED) {
  cells <- cells %>%
    mutate(IAV_state = factor(
      if_else(as.logical(infected), "Infected", "Bystander"),
      levels = c("Bystander", "Infected")))
}

cells$Antibody  <- droplevels(cells$Antibody)
cells$condition <- droplevels(cells$condition)
abs_present     <- levels(cells$Antibody)
conds_present   <- levels(cells$condition)
n_reps_global   <- n_distinct(cells$replicate)

message("Antibody levels present : ", paste(abs_present, collapse = ", "))
message("Condition levels present: ", paste(conds_present, collapse = ", "))
message("Biological replicates   : ", n_reps_global,
        if (n_reps_global == 1) "  <- EXPLORATORY mode (single replicate)" else "")

if (!is.null(exclude_csv) && file.exists(exclude_csv)) {
  excl <- read_csv(exclude_csv, show_col_types = FALSE)
  n_before <- nrow(cells)
  cells <- cells %>% filter(!image_id %in% excl$image_id)
  message("Excluded ", n_before - nrow(cells), " cells from ",
          length(unique(excl$image_id)), " flagged FOVs")
}

if (!is.null(master_coloc_csv) && file.exists(master_coloc_csv)) {
  coloc <- read_csv(master_coloc_csv, show_col_types = FALSE)
  message("Loaded ", nrow(coloc), " coloc rows; merging into cells")
  coloc_cols <- setdiff(colnames(coloc),
                        c(colnames(cells), "nucleus_id", "image_id"))
  cells <- cells %>%
    left_join(coloc %>% select(nucleus_id, image_id, all_of(coloc_cols)),
              by = c("nucleus_id", "image_id"))
}

HAS_COSTES <- "nuc_asyn_iav_costes_p" %in% colnames(cells)
if (!HAS_COSTES) {
  warning("Column 'nuc_asyn_iav_costes_p' not found. ",
          "Module 3C will use |Pearson r| > 0.3 as the colocalisation criterion.")
}

message("\nFinal cells table: ", nrow(cells), " cells x ",
        n_distinct(cells$image_id), " FOVs")
message("Breakdown:")
print(cells %>%
        group_by(Antibody, condition) %>%
        summarize(n_cells = n(),
                  n_FOVs = n_distinct(image_id),
                  .groups = "drop"))

if (HAS_INFECTED) {
  message("\nInfection breakdown (V1, V2 only):")
  print(cells %>%
          filter(condition %in% c("V1", "V2")) %>%
          group_by(Antibody, condition, IAV_state) %>%
          summarize(n_cells = n(),
                    n_FOVs = n_distinct(image_id),
                    .groups = "drop"))
}

# ==============================================================================
# METRIC CONFIGS
# ==============================================================================

metric_configs <- list(
  # ─── Abundance / morphology metrics ───
  # log = TRUE: analysed as log(x + 1e-9)
  NC_ratio = list(col = "NC_ratio", log = TRUE,
                  label = "N/P ratio (aSyn)",
                  type = "abundance",
                  description = "Nuclear/perinuclear aSyn intensity (log)"),
  n_puncta = list(col = "n_nuclear_puncta", log = TRUE,
                  label = "Nuclear puncta per cell",
                  type = "abundance"),
  asyn_nuc = list(col = "asyn_nuclear_mean", log = FALSE,
                  label = "Nuclear aSyn (mean intensity)",
                  type = "abundance",
                  description = "Per-antibody only — NOT cross-antibody comparable"),
  asyn_cyto = list(col = "asyn_cytoplasmic_mean", log = FALSE,
                   label = "Perinuclear aSyn (mean intensity)",
                   type = "abundance",
                   description = "Per-antibody only"),
  # ─── aSyn × DAPI correlation: defined for every nucleus, so run in all abundance modules ───
  coloc_dapi_r = list(col = "nuc_asyn_dapi_pearson_r", log = FALSE,
                      label = "aSyn x DAPI chromatin Pearson r",
                      type = "abundance"),
  # ─── Colocalisation metrics — Module 3C only ───
  coloc_iav_r  = list(col = "nuc_asyn_iav_pearson_r", log = FALSE,
                      label = "aSyn x IAV nuclear Pearson r",
                      type = "coloc"),
  coloc_iav_M1 = list(col = "nuc_asyn_iav_M1", log = FALSE,
                      label = "aSyn x IAV Manders M1",
                      type = "coloc")
  
)
metric_configs <- metric_configs[
  sapply(metric_configs, function(m) m$col %in% colnames(cells))
]
message("\nMetrics found in data: ",
        paste(names(metric_configs), collapse = ", "))

abundance_metrics <- names(Filter(function(m) m$type == "abundance",
                                  metric_configs))
coloc_metrics     <- names(Filter(function(m) m$type == "coloc",
                                  metric_configs))

# ==============================================================================
# CORE FITTER — LMER ADAPTIVE WITH FISHER Z FOR CORRELATIONS
# ==============================================================================

fit_lmer_adaptive <- function(df, vcol, fixed_effects_str,
                              log_transform = FALSE) {
  d <- df %>% filter(!is.na(.data[[vcol]]), is.finite(.data[[vcol]]))
  if (nrow(d) < MIN_CELLS_FOR_LMER) return(NULL)
  
  is_corr <- is_correlation_metric(vcol)
  if (is_corr) {
    # Fisher z instead of log
    d$.y <- fisher_z(d[[vcol]])
  } else {
    d$.y <- if (log_transform) log(d[[vcol]] + 1e-9) else d[[vcol]]
  }
  
  n_fovs <- n_distinct(d$image_id)
  f_lmer <- as.formula(paste(".y ~", fixed_effects_str, "+ (1|image_id)"))
  f_lm   <- as.formula(paste(".y ~", fixed_effects_str))
  
  fit <- if (n_fovs >= 2) {
    suppressWarnings(suppressMessages(
      tryCatch(lmerTest::lmer(f_lmer, data = d, REML = TRUE),
               error = function(e) NULL)
    ))
  } else NULL
  
  used_lm <- FALSE
  if (is.null(fit)) {
    fit <- suppressWarnings(suppressMessages(
      tryCatch(lm(f_lm, data = d), error = function(e) NULL)
    ))
    used_lm <- !is.null(fit)
  }
  if (is.null(fit)) return(NULL)
  list(fit = fit, data = d, used_lm = used_lm, is_correlation = is_corr)
}

# Back-transform Fisher-z emmeans output to the r scale
back_transform_emm <- function(emm_df, is_corr) {
  if (!is_corr) return(emm_df)
  for (col in c("emmean", "lower.CL", "upper.CL")) {
    if (col %in% colnames(emm_df)) emm_df[[col]] <- inv_fisher_z(emm_df[[col]])
  }
  emm_df
}

# Tidy ANOVA table
clean_aov <- function(fit, metric_name, module_label) {
  aov_tab <- tryCatch({
    if (inherits(fit, "merMod")) {
      as.data.frame(anova(fit, type = 3, ddf = "Satterthwaite"))
    } else {
      opts_old <- options(contrasts = c("contr.sum", "contr.poly"))
      a <- as.data.frame(car::Anova(fit, type = 3))
      options(opts_old)
      a <- a[!rownames(a) %in% c("(Intercept)", "Residuals"), , drop = FALSE]
      a
    }
  }, error = function(e) NULL)
  if (is.null(aov_tab) || nrow(aov_tab) == 0) return(NULL)
  aov_tab$Term         <- rownames(aov_tab)
  aov_tab$Metric       <- metric_name
  aov_tab$Module       <- module_label
  p_col <- if ("Pr(>F)" %in% colnames(aov_tab)) "Pr(>F)" else "Pr(>Chisq)"
  aov_tab$Significance <- sapply(aov_tab[[p_col]], p_to_stars)
  aov_tab
}

# Helper function to add model transparency columns to data frames
add_transparency_cols <- function(df, model_obj, include_cells = FALSE) {
  if (is.null(df) || is.null(model_obj)) return(df)
  model_type <- if (model_obj$used_lm) "LM" else "LMER"
  n_fovs <- n_distinct(model_obj$data$image_id)
  
  df$Model_Type <- model_type
  df$n_FOVs     <- n_fovs
  if (include_cells) {
    df$n_cells_used <- nrow(model_obj$data)
  }
  df
}

# Stratified/Standard Cluster Bootstrap with Cell Fallback
cluster_bootstrap_contrast <- function(df, vcol, fixed_effects_str, 
                                      contrast_fn, log_transform = FALSE,
                                      n_boot = N_BOOTSTRAP, stratify_col = NULL) {
  # Pre-filter NA values
  d_clean <- df %>% filter(!is.na(.data[[vcol]]), is.finite(.data[[vcol]]))
  if (nrow(d_clean) < 2) return(NULL)
  
  # Determine overall boot mode
  is_corr <- is_correlation_metric(vcol)
  
  # Check unique FOVs per stratum
  if (!is.null(stratify_col)) {
    strata_levels <- unique(d_clean[[stratify_col]])
    fovs_per_stratum <- sapply(strata_levels, function(lvl) {
      n_distinct(d_clean$image_id[d_clean[[stratify_col]] == lvl])
    })
    any_single_fov <- any(fovs_per_stratum == 1)
    boot_mode <- if (any_single_fov) "cell" else "cluster"
  } else {
    n_fovs_total <- n_distinct(d_clean$image_id)
    boot_mode <- if (n_fovs_total <= 1) "cell" else "cluster"
  }
  
  # We will collect bootstrap estimates
  boot_list <- vector("list", n_boot)
  
  message("  Running bootstrap (n=", n_boot, ", mode=", boot_mode, ")...")
  
  for (i in 1:n_boot) {
    if (i %% 1000 == 0) {
      message("    Completed ", i, " / ", n_boot, " iterations...")
    }
    
    boot_df <- tryCatch({
      if (boot_mode == "cell") {
        # Cell-level bootstrap: sample rows with replacement
        if (!is.null(stratify_col)) {
          # Stratified cell-level bootstrap
          d_clean %>%
            group_by(.data[[stratify_col]]) %>%
            slice_sample(prop = 1, replace = TRUE) %>%
            ungroup()
        } else {
          d_clean %>% slice_sample(prop = 1, replace = TRUE)
        }
      } else {
        # Cluster-level bootstrap
        if (!is.null(stratify_col)) {
          # Stratified cluster bootstrap
          strata_levels <- unique(d_clean[[stratify_col]])
          boot_parts <- lapply(seq_along(strata_levels), function(idx) {
            lvl <- strata_levels[idx]
            sub_df <- d_clean %>% filter(.data[[stratify_col]] == lvl)
            st_image_ids <- unique(sub_df$image_id)
            sampled_ids <- sample(st_image_ids, size = length(st_image_ids), replace = TRUE)
            
            lapply(seq_along(sampled_ids), function(j) {
              id <- sampled_ids[j]
              clus_df <- sub_df %>% filter(image_id == id)
              clus_df$image_id <- paste0(id, "_boot_", idx, "_", j)
              clus_df
            }) %>% bind_rows()
          })
          bind_rows(boot_parts)
        } else {
          # Standard cluster bootstrap
          image_ids <- unique(d_clean$image_id)
          sampled_ids <- sample(image_ids, size = length(image_ids), replace = TRUE)
          lapply(seq_along(sampled_ids), function(j) {
            id <- sampled_ids[j]
            clus_df <- d_clean %>% filter(image_id == id)
            clus_df$image_id <- paste0(id, "_boot_", j)
            clus_df
          }) %>% bind_rows()
        }
      }
    }, error = function(e) NULL)
    
    if (is.null(boot_df) || nrow(boot_df) == 0) {
      boot_list[[i]] <- NULL
      next
    }
    
    # Fit model on bootstrap sample
    model_obj <- tryCatch({
      fit_lmer_adaptive(boot_df, vcol, fixed_effects_str, log_transform)
    }, error = function(e) NULL)
    
    if (is.null(model_obj) || is.null(model_obj$fit)) {
      boot_list[[i]] <- NULL
      next
    }
    
    # Extract statistic
    est <- suppressWarnings(suppressMessages(
      tryCatch({
        contrast_fn(model_obj$fit, boot_df)
      }, error = function(e) NULL)
    ))
    
    boot_list[[i]] <- est
  }
  
  # Remove null/failed iterations
  valid_idx <- sapply(boot_list, function(x) !is.null(x) && !any(is.na(x)))
  n_successful <- sum(valid_idx)
  
  if (n_successful == 0) {
    return(NULL)
  }
  
  # Extract names if they are vector statistics
  first_valid <- boot_list[valid_idx][[1]]
  K <- length(first_valid)
  stat_names <- names(first_valid)
  
  # Compile estimates
  estimates_mat <- matrix(NA_real_, nrow = n_boot, ncol = K)
  for (i in which(valid_idx)) {
    estimates_mat[i, ] <- boot_list[[i]]
  }
  
  # Back-transform if correlation
  if (is_corr) {
    # Apply inv_fisher_z to estimates
    estimates_mat <- inv_fisher_z(estimates_mat)
  }
  
  results <- list()
  for (k in 1:K) {
    k_est <- estimates_mat[, k]
    k_valid <- k_est[!is.na(k_est)]
    
    if (length(k_valid) == 0) {
      results[[k]] <- list(
        mean = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        p_boot = NA_real_, n_successful = 0, boot_mode = boot_mode
      )
      next
    }
    
    boot_mean   <- mean(k_valid)
    boot_median <- median(k_valid)
    ci_low      <- quantile(k_valid, 0.025, na.rm = TRUE)
    ci_high     <- quantile(k_valid, 0.975, na.rm = TRUE)
    
    # P_boot: 2 * min(mean(boot <= 0), mean(boot >= 0))
    p_val <- 2 * min(mean(k_valid <= 0, na.rm = TRUE), mean(k_valid >= 0, na.rm = TRUE))
    p_val <- min(1, p_val)
    
    results[[k]] <- list(
      mean = boot_mean,
      median = boot_median,
      CI_low = as.numeric(ci_low),
      CI_high = as.numeric(ci_high),
      p_boot = p_val,
      n_successful = n_successful,
      boot_mode = boot_mode
    )
  }
  
  if (!is.null(stat_names)) {
    names(results) <- stat_names
  }
  results
}

# ==============================================================================
# MODULE 1 — TREATMENT LEVEL [EXPLORATORY]
# ==============================================================================

run_module1 <- function(cells_df, vcol, log_transform = FALSE,
                        metric_name = "metric") {
  n_ab   <- n_distinct(cells_df$Antibody)
  n_cond <- n_distinct(cells_df$condition)
  if (n_cond < 2) return(NULL)
  
  fixed_str <- if (n_ab >= 2) "condition * Antibody" else "condition"
  
  model_obj <- fit_lmer_adaptive(cells_df, vcol, fixed_str, log_transform)
  if (is.null(model_obj)) return(NULL)
  fit <- model_obj$fit; d <- model_obj$data; is_corr <- model_obj$is_correlation
  
  aov_tab <- add_transparency_cols(clean_aov(fit, metric_name, "1_Treatment_EXPLORATORY"), model_obj)
  
  emm_spec <- if (n_ab >= 2) ~ condition * Antibody else ~ condition
  emm_obj  <- tryCatch(emmeans(fit, emm_spec), error = function(e) NULL)
  emm_df <- if (!is.null(emm_obj)) {
    add_transparency_cols(
      back_transform_emm(as.data.frame(emm_obj), is_corr) %>%
        mutate(Metric = metric_name, Module = "1_Treatment_EXPLORATORY"),
      model_obj
    )
  } else NULL
  
  # Planned contrasts
  pw_df <- NULL
  emm_cond <- if (n_ab >= 2) {
    tryCatch(emmeans(fit, ~ condition | Antibody), error = function(e) NULL)
  } else {
    tryCatch(emmeans(fit, ~ condition), error = function(e) NULL)
  }
  if (!is.null(emm_cond)) {
    present <- levels(d$condition)
    valid_comps <- Filter(function(p) all(p %in% present), my_comparisons_M1)
    if (length(valid_comps) > 0) {
      contrast_list <- lapply(valid_comps, function(pair) {
        vec <- rep(0, length(present))
        vec[present == pair[1]] <- -1
        vec[present == pair[2]] <-  1
        vec
      })
      names(contrast_list) <- sapply(valid_comps,
                                     function(p) paste(p[2], "vs", p[1]))
      ctr_raw <- tryCatch(as.data.frame(summary(
        contrast(emm_cond, method = contrast_list, adjust = "none"))),
        error = function(e) NULL)
      ctr_adj <- tryCatch(as.data.frame(summary(
        contrast(emm_cond, method = contrast_list,
                 adjust = ADJUST_METHOD))),
        error = function(e) NULL)
      if (!is.null(ctr_raw) && !is.null(ctr_adj)) {
        ab_col <- if ("Antibody" %in% colnames(ctr_raw))
          as.character(ctr_raw$Antibody) else rep("ALL", nrow(ctr_raw))
        pw_df <- data.frame(
          Metric       = metric_name,
          Module       = "1_Treatment_EXPLORATORY",
          Antibody     = ab_col,
          Contrast     = as.character(ctr_raw$contrast),
          Estimate     = ctr_raw$estimate,
          SE           = ctr_raw$SE,
          t_ratio      = ctr_raw$t.ratio,
          df           = ctr_raw$df,
          P_raw        = ctr_raw$p.value,
          P_adjusted   = ctr_adj$p.value,
          Significance = sapply(ctr_adj$p.value, p_to_stars),
          stringsAsFactors = FALSE
        )
        pw_df <- add_transparency_cols(pw_df, model_obj)
      }
    }
  }
  
  list(anova = aov_tab, marginals = emm_df, pairwise = pw_df,
       n_cells = nrow(d), n_fovs = n_distinct(d$image_id),
       used_lm = model_obj$used_lm, is_correlation = is_corr)
}

# ==============================================================================
# MODULE 2 — WITHIN-FOV INFECTION CONTRAST [PRIMARY - MAIN]
# ==============================================================================
# Run separately for V1 and V2 (Control has no infected cells).
# Not run for aSyn × IAV colocalisation metrics, which are uninformative in
# bystander cells (no IAV signal).
# ==============================================================================

run_module2 <- function(cells_df, vcol, log_transform = FALSE,
                        metric_name = "metric", virus_label = "V1") {
  sub_df <- cells_df %>% filter(condition == virus_label,
                                !is.na(.data[[vcol]]), !is.na(IAV_state))
  if (nrow(sub_df) < MIN_CELLS_FOR_LMER) return(NULL)
  
  state_counts <- sub_df %>% count(IAV_state)
  if (nrow(state_counts) < 2 || any(state_counts$n < MIN_CELLS_PER_GROUP))
    return(NULL)
  
  n_ab <- n_distinct(sub_df$Antibody)
  fixed_str <- if (n_ab >= 2) "IAV_state * Antibody" else "IAV_state"
  
  model_obj <- fit_lmer_adaptive(sub_df, vcol, fixed_str, log_transform)
  if (is.null(model_obj)) return(NULL)
  fit <- model_obj$fit; d <- model_obj$data; is_corr <- model_obj$is_correlation
  
  aov_tab <- add_transparency_cols(clean_aov(fit, metric_name, "2_WithinFOV_PRIMARY"), model_obj)
  if (!is.null(aov_tab)) aov_tab$Virus <- virus_label
  
  emm_spec <- if (n_ab >= 2) ~ IAV_state * Antibody else ~ IAV_state
  emm_obj  <- tryCatch(emmeans(fit, emm_spec), error = function(e) NULL)
  emm_df <- if (!is.null(emm_obj)) {
    add_transparency_cols(
      back_transform_emm(as.data.frame(emm_obj), is_corr) %>%
        mutate(Metric = metric_name, Virus = virus_label,
               Module = "2_WithinFOV_PRIMARY"),
      model_obj
    )
  } else NULL
  
  emm_pair <- if (n_ab >= 2) {
    tryCatch(emmeans(fit, ~ IAV_state | Antibody), error = function(e) NULL)
  } else {
    tryCatch(emmeans(fit, ~ IAV_state), error = function(e) NULL)
  }
  pw_df <- NULL
  if (!is.null(emm_pair)) {
    ctr_raw <- tryCatch(as.data.frame(summary(
      contrast(emm_pair, method = "pairwise", adjust = "none"))),
      error = function(e) NULL)
    ctr_adj <- tryCatch(as.data.frame(summary(
      contrast(emm_pair, method = "pairwise", adjust = ADJUST_METHOD))),
      error = function(e) NULL)
    if (!is.null(ctr_raw) && !is.null(ctr_adj)) {
      ab_col <- if ("Antibody" %in% colnames(ctr_raw))
        as.character(ctr_raw$Antibody) else rep("ALL", nrow(ctr_raw))
      pw_df <- data.frame(
        Metric       = metric_name,
        Module       = "2_WithinFOV_PRIMARY",
        Virus        = virus_label,
        Antibody     = ab_col,
        Contrast     = as.character(ctr_raw$contrast),
        Estimate     = if (is_corr) inv_fisher_z(ctr_raw$estimate) else ctr_raw$estimate,
        SE           = ctr_raw$SE,
        t_ratio      = ctr_raw$t.ratio,
        df           = ctr_raw$df,
        P_raw        = ctr_raw$p.value,
        P_adjusted   = ctr_adj$p.value,
        Significance = sapply(ctr_adj$p.value, p_to_stars),
        stringsAsFactors = FALSE
      )
      pw_df <- add_transparency_cols(pw_df, model_obj)
    }
  }
  
  # Bootstrap supplement
  if ("M2" %in% BOOTSTRAP_MODULES) {
    contrast_fn_m2 <- function(model_fit, model_data) {
      n_ab_fit <- n_distinct(model_data$Antibody)
      emm_pair_fit <- if (n_ab_fit >= 2) {
        emmeans::emmeans(model_fit, ~ IAV_state | Antibody)
      } else {
        emmeans::emmeans(model_fit, ~ IAV_state)
      }
      ctr <- tryCatch(
        as.data.frame(summary(contrast(emm_pair_fit, method = "pairwise", adjust = "none"))),
        error = function(e) NULL
      )
      if (is.null(ctr)) return(NULL)
      if ("Antibody" %in% colnames(ctr)) {
        est_vec <- ctr$estimate
        names(est_vec) <- as.character(ctr$Antibody)
        est_vec
      } else {
        c("ALL" = ctr$estimate[1])
      }
    }
    
    boot_res <- cluster_bootstrap_contrast(
      df = sub_df, vcol = vcol, fixed_effects_str = fixed_str,
      contrast_fn = contrast_fn_m2, log_transform = log_transform,
      n_boot = N_BOOTSTRAP
    )
    
    if (!is.null(boot_res) && !is.null(pw_df) && nrow(pw_df) > 0) {
      pw_df$Boot_Estimate     <- NA_real_
      pw_df$Boot_CI_Low       <- NA_real_
      pw_df$Boot_CI_High      <- NA_real_
      pw_df$P_boot            <- NA_real_
      pw_df$Boot_Mode         <- NA_character_
      pw_df$Boot_N_Successful <- NA_integer_
      
      for (i in 1:nrow(pw_df)) {
        ab <- pw_df$Antibody[i]
        key <- if (ab %in% names(boot_res)) ab else "ALL"
        if (key %in% names(boot_res)) {
          pw_df$Boot_Estimate[i]     <- boot_res[[key]]$mean
          pw_df$Boot_CI_Low[i]       <- boot_res[[key]]$CI_low
          pw_df$Boot_CI_High[i]      <- boot_res[[key]]$CI_high
          pw_df$P_boot[i]            <- boot_res[[key]]$p_boot
          pw_df$Boot_Mode[i]         <- boot_res[[key]]$boot_mode
          pw_df$Boot_N_Successful[i] <- boot_res[[key]]$n_successful
        }
      }
    }
  }
  
  list(anova = aov_tab, marginals = emm_df, pairwise = pw_df,
       n_cells = nrow(d), n_fovs = n_distinct(d$image_id),
       n_infected = sum(d$IAV_state == "Infected"),
       n_bystander = sum(d$IAV_state == "Bystander"),
       virus = virus_label, used_lm = model_obj$used_lm,
       is_correlation = is_corr)
}

# ==============================================================================
# MODULE 3A — PARACRINE EFFECT 
# ==============================================================================

run_module3a_paracrine <- function(cells_df, vcol, log_transform = FALSE,
                                   metric_name = "metric") {
  para_df <- cells_df %>%
    mutate(para_group = case_when(
      condition == "Control" ~ "Control",
      condition == "V1" & IAV_state == "Bystander" ~ "V1.Bystander",
      condition == "V2" & IAV_state == "Bystander" ~ "V2.Bystander",
      TRUE ~ NA_character_
    )) %>% filter(!is.na(para_group),
                  !is.na(.data[[vcol]]),
                  is.finite(.data[[vcol]]))
  
  para_df$para_group <- factor(para_df$para_group,
                               levels = c("Control", "V1.Bystander",
                                          "V2.Bystander"))
  para_df$para_group <- droplevels(para_df$para_group)
  if (nlevels(para_df$para_group) < 2) return(NULL)
  if (nrow(para_df) < MIN_CELLS_FOR_LMER) return(NULL)
  
  n_ab <- n_distinct(para_df$Antibody)
  fixed_str <- if (n_ab >= 2) "para_group * Antibody" else "para_group"
  
  model_obj <- fit_lmer_adaptive(para_df, vcol, fixed_str, log_transform)
  if (is.null(model_obj)) return(NULL)
  fit <- model_obj$fit; d <- model_obj$data; is_corr <- model_obj$is_correlation
  
  aov_tab <- add_transparency_cols(clean_aov(fit, metric_name, "3A_Paracrine"), model_obj)
  
  emm_spec <- if (n_ab >= 2) ~ para_group * Antibody else ~ para_group
  emm_obj  <- tryCatch(emmeans(fit, emm_spec), error = function(e) NULL)
  emm_df <- if (!is.null(emm_obj)) {
    add_transparency_cols(
      back_transform_emm(as.data.frame(emm_obj), is_corr) %>%
        mutate(Metric = metric_name, Module = "3A_Paracrine"),
      model_obj
    )
  } else NULL
  
  emm_pair <- if (n_ab >= 2) {
    tryCatch(emmeans(fit, ~ para_group | Antibody), error = function(e) NULL)
  } else {
    tryCatch(emmeans(fit, ~ para_group), error = function(e) NULL)
  }
  pw_df <- NULL
  if (!is.null(emm_pair)) {
    ctr_raw <- tryCatch(as.data.frame(summary(
      contrast(emm_pair, method = "pairwise", adjust = "none"))),
      error = function(e) NULL)
    ctr_adj <- tryCatch(as.data.frame(summary(
      contrast(emm_pair, method = "pairwise", adjust = ADJUST_METHOD))),
      error = function(e) NULL)
    if (!is.null(ctr_raw) && !is.null(ctr_adj)) {
      ab_col <- if ("Antibody" %in% colnames(ctr_raw))
        as.character(ctr_raw$Antibody) else rep("ALL", nrow(ctr_raw))
      pw_df <- data.frame(
        Metric       = metric_name,
        Module       = "3A_Paracrine",
        Antibody     = ab_col,
        Contrast     = as.character(ctr_raw$contrast),
        Estimate     = ctr_raw$estimate,
        SE           = ctr_raw$SE,
        t_ratio      = ctr_raw$t.ratio,
        df           = ctr_raw$df,
        P_raw        = ctr_raw$p.value,
        P_adjusted   = ctr_adj$p.value,
        Significance = sapply(ctr_adj$p.value, p_to_stars),
        stringsAsFactors = FALSE
      )
      pw_df <- add_transparency_cols(pw_df, model_obj)
    }
  }
  
  list(anova = aov_tab, marginals = emm_df, pairwise = pw_df,
       n_cells = nrow(d), n_fovs = n_distinct(d$image_id),
       used_lm = model_obj$used_lm, is_correlation = is_corr)
}

# ==============================================================================
# MODULE 3B — STRAIN INTERACTION (V1 vs V2, condition * IAV_state)
# ==============================================================================
# V1 and V2 only (Control has no infected cells).
# Tests whether the infected-vs-bystander difference depends on the strain.
# ==============================================================================

run_module3b_strain <- function(cells_df, vcol, log_transform = FALSE,
                                metric_name = "metric") {
  d_filt <- cells_df %>%
    filter(condition %in% c("V1", "V2"),
           !is.na(.data[[vcol]]), is.finite(.data[[vcol]]),
           !is.na(IAV_state))
  d_filt$condition <- droplevels(d_filt$condition)
  d_filt$IAV_state <- droplevels(d_filt$IAV_state)
  
  if (n_distinct(d_filt$condition) < 2) return(NULL)
  if (n_distinct(d_filt$IAV_state) < 2) return(NULL)
  if (nrow(d_filt) < MIN_CELLS_FOR_LMER) return(NULL)
  
  n_ab <- n_distinct(d_filt$Antibody)
  fixed_str <- if (n_ab >= 2)
    "condition * IAV_state + Antibody"
  else
    "condition * IAV_state"
  
  model_obj <- fit_lmer_adaptive(d_filt, vcol, fixed_str, log_transform)
  if (is.null(model_obj)) return(NULL)
  fit <- model_obj$fit; d <- model_obj$data; is_corr <- model_obj$is_correlation
  
  aov_tab <- add_transparency_cols(clean_aov(fit, metric_name, "3B_StrainInteraction"), model_obj)
  
  emm_spec <- if (n_ab >= 2) ~ condition * IAV_state * Antibody
  else            ~ condition * IAV_state
  emm_obj  <- tryCatch(emmeans(fit, emm_spec), error = function(e) NULL)
  emm_df <- if (!is.null(emm_obj)) {
    add_transparency_cols(
      back_transform_emm(as.data.frame(emm_obj), is_corr) %>%
        mutate(Metric = metric_name, Module = "3B_StrainInteraction"),
      model_obj
    )
  } else NULL
  
  # Targeted pairwise: V1.Inf vs V2.Inf, V1.Bys vs V2.Bys, Inf vs Bys within each strain
  d$cond_state <- factor(paste(d$condition, d$IAV_state, sep = "."))
  pw_df <- NULL
  fit_cs <- tryCatch({
    if (inherits(fit, "merMod")) {
      f_cs <- if (n_ab >= 2)
        as.formula(".y ~ cond_state + Antibody + (1|image_id)")
      else
        as.formula(".y ~ cond_state + (1|image_id)")
      lmerTest::lmer(f_cs, data = d, REML = TRUE)
    } else {
      f_cs <- if (n_ab >= 2)
        as.formula(".y ~ cond_state + Antibody")
      else
        as.formula(".y ~ cond_state")
      lm(f_cs, data = d)
    }
  }, error = function(e) NULL)
  
  if (!is.null(fit_cs)) {
    emm_cs <- tryCatch(emmeans(fit_cs, ~ cond_state), error = function(e) NULL)
    if (!is.null(emm_cs)) {
      present_cs <- levels(d$cond_state)
      wanted <- list(
        c("V1.Bystander", "V2.Bystander"),
        c("V1.Infected",  "V2.Infected"),
        c("V1.Infected",  "V1.Bystander"),
        c("V2.Infected",  "V2.Bystander")
      )
      valid <- Filter(function(p) all(p %in% present_cs), wanted)
      if (length(valid) > 0) {
        contrast_list <- lapply(valid, function(pair) {
          vec <- rep(0, length(present_cs))
          vec[present_cs == pair[1]] <-  1
          vec[present_cs == pair[2]] <- -1
          vec
        })
        names(contrast_list) <- sapply(valid,
                                       function(p) paste(p[1], "vs", p[2]))
        ctr_raw <- tryCatch(as.data.frame(summary(
          contrast(emm_cs, method = contrast_list, adjust = "none"))),
          error = function(e) NULL)
        ctr_adj <- tryCatch(as.data.frame(summary(
          contrast(emm_cs, method = contrast_list,
                   adjust = ADJUST_METHOD))),
          error = function(e) NULL)
        if (!is.null(ctr_raw) && !is.null(ctr_adj)) {
          pw_df <- data.frame(
            Metric       = metric_name,
            Module       = "3B_StrainInteraction",
            Contrast     = as.character(ctr_raw$contrast),
            Estimate     = if (is_corr) inv_fisher_z(ctr_raw$estimate) else ctr_raw$estimate,
            SE           = ctr_raw$SE,
            t_ratio      = ctr_raw$t.ratio,
            df           = ctr_raw$df,
            P_raw        = ctr_raw$p.value,
            P_adjusted   = ctr_adj$p.value,
            Significance = sapply(ctr_adj$p.value, p_to_stars),
            stringsAsFactors = FALSE
          )
          pw_df <- add_transparency_cols(pw_df, model_obj)
        }
      }
    }
  }
  
  if ("M3B" %in% BOOTSTRAP_MODULES) {
    contrast_fn_m3b <- function(model_fit, model_data) {
      emm_cs_fit <- tryCatch(
        emmeans::emmeans(model_fit, ~ cond_state),
        error = function(e) NULL
      )
      if (is.null(emm_cs_fit)) return(c(
        "V1.Bystander vs V2.Bystander" = NA_real_,
        "V1.Infected vs V2.Infected"   = NA_real_,
        "V1.Infected vs V1.Bystander"  = NA_real_,
        "V2.Infected vs V2.Bystander"  = NA_real_
      ))
      
      present_cs <- levels(model_data$cond_state)
      wanted <- list(
        "V1.Bystander vs V2.Bystander" = c("V1.Bystander", "V2.Bystander"),
        "V1.Infected vs V2.Infected"   = c("V1.Infected",  "V2.Infected"),
        "V1.Infected vs V1.Bystander"  = c("V1.Infected",  "V1.Bystander"),
        "V2.Infected vs V2.Bystander"  = c("V2.Infected",  "V2.Bystander")
      )
      
      est_vec <- numeric(0)
      for (c_name in names(wanted)) {
        pair <- wanted[[c_name]]
        if (all(pair %in% present_cs)) {
          vec <- rep(0, length(present_cs))
          vec[present_cs == pair[1]] <-  1
          vec[present_cs == pair[2]] <- -1
          ctr_summary <- tryCatch(
            as.data.frame(summary(contrast(emm_cs_fit, method = list(vec), adjust = "none"))),
            error = function(e) NULL
          )
          est_vec[c_name] <- if (!is.null(ctr_summary)) ctr_summary$estimate[1] else NA_real_
        } else {
          est_vec[c_name] <- NA_real_
        }
      }
      est_vec
    }
    
    fixed_str_cs <- if (n_ab >= 2) "cond_state + Antibody" else "cond_state"
    boot_res <- cluster_bootstrap_contrast(
      df = d, vcol = vcol, fixed_effects_str = fixed_str_cs,
      contrast_fn = contrast_fn_m3b, log_transform = log_transform,
      n_boot = N_BOOTSTRAP, stratify_col = "condition"
    )
    
    if (!is.null(boot_res) && !is.null(pw_df) && nrow(pw_df) > 0) {
      pw_df$Boot_Estimate     <- NA_real_
      pw_df$Boot_CI_Low       <- NA_real_
      pw_df$Boot_CI_High      <- NA_real_
      pw_df$P_boot            <- NA_real_
      pw_df$Boot_Mode         <- NA_character_
      pw_df$Boot_N_Successful <- NA_integer_
      
      for (i in 1:nrow(pw_df)) {
        contrast_name <- pw_df$Contrast[i]
        if (contrast_name %in% names(boot_res)) {
          pw_df$Boot_Estimate[i]     <- boot_res[[contrast_name]]$mean
          pw_df$Boot_CI_Low[i]       <- boot_res[[contrast_name]]$CI_low
          pw_df$Boot_CI_High[i]      <- boot_res[[contrast_name]]$CI_high
          pw_df$P_boot[i]            <- boot_res[[contrast_name]]$p_boot
          pw_df$Boot_Mode[i]         <- boot_res[[contrast_name]]$boot_mode
          pw_df$Boot_N_Successful[i] <- boot_res[[contrast_name]]$n_successful
        }
      }
    }
  }
  
  list(anova = aov_tab, marginals = emm_df, pairwise = pw_df,
       n_cells = nrow(d), n_fovs = n_distinct(d$image_id),
       used_lm = model_obj$used_lm, is_correlation = is_corr)
}

# ==============================================================================
# MODULE 3C — COLOCALIZATION PROFILING
# ==============================================================================

# "true_coloc" flag: Costes p < threshold (fallback |r| > 0.3 if Costes is unavailable)
make_true_coloc_flag <- function(cells_df) {
  if (HAS_COSTES) {
    cells_df %>%
      mutate(true_coloc = !is.na(nuc_asyn_iav_costes_p) &
               nuc_asyn_iav_costes_p < COSTES_P_THRESHOLD)
  } else {
    cells_df %>%
      mutate(true_coloc = !is.na(nuc_asyn_iav_pearson_r) &
               abs(nuc_asyn_iav_pearson_r) > 0.3)
  }
}

# ─── 3C.a — Frequency test of "true colocalization" ───
run_module3c_frequency <- function(cells_df) {
  d <- cells_df %>%
    filter(condition %in% c("V1", "V2"), IAV_state == "Infected") %>%
    make_true_coloc_flag()
  
  if (nrow(d) < 10) return(NULL)
  
  # Frequency table
  freq_summary <- d %>%
    group_by(condition, Antibody) %>%
    summarize(
      Total_Infected_Cells = n(),
      True_Coloc_Cells     = sum(true_coloc, na.rm = TRUE),
      Pct_True_Coloc       = 100 * True_Coloc_Cells / Total_Infected_Cells,
      .groups = "drop"
    )
  
  res_list <- list()
  
  # Chi-square / Fisher: V1 vs V2 (within each antibody)
  for (ab in levels(d$Antibody)) {
    sub <- d %>% filter(Antibody == ab)
    if (n_distinct(sub$condition) < 2) next
    tab <- table(sub$condition, sub$true_coloc)
    if (any(dim(tab) < 2)) next
    test <- tryCatch({
      if (min(tab) < 5) fisher.test(tab) else chisq.test(tab)
    }, error = function(e) NULL)
    
    if (!is.null(test)) {
      res_list[[length(res_list) + 1]] <- tibble(
        Comparison   = "V1 vs V2",
        Antibody     = ab,
        Method       = as.character(test$method[1]),
        Statistic    = if (!is.null(test$statistic)) as.numeric(test$statistic[1]) else NA_real_,
        df           = if (!is.null(test$parameter)) as.numeric(test$parameter[1]) else NA_real_,
        P_raw        = as.numeric(test$p.value[1]),
        Significance = p_to_stars(test$p.value[1])
      )
    }
  }
  
  for (cnd in c("V1", "V2")) {
    sub <- d %>% filter(condition == cnd)
    if (n_distinct(sub$Antibody) < 2) next
    tab <- table(sub$Antibody, sub$true_coloc)
    if (any(dim(tab) < 2)) next
    test <- tryCatch({
      if (min(tab) < 5) fisher.test(tab) else chisq.test(tab)
    }, error = function(e) NULL)
    
    if (!is.null(test)) {
      res_list[[length(res_list) + 1]] <- tibble(
        Comparison   = "Syn1 vs EPR20535",
        Antibody     = cnd,
        Method       = as.character(test$method[1]),
        Statistic    = if (!is.null(test$statistic)) as.numeric(test$statistic[1]) else NA_real_,
        df           = if (!is.null(test$parameter)) as.numeric(test$parameter[1]) else NA_real_,
        P_raw        = as.numeric(test$p.value[1]),
        Significance = p_to_stars(test$p.value[1])
      )
    }
  }
  
  tests_df <- NULL
  if (length(res_list) > 0) {
    tests_df <- bind_rows(res_list)
    # Adjust only if P_raw is present
    if ("P_raw" %in% names(tests_df) && nrow(tests_df) > 0) {
      tests_df$P_adjusted <- p.adjust(tests_df$P_raw, method = ADJUST_METHOD)
      tests_df$Significance <- sapply(tests_df$P_adjusted, p_to_stars)
      colnames(tests_df)[colnames(tests_df) == "Antibody"] <- "Stratum"
    }
  }
  
  list(frequency = freq_summary, tests = tests_df)
}

# ─── 3C.b — LME on cells passing Costes filter ───
run_module3c_coloc_lme <- function(cells_df, vcol, metric_name) {
  filt_df <- cells_df %>%
    filter(condition %in% c("V1", "V2"), IAV_state == "Infected") %>%
    make_true_coloc_flag() %>%
    filter(true_coloc, !is.na(.data[[vcol]]), is.finite(.data[[vcol]]))
  filt_df$condition <- droplevels(filt_df$condition)
  
  if (nrow(filt_df) < MIN_CELLS_FOR_COLOC) return(NULL)
  if (n_distinct(filt_df$condition) < 2 &&
      n_distinct(filt_df$Antibody)  < 2) return(NULL)
  
  n_ab <- n_distinct(filt_df$Antibody)
  n_cd <- n_distinct(filt_df$condition)
  
  fixed_str <- if (n_ab >= 2 && n_cd >= 2) "condition * Antibody"
  else if (n_cd >= 2)         "condition"
  else                        "Antibody"
  
  model_obj <- fit_lmer_adaptive(filt_df, vcol, fixed_str,
                                 log_transform = FALSE)
  if (is.null(model_obj)) return(NULL)
  fit <- model_obj$fit; d <- model_obj$data; is_corr <- model_obj$is_correlation
  
  aov_tab <- add_transparency_cols(clean_aov(fit, metric_name, "3C_ColocLME_Filtered"), model_obj)
  
  emm_spec <- if (n_ab >= 2 && n_cd >= 2) ~ condition * Antibody
  else if (n_cd >= 2)         ~ condition
  else                        ~ Antibody
  emm_obj <- tryCatch(emmeans(fit, emm_spec), error = function(e) NULL)
  emm_df <- if (!is.null(emm_obj)) {
    add_transparency_cols(
      back_transform_emm(as.data.frame(emm_obj), is_corr) %>%
        mutate(Metric = metric_name, Module = "3C_ColocLME_Filtered"),
      model_obj
    )
  } else NULL
  
  pw_list <- list()
  
  if (n_cd >= 2) {
    emm_cond <- if (n_ab >= 2) {
      tryCatch(emmeans(fit, ~ condition | Antibody), error = function(e) NULL)
    } else {
      tryCatch(emmeans(fit, ~ condition), error = function(e) NULL)
    }
    if (!is.null(emm_cond)) {
      ctr_raw <- tryCatch(as.data.frame(summary(contrast(emm_cond,
                                                         method = "pairwise", adjust = "none"))), error = function(e) NULL)
      ctr_adj <- tryCatch(as.data.frame(summary(contrast(emm_cond,
                                                         method = "pairwise", adjust = ADJUST_METHOD))), error = function(e) NULL)
      if (!is.null(ctr_raw) && !is.null(ctr_adj)) {
        ab_col <- if ("Antibody" %in% colnames(ctr_raw))
          as.character(ctr_raw$Antibody) else rep("ALL", nrow(ctr_raw))
        pw_list[["condition"]] <- data.frame(
          Metric       = metric_name,
          Module       = "3C_ColocLME_Filtered",
          ContrastType = "V1 vs V2",
          Antibody     = ab_col,
          Contrast     = as.character(ctr_raw$contrast),
          Estimate     = ctr_raw$estimate,
          SE           = ctr_raw$SE,
          t_ratio      = ctr_raw$t.ratio,
          df           = ctr_raw$df,
          P_raw        = ctr_raw$p.value,
          P_adjusted   = ctr_adj$p.value,
          Significance = sapply(ctr_adj$p.value, p_to_stars),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (n_ab >= 2) {
    emm_ab <- if (n_cd >= 2) {
      tryCatch(emmeans(fit, ~ Antibody | condition), error = function(e) NULL)
    } else {
      tryCatch(emmeans(fit, ~ Antibody), error = function(e) NULL)
    }
    if (!is.null(emm_ab)) {
      ctr_raw <- tryCatch(as.data.frame(summary(contrast(emm_ab,
                                                         method = "pairwise", adjust = "none"))), error = function(e) NULL)
      ctr_adj <- tryCatch(as.data.frame(summary(contrast(emm_ab,
                                                         method = "pairwise", adjust = ADJUST_METHOD))), error = function(e) NULL)
      if (!is.null(ctr_raw) && !is.null(ctr_adj)) {
        cd_col <- if ("condition" %in% colnames(ctr_raw))
          as.character(ctr_raw$condition) else rep("ALL", nrow(ctr_raw))
        pw_list[["antibody"]] <- data.frame(
          Metric       = metric_name,
          Module       = "3C_ColocLME_Filtered",
          ContrastType = "Syn1 vs EPR20535",
          Antibody     = cd_col,   # here the Antibody column holds the condition stratum
          Contrast     = as.character(ctr_raw$contrast),
          Estimate     = ctr_raw$estimate,
          SE           = ctr_raw$SE,
          t_ratio      = ctr_raw$t.ratio,
          df           = ctr_raw$df,
          P_raw        = ctr_raw$p.value,
          P_adjusted   = ctr_adj$p.value,
          Significance = sapply(ctr_adj$p.value, p_to_stars),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  pw_df <- bind_rows(pw_list)
  if (nrow(pw_df) == 0) pw_df <- NULL
  if (!is.null(pw_df)) pw_df <- add_transparency_cols(pw_df, model_obj)
  
  list(anova = aov_tab, marginals = emm_df, pairwise = pw_df,
       n_cells = nrow(d), n_fovs = n_distinct(d$image_id),
       used_lm = model_obj$used_lm, is_correlation = is_corr,
       filtered_data = filt_df)
}

# ==============================================================================
# MODULE 4 — CONTINUOUS DOSE-RESPONSE (IAV intensity → aSyn)
# ==============================================================================
extract_iav_slope <- function(fit, term_name = "log_iav_mean") {
  coefs <- tryCatch({
    if (inherits(fit, "merMod")) {
      as.data.frame(summary(fit)$coefficients)
    } else {
      as.data.frame(summary(fit)$coefficients)
    }
  }, error = function(e) NULL)
  if (is.null(coefs) || nrow(coefs) == 0) return(NULL)
  coefs$Term <- rownames(coefs)
  rownames(coefs) <- NULL
  coefs
}

# ─── 4a — Per-strain slope ───
run_module4a_per_strain <- function(cells_df, vcol, log_transform = FALSE,
                                    metric_name = "metric", virus_label = "V1") {
  sub <- cells_df %>%
    filter(condition == virus_label,
           IAV_state == "Infected",
           !is.na(.data[[vcol]]), is.finite(.data[[vcol]]),
           !is.na(iav_mean), is.finite(iav_mean), iav_mean > 0)
  if (nrow(sub) < MIN_CELLS_FOR_LMER) return(NULL)
  
  sub$log_iav_mean <- log10(sub$iav_mean)
  
  n_ab <- n_distinct(sub$Antibody)
  fixed_str <- if (n_ab >= 2) "log_iav_mean + Antibody" else "log_iav_mean"
  
  model_obj <- fit_lmer_adaptive(sub, vcol, fixed_str, log_transform)
  if (is.null(model_obj)) return(NULL)
  fit <- model_obj$fit; d <- model_obj$data; is_corr <- model_obj$is_correlation
  n_fovs <- n_distinct(d$image_id)
  
  coefs <- extract_iav_slope(fit)
  if (is.null(coefs)) return(NULL)
  coefs$Metric <- metric_name
  coefs$Virus  <- virus_label
  coefs$Module <- "4A_PerStrain_Slope"
  coefs$is_correlation <- is_corr
  
  if (is_corr) {
    coefs$Estimate <- inv_fisher_z(coefs$Estimate)
  }
  
  # Pull out slope row for iav_mean and significance
  slope_row <- coefs[coefs$Term == "log_iav_mean", , drop = FALSE]
  if (nrow(slope_row) > 0) {
    p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(slope_row))[1]
    if (!is.na(p_col)) slope_row$Significance <- sapply(slope_row[[p_col]],
                                                        p_to_stars)
  }
  
  coefs     <- add_transparency_cols(coefs, model_obj, include_cells = TRUE)
  slope_row <- add_transparency_cols(slope_row, model_obj, include_cells = TRUE)
  
  # Bootstrap supplement
  if ("M4" %in% BOOTSTRAP_MODULES) {
    contrast_fn_m4a <- function(model_fit, model_data) {
      coefs_fit <- extract_iav_slope(model_fit)
      if (!is.null(coefs_fit)) {
        val <- coefs_fit$Estimate[coefs_fit$Term == "log_iav_mean"]
        c("log_iav_mean" = val[1])
      } else {
        c("log_iav_mean" = NA_real_)
      }
    }
    
    boot_res <- cluster_bootstrap_contrast(
      df = sub, vcol = vcol, fixed_effects_str = fixed_str,
      contrast_fn = contrast_fn_m4a, log_transform = log_transform,
      n_boot = N_BOOTSTRAP
    )
    
    if (!is.null(boot_res) && !is.null(slope_row) && nrow(slope_row) > 0) {
      slope_row$Boot_Estimate     <- boot_res[["log_iav_mean"]]$mean
      slope_row$Boot_CI_Low       <- boot_res[["log_iav_mean"]]$CI_low
      slope_row$Boot_CI_High      <- boot_res[["log_iav_mean"]]$CI_high
      slope_row$P_boot            <- boot_res[["log_iav_mean"]]$p_boot
      slope_row$Boot_Mode         <- boot_res[["log_iav_mean"]]$boot_mode
      slope_row$Boot_N_Successful <- boot_res[["log_iav_mean"]]$n_successful
    }
  }
  
  list(coefs = coefs, slope = slope_row, data = d, virus = virus_label,
       is_correlation = is_corr, used_lm = model_obj$used_lm, n_fovs = n_fovs)
}
 
# ─── 4b — Strain dynamics (interaction) ───
run_module4b_strain_dynamics <- function(cells_df, vcol, log_transform = FALSE,
                                         metric_name = "metric") {
  sub <- cells_df %>%
    filter(condition %in% c("V1", "V2"),
           IAV_state == "Infected",
           !is.na(.data[[vcol]]), is.finite(.data[[vcol]]),
           !is.na(iav_mean), is.finite(iav_mean), iav_mean > 0)
  if (nrow(sub) < MIN_CELLS_FOR_LMER) return(NULL)
  sub$condition <- droplevels(sub$condition)
  if (n_distinct(sub$condition) < 2) return(NULL)
  
  sub$log_iav_mean <- log10(sub$iav_mean)
  n_ab <- n_distinct(sub$Antibody)
  fixed_str <- if (n_ab >= 2)
    "log_iav_mean * condition + Antibody"
  else
    "log_iav_mean * condition"
  
  model_obj <- fit_lmer_adaptive(sub, vcol, fixed_str, log_transform)
  if (is.null(model_obj)) return(NULL)
  fit <- model_obj$fit; d <- model_obj$data; is_corr <- model_obj$is_correlation
  n_fovs <- n_distinct(d$image_id)
  
  aov_tab <- clean_aov(fit, metric_name, "4B_StrainDynamics")
  
  coefs <- extract_iav_slope(fit)
  if (!is.null(coefs)) {
    coefs$Metric <- metric_name
    coefs$Module <- "4B_StrainDynamics"
    coefs$is_correlation <- is_corr
    if (is_corr) {
      coefs$Estimate <- inv_fisher_z(coefs$Estimate)
    }
    p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coefs))[1]
    if (!is.na(p_col)) coefs$Significance <- sapply(coefs[[p_col]], p_to_stars)
  }
  
  # V1 and V2 slopes via emmeans::emtrends
  trends_df <- NULL
  emt <- tryCatch(
    emmeans::emtrends(fit, ~ condition, var = "log_iav_mean"),
    error = function(e) NULL)
  if (!is.null(emt)) {
    trends_df <- as.data.frame(emt) %>%
      mutate(Metric = metric_name, Module = "4B_StrainDynamics")
    if (is_corr) {
      trend_col <- grep("trend", colnames(trends_df), value = TRUE)[1]
      if (!is.null(trend_col) && length(trend_col) > 0) {
        trends_df[[trend_col]] <- inv_fisher_z(trends_df[[trend_col]])
      }
    }
  }
  
  # Test of V1 vs V2 slope difference
  slope_diff_df <- NULL
  if (!is.null(emt)) {
    diff_test <- tryCatch(as.data.frame(summary(
      contrast(emt, method = "pairwise", adjust = ADJUST_METHOD))),
      error = function(e) NULL)
    if (!is.null(diff_test)) {
      slope_diff_df <- diff_test %>%
        mutate(Metric = metric_name, Module = "4B_StrainDynamics") %>%
        rename(Contrast = contrast, P_adjusted = p.value)
      if (is_corr) {
        slope_diff_df$estimate <- inv_fisher_z(slope_diff_df$estimate)
      }
      slope_diff_df$Significance <- sapply(slope_diff_df$P_adjusted,
                                           p_to_stars)
    }
  }
  
  aov_tab       <- add_transparency_cols(aov_tab, model_obj, include_cells = TRUE)
  coefs         <- add_transparency_cols(coefs, model_obj, include_cells = TRUE)
  trends_df     <- add_transparency_cols(trends_df, model_obj, include_cells = TRUE)
  slope_diff_df <- add_transparency_cols(slope_diff_df, model_obj, include_cells = TRUE)
  
  # Bootstrap supplement
  if ("M4" %in% BOOTSTRAP_MODULES) {
    contrast_fn_m4b <- function(model_fit, model_data) {
      emt_fit <- tryCatch(
        emmeans::emtrends(model_fit, ~ condition, var = "log_iav_mean"),
        error = function(e) NULL
      )
      if (is.null(emt_fit)) {
        return(c(V1 = NA_real_, V2 = NA_real_, diff = NA_real_))
      }
      trends_df_fit <- as.data.frame(emt_fit)
      diff_df_fit <- tryCatch(
        as.data.frame(summary(contrast(emt_fit, method = "pairwise"))),
        error = function(e) NULL
      )
      
      trend_col <- grep("trend", colnames(trends_df_fit), value = TRUE)[1]
      v1_val <- NA_real_
      v2_val <- NA_real_
      if (!is.null(trend_col) && length(trend_col) > 0) {
        v1_val <- trends_df_fit[[trend_col]][trends_df_fit$condition == "V1"][1]
        v2_val <- trends_df_fit[[trend_col]][trends_df_fit$condition == "V2"][1]
      }
      
      diff_val <- NA_real_
      if (!is.null(diff_df_fit) && nrow(diff_df_fit) > 0) {
        diff_val <- diff_df_fit$estimate[1]
      }
      
      c(V1 = v1_val, V2 = v2_val, diff = diff_val)
    }
    
    boot_res <- cluster_bootstrap_contrast(
      df = sub, vcol = vcol, fixed_effects_str = fixed_str,
      contrast_fn = contrast_fn_m4b, log_transform = log_transform,
      n_boot = N_BOOTSTRAP
    )
    
    if (!is.null(boot_res)) {
      if (!is.null(trends_df) && nrow(trends_df) > 0) {
        trends_df$Boot_Estimate     <- NA_real_
        trends_df$Boot_CI_Low       <- NA_real_
        trends_df$Boot_CI_High      <- NA_real_
        trends_df$P_boot            <- NA_real_
        trends_df$Boot_Mode         <- NA_character_
        trends_df$Boot_N_Successful <- NA_integer_
        
        for (i in 1:nrow(trends_df)) {
          cond <- as.character(trends_df$condition[i])
          if (cond %in% names(boot_res)) {
            trends_df$Boot_Estimate[i]     <- boot_res[[cond]]$mean
            trends_df$Boot_CI_Low[i]       <- boot_res[[cond]]$CI_low
            trends_df$Boot_CI_High[i]      <- boot_res[[cond]]$CI_high
            trends_df$P_boot[i]            <- boot_res[[cond]]$p_boot
            trends_df$Boot_Mode[i]         <- boot_res[[cond]]$boot_mode
            trends_df$Boot_N_Successful[i] <- boot_res[[cond]]$n_successful
          }
        }
      }
      
      if (!is.null(slope_diff_df) && nrow(slope_diff_df) > 0) {
        slope_diff_df$Boot_Estimate     <- boot_res[["diff"]]$mean
        slope_diff_df$Boot_CI_Low       <- boot_res[["diff"]]$CI_low
        slope_diff_df$Boot_CI_High      <- boot_res[["diff"]]$CI_high
        slope_diff_df$P_boot            <- boot_res[["diff"]]$p_boot
        slope_diff_df$Boot_Mode         <- boot_res[["diff"]]$boot_mode
        slope_diff_df$Boot_N_Successful <- boot_res[["diff"]]$n_successful
      }
    }
  }
  
  list(anova = aov_tab, coefs = coefs, trends = trends_df,
       slope_diff = slope_diff_df, data = d, is_correlation = is_corr,
       used_lm = model_obj$used_lm, n_fovs = n_fovs)
}

# ==============================================================================
# EXECUTE PIPELINE
# ==============================================================================

results_M1     <- list()
results_M2_V1  <- list()
results_M2_V2  <- list()
results_M3A    <- list()
results_M3B    <- list()
results_M3C    <- list()
results_M4A_V1 <- list()
results_M4A_V2 <- list()
results_M4B    <- list()

# ─── Module 3C.a — frequency test (run once) ───
m3c_freq <- if (HAS_INFECTED) run_module3c_frequency(cells) else NULL
if (!is.null(m3c_freq)) {
  message("\n=== Module 3C.a — Frequency of true colocalization ===")
  print(m3c_freq$frequency)
  if (!is.null(m3c_freq$tests) && nrow(m3c_freq$tests) > 0) {
    message("\nFrequency test results:")
    print(m3c_freq$tests)
  }
  write_csv(m3c_freq$frequency,
            file.path(out_dir_data, "M3C_coloc_frequency.csv"))
  if (!is.null(m3c_freq$tests))
    write_csv(m3c_freq$tests,
              file.path(out_dir_data, "M3C_coloc_frequency_tests.csv"))
}

# ─── Loop over metrics ───
for (mname in names(metric_configs)) {
  cfg <- metric_configs[[mname]]
  is_absolute <- mname %in% c("asyn_nuc", "asyn_cyto")
  is_coloc    <- cfg$type == "coloc"
  
  message("\n--- Metric: ", mname, " (", cfg$type, ") ---")
  
  # ── Module 1 — abundance metrics; absolute intensities split by antibody ──
  if (!is_coloc) {
    if (is_absolute) {
      for (ab in abs_present) {
        sub <- cells %>% filter(Antibody == ab)
        key <- paste0(mname, "__", ab)
        r <- run_module1(sub, cfg$col, cfg$log,
                         metric_name = paste0(mname, " [", ab, "]"))
        if (!is.null(r)) {
          results_M1[[key]] <- r
          message("  OK M1 ", key, ": n_cells=", r$n_cells)
        }
      }
    } else {
      r <- run_module1(cells, cfg$col, cfg$log, mname)
      if (!is.null(r)) {
        results_M1[[mname]] <- r
        message("  OK M1 ", mname, ": n_cells=", r$n_cells)
      }
    }
  }
  
  # ── Module 2 — abundance metrics, no Control ──
  if (!is_coloc && HAS_INFECTED) {
    for (vir in c("V1", "V2")) {
      if (!vir %in% conds_present) next
      if (is_absolute) {
        for (ab in abs_present) {
          sub <- cells %>% filter(Antibody == ab)
          key <- paste0(mname, "__", ab)
          r <- run_module2(sub, cfg$col, cfg$log,
                           paste0(mname, " [", ab, "]"), vir)
          if (!is.null(r)) {
            if (vir == "V1") results_M2_V1[[key]] <- r
            else             results_M2_V2[[key]] <- r
            message("  OK M2 ", vir, " ", key,
                    ": n_cells=", r$n_cells,
                    " (inf=", r$n_infected, ", bys=", r$n_bystander, ")")
          }
        }
      } else {
        r <- run_module2(cells, cfg$col, cfg$log, mname, vir)
        if (!is.null(r)) {
          if (vir == "V1") results_M2_V1[[mname]] <- r
          else             results_M2_V2[[mname]] <- r
          message("  OK M2 ", vir, " ", mname,
                  ": n_cells=", r$n_cells,
                  " (inf=", r$n_infected, ", bys=", r$n_bystander, ")")
        }
      }
    }
  }
  
  # ── Module 3A (Paracrine) — abundance metrics ──
  if (!is_coloc && HAS_INFECTED) {
    if (is_absolute) {
      for (ab in abs_present) {
        sub <- cells %>% filter(Antibody == ab)
        key <- paste0(mname, "__", ab)
        r <- run_module3a_paracrine(sub, cfg$col, cfg$log,
                                    paste0(mname, " [", ab, "]"))
        if (!is.null(r)) {
          results_M3A[[key]] <- r
          message("  OK M3A ", key, ": n_cells=", r$n_cells)
        }
      }
    } else {
      r <- run_module3a_paracrine(cells, cfg$col, cfg$log, mname)
      if (!is.null(r)) {
        results_M3A[[mname]] <- r
        message("  OK M3A ", mname, ": n_cells=", r$n_cells)
      }
    }
  }
  
  # ── Module 3B (Strain interaction) — abundance metrics ──
  if (!is_coloc && HAS_INFECTED) {
    if (is_absolute) {
      for (ab in abs_present) {
        sub <- cells %>% filter(Antibody == ab)
        key <- paste0(mname, "__", ab)
        r <- run_module3b_strain(sub, cfg$col, cfg$log,
                                 paste0(mname, " [", ab, "]"))
        if (!is.null(r)) {
          results_M3B[[key]] <- r
          message("  OK M3B ", key, ": n_cells=", r$n_cells)
        }
      }
    } else {
      r <- run_module3b_strain(cells, cfg$col, cfg$log, mname)
      if (!is.null(r)) {
        results_M3B[[mname]] <- r
        message("  OK M3B ", mname, ": n_cells=", r$n_cells)
      }
    }
  }
  
  # ── Module 3C.b (Coloc LME on Costes-positive cells) — colocalisation metrics ──
  if (is_coloc && HAS_INFECTED) {
    r <- run_module3c_coloc_lme(cells, cfg$col, mname)
    if (!is.null(r)) {
      results_M3C[[mname]] <- r
      message("  OK M3C ", mname, ": n_cells=", r$n_cells,
              " (filtered by Costes p<", COSTES_P_THRESHOLD, ")")
    }
  }
  
  # ── Module 4 (Viral-load dependence) — abundance metrics ──
  if (!is_coloc && HAS_INFECTED && "iav_mean" %in% colnames(cells)) {
    if (is_absolute) {
      for (ab in abs_present) {
        sub <- cells %>% filter(Antibody == ab)
        key <- paste0(mname, "__", ab)
        r1 <- run_module4a_per_strain(sub, cfg$col, cfg$log,
                                      paste0(mname, " [", ab, "]"), "V1")
        if (!is.null(r1)) results_M4A_V1[[key]] <- r1
        r2 <- run_module4a_per_strain(sub, cfg$col, cfg$log,
                                      paste0(mname, " [", ab, "]"), "V2")
        if (!is.null(r2)) results_M4A_V2[[key]] <- r2
        rB <- run_module4b_strain_dynamics(sub, cfg$col, cfg$log,
                                           paste0(mname, " [", ab, "]"))
        if (!is.null(rB)) results_M4B[[key]] <- rB
      }
    } else {
      r1 <- run_module4a_per_strain(cells, cfg$col, cfg$log, mname, "V1")
      if (!is.null(r1)) results_M4A_V1[[mname]] <- r1
      r2 <- run_module4a_per_strain(cells, cfg$col, cfg$log, mname, "V2")
      if (!is.null(r2)) results_M4A_V2[[mname]] <- r2
      rB <- run_module4b_strain_dynamics(cells, cfg$col, cfg$log, mname)
      if (!is.null(rB)) results_M4B[[mname]] <- rB
      if (!is.null(r1) || !is.null(r2) || !is.null(rB))
        message("  OK M4 ", mname)
    }
  }
}

# ==============================================================================
# WRITE COMPREHENSIVE EXCEL LOGS
# ==============================================================================

stats_path <- file.path(out_dir_data,
                        "Confocal_Statistics_LMER.xlsx")

wb <- createWorkbook()

methods_text <- data.frame(
  Section = c("Design & Context", "Replicates", "Strict Stats (Fisher Z)",
              "Module 1", "Module 2 (Main Story)",
              "Module 3A (Paracrine)", "Module 3B (Strain Interaction)",
              "Module 3C (Coloc Profiling)", "Module 4 (Dose-Response)",
              "Random effects", "Multiple comparisons", "Software"),
  Description = c(
    sprintf("Conditions: %s. Antibody: %s. PFF & high_mag excluded upstream.",
            paste(conds_present, collapse = " / "),
            paste(abs_present,   collapse = " / ")),
    sprintf("n=%d biological replicate(s). %s",
            n_reps_global,
            if (n_reps_global == 1)
              "Single-replicate — Module 1 EXPLORATORY only. Modules 2, 3A, 3B leverage within-FOV / within-condition cell-level contrasts."
            else "Replicate variance estimable."),
    "Fisher's z transformation z = 0.5*ln((1+r)/(1-r)) is applied to all metrics ending in '_r' or containing 'pearson_r' before fitting. Marginal means and CIs are back-transformed to the r scale for reporting.",
    "Module 1: lmer(metric ~ condition * Antibody + (1|image_id)). Treatment-level comparison; EXPLORATORY with n = 1 culture (treatment confounded with well).",
    "Module 2: lmer(metric ~ IAV_state * Antibody + (1|image_id)), fitted separately for V1 and V2. Within-FOV contrast, robust to well/batch effects. Colocalisation metrics excluded (no IAV signal in bystanders). No Control.",
    "Module 3A: lmer(metric ~ para_group * Antibody + (1|image_id)) with para_group = {Control, V1.Bystander, V2.Bystander}. Tests a paracrine effect on uninfected bystander cells.",
    "Module 3B: lmer(metric ~ condition * IAV_state + Antibody + (1|image_id)) on V1 and V2 only (Control has no infected cells). Tests whether the infection effect differs between strains.",
    sprintf("Module 3C: (a) Frequency: chi-square / Fisher exact tests of the %% of infected cells with true colocalisation (Costes p<%.2f), V1 vs V2 and Syn1 vs EPR20535. (b) LME restricted to Costes-positive cells, comparing colocalisation strength (Pearson r, Fisher z) between strains and antibodies.", COSTES_P_THRESHOLD),
    "Module 4: (a) Per-strain slope: lmer(metric ~ log10(iav_mean) + Antibody + (1|image_id)), V1 and V2 separately. (b) Strain dynamics: lmer(metric ~ log10(iav_mean) * condition + Antibody + (1|image_id)) on V1+V2; slope difference V1 vs V2 via emtrends + pairwise contrast.",
    "Field of view (image_id) as random intercept; lm fallback when only one FOV is available.",
    sprintf("Pairwise contrasts adjusted by '%s' within each metric.", ADJUST_METHOD),
    sprintf("R %s.%s; lme4, lmerTest, emmeans, car, ggsignif, openxlsx",
            R.version$major, R.version$minor)
  )
)
addWorksheet(wb, "0_Methods"); writeData(wb, "0_Methods", methods_text)

# Group summary
group_summary <- cells %>%
  group_by(Antibody, condition) %>%
  summarize(n_cells = n(),
            n_FOVs  = n_distinct(image_id),
            n_replicates = n_distinct(replicate),
            .groups = "drop")
addWorksheet(wb, "1_Group_Summary")
writeData(wb, "1_Group_Summary", group_summary)

if (HAS_INFECTED) {
  infection_summary <- cells %>%
    filter(condition %in% c("V1", "V2")) %>%
    group_by(Antibody, condition, IAV_state) %>%
    summarize(n_cells = n(), n_FOVs = n_distinct(image_id), .groups = "drop")
  addWorksheet(wb, "1b_Infection_Summary")
  writeData(wb, "1b_Infection_Summary", infection_summary)
}

collect_field <- function(res_list, fld) {
  bind_rows(lapply(res_list, function(r) {
    if (!is.null(r) && !is.null(r[[fld]])) r[[fld]] else NULL
  }))
}

# Module 1 sheets
m1_anova <- collect_field(results_M1, "anova")
m1_pw    <- collect_field(results_M1, "pairwise")
m1_emm   <- collect_field(results_M1, "marginals")
if (nrow(m1_anova) > 0) { addWorksheet(wb, "M1_ANOVA");     writeData(wb, "M1_ANOVA",     m1_anova) }
if (nrow(m1_pw)    > 0) { addWorksheet(wb, "M1_Pairwise");  writeData(wb, "M1_Pairwise",  m1_pw) }
if (nrow(m1_emm)   > 0) { addWorksheet(wb, "M1_Marginals"); writeData(wb, "M1_Marginals", m1_emm) }

# Module 2 sheets (V1 and V2 separately)
m2v1_anova <- collect_field(results_M2_V1, "anova")
m2v1_pw    <- collect_field(results_M2_V1, "pairwise")
m2v1_emm   <- collect_field(results_M2_V1, "marginals")
m2v2_anova <- collect_field(results_M2_V2, "anova")
m2v2_pw    <- collect_field(results_M2_V2, "pairwise")
m2v2_emm   <- collect_field(results_M2_V2, "marginals")
if (nrow(m2v1_anova) > 0) { addWorksheet(wb, "M2_V1_ANOVA");     writeData(wb, "M2_V1_ANOVA",     m2v1_anova) }
if (nrow(m2v1_pw)    > 0) { addWorksheet(wb, "M2_V1_Pairwise");  writeData(wb, "M2_V1_Pairwise",  m2v1_pw) }
if (nrow(m2v1_emm)   > 0) { addWorksheet(wb, "M2_V1_Marginals"); writeData(wb, "M2_V1_Marginals", m2v1_emm) }
if (nrow(m2v2_anova) > 0) { addWorksheet(wb, "M2_V2_ANOVA");     writeData(wb, "M2_V2_ANOVA",     m2v2_anova) }
if (nrow(m2v2_pw)    > 0) { addWorksheet(wb, "M2_V2_Pairwise");  writeData(wb, "M2_V2_Pairwise",  m2v2_pw) }
if (nrow(m2v2_emm)   > 0) { addWorksheet(wb, "M2_V2_Marginals"); writeData(wb, "M2_V2_Marginals", m2v2_emm) }

# Module 3A
m3a_anova <- collect_field(results_M3A, "anova")
m3a_pw    <- collect_field(results_M3A, "pairwise")
m3a_emm   <- collect_field(results_M3A, "marginals")
if (nrow(m3a_anova) > 0) { addWorksheet(wb, "M3A_ANOVA");     writeData(wb, "M3A_ANOVA",     m3a_anova) }
if (nrow(m3a_pw)    > 0) { addWorksheet(wb, "M3A_Pairwise");  writeData(wb, "M3A_Pairwise",  m3a_pw) }
if (nrow(m3a_emm)   > 0) { addWorksheet(wb, "M3A_Marginals"); writeData(wb, "M3A_Marginals", m3a_emm) }

# Module 3B
m3b_anova <- collect_field(results_M3B, "anova")
m3b_pw    <- collect_field(results_M3B, "pairwise")
m3b_emm   <- collect_field(results_M3B, "marginals")
if (nrow(m3b_anova) > 0) { addWorksheet(wb, "M3B_ANOVA");     writeData(wb, "M3B_ANOVA",     m3b_anova) }
if (nrow(m3b_pw)    > 0) { addWorksheet(wb, "M3B_Pairwise");  writeData(wb, "M3B_Pairwise",  m3b_pw) }
if (nrow(m3b_emm)   > 0) { addWorksheet(wb, "M3B_Marginals"); writeData(wb, "M3B_Marginals", m3b_emm) }

# Module 3C
if (!is.null(m3c_freq)) {
  addWorksheet(wb, "M3C_Freq_Summary")
  writeData(wb, "M3C_Freq_Summary", m3c_freq$frequency)
  if (!is.null(m3c_freq$tests) && nrow(m3c_freq$tests) > 0) {
    addWorksheet(wb, "M3C_Freq_Tests")
    writeData(wb, "M3C_Freq_Tests", m3c_freq$tests)
  }
}
m3c_anova <- collect_field(results_M3C, "anova")
m3c_pw    <- collect_field(results_M3C, "pairwise")
m3c_emm   <- collect_field(results_M3C, "marginals")
if (nrow(m3c_anova) > 0) { addWorksheet(wb, "M3C_LME_ANOVA");     writeData(wb, "M3C_LME_ANOVA",     m3c_anova) }
if (nrow(m3c_pw)    > 0) { addWorksheet(wb, "M3C_LME_Pairwise");  writeData(wb, "M3C_LME_Pairwise",  m3c_pw) }
if (nrow(m3c_emm)   > 0) { addWorksheet(wb, "M3C_LME_Marginals"); writeData(wb, "M3C_LME_Marginals", m3c_emm) }

# Module 4
m4a_v1_slopes <- collect_field(results_M4A_V1, "slope")
m4a_v2_slopes <- collect_field(results_M4A_V2, "slope")
m4a_v1_coefs  <- collect_field(results_M4A_V1, "coefs")
m4a_v2_coefs  <- collect_field(results_M4A_V2, "coefs")
m4b_anova     <- collect_field(results_M4B, "anova")
m4b_coefs     <- collect_field(results_M4B, "coefs")
m4b_trends    <- collect_field(results_M4B, "trends")
m4b_slopediff <- collect_field(results_M4B, "slope_diff")
if (nrow(m4a_v1_slopes) > 0) { addWorksheet(wb, "M4A_V1_Slopes"); writeData(wb, "M4A_V1_Slopes", m4a_v1_slopes) }
if (nrow(m4a_v2_slopes) > 0) { addWorksheet(wb, "M4A_V2_Slopes"); writeData(wb, "M4A_V2_Slopes", m4a_v2_slopes) }
if (nrow(m4a_v1_coefs)  > 0) { addWorksheet(wb, "M4A_V1_AllCoefs"); writeData(wb, "M4A_V1_AllCoefs", m4a_v1_coefs) }
if (nrow(m4a_v2_coefs)  > 0) { addWorksheet(wb, "M4A_V2_AllCoefs"); writeData(wb, "M4A_V2_AllCoefs", m4a_v2_coefs) }
if (nrow(m4b_anova)     > 0) { addWorksheet(wb, "M4B_ANOVA");     writeData(wb, "M4B_ANOVA",     m4b_anova) }
if (nrow(m4b_coefs)     > 0) { addWorksheet(wb, "M4B_Coefs");     writeData(wb, "M4B_Coefs",     m4b_coefs) }
if (nrow(m4b_trends)    > 0) { addWorksheet(wb, "M4B_Trends");    writeData(wb, "M4B_Trends",    m4b_trends) }
if (nrow(m4b_slopediff) > 0) { addWorksheet(wb, "M4B_SlopeDiff"); writeData(wb, "M4B_SlopeDiff", m4b_slopediff) }

saveWorkbook(wb, stats_path, overwrite = TRUE)
message("\nStats log saved: ", stats_path)

# ==============================================================================
# PLOTTING
# ==============================================================================

# Generic bar plot (mean + SEM + jitter + significance stars)
render_bar_plot <- function(df, vcol, gcol, metric_label, subtitle_str,
                            pw_df, color_map, anova_tab = NULL) {
  d_cells <- df %>% filter(!is.na(.data[[vcol]]), is.finite(.data[[vcol]]),
                           !is.na(.data[[gcol]]))
  if (nrow(d_cells) == 0) return(NULL)
  
  ylim_clip <- compute_clip_ylim(d_cells, vcol, gcol)
  has_multi_ab <- "Antibody" %in% colnames(d_cells) &&
    n_distinct(d_cells$Antibody) >= 2
  
  sdf <- if (has_multi_ab) {
    d_cells %>%
      group_by(.data[[gcol]], Antibody) %>%
      summarize(m = mean(.data[[vcol]], na.rm = TRUE),
                s = sd(.data[[vcol]], na.rm = TRUE) /
                  sqrt(sum(!is.na(.data[[vcol]]))),
                .groups = "drop")
  } else {
    d_cells %>%
      group_by(.data[[gcol]]) %>%
      summarize(m = mean(.data[[vcol]], na.rm = TRUE),
                s = sd(.data[[vcol]], na.rm = TRUE) /
                  sqrt(sum(!is.na(.data[[vcol]]))),
                .groups = "drop")
  }
  
  anova_caption <- if (!is.null(anova_tab) && nrow(anova_tab) > 0) {
    p_col <- if ("Pr(>F)" %in% colnames(anova_tab)) "Pr(>F)" else "Pr(>Chisq)"
    f_col <- if ("F value" %in% colnames(anova_tab)) "F value" else "Chisq"
    anova_tab %>% filter(!is.na(.data[[p_col]])) %>%
      mutate(text = sprintf("%s: F=%.2f, p=%.3g %s",
                            gsub(":", " x ", Term),
                            .data[[f_col]], .data[[p_col]],
                            Significance)) %>%
      pull(text) %>% paste(collapse = "  |  ")
  } else ""
  
  p <- ggplot(d_cells, aes(x = .data[[gcol]], y = .data[[vcol]])) +
    geom_col(data = sdf,
             aes(x = .data[[gcol]], y = m, fill = .data[[gcol]]),
             alpha = 0.55, width = 0.7, color = NA) +
    geom_errorbar(data = sdf,
                  aes(x = .data[[gcol]], y = m,
                      ymin = m - s, ymax = m + s,
                      color = .data[[gcol]]),
                  width = 0.22, linewidth = 0.6, inherit.aes = FALSE) +
    geom_jitter(aes(color = .data[[gcol]]),
                shape = 21, fill = "white", size = 1.5, stroke = 0.4,
                alpha = 0.5, width = 0.18) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_color_manual(values = color_map, drop = FALSE) +
    labs(y = metric_label, x = NULL, title = metric_label,
         subtitle = subtitle_str, caption = anova_caption) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          plot.caption  = element_text(size = 8, color = "grey25",
                                       hjust = 0, face = "italic"),
          legend.position = "none")
  if (has_multi_ab) p <- p + facet_wrap(~ Antibody, scales = "free_y")
  
  if (!is.null(ylim_clip) && all(is.finite(ylim_clip))) {
    span <- diff(ylim_clip)
    p <- p + coord_cartesian(
      ylim = c(ylim_clip[1], ylim_clip[2] + 0.30 * span), clip = "off")
  } else {
    p <- p + coord_cartesian(clip = "off")
  }
  
  if (!is.null(pw_df) && nrow(pw_df) > 0) {
    if (has_multi_ab) {
      for (ab in abs_present) {
        sig <- get_sig_pairs(pw_df, antibody_filter = ab)
        if (length(sig$pairs) > 0) {
          p <- p + ggsignif::geom_signif(
            comparisons = sig$pairs, annotations = sig$labels,
            step_increase = 0.09, tip_length = 0.01, vjust = 0.4,
            textsize = 3.5,
            data = subset(d_cells, Antibody == ab))
        }
      }
    } else {
      sig <- get_sig_pairs(pw_df)
      if (length(sig$pairs) > 0) {
        p <- p + ggsignif::geom_signif(
          comparisons = sig$pairs, annotations = sig$labels,
          step_increase = 0.09, tip_length = 0.01, vjust = 0.4,
          textsize = 3.5)
      }
    }
  }
  p
}

message("\n", paste(rep("=", 70), collapse = ""))
message("RENDERING FIGURES")
message(paste(rep("=", 70), collapse = ""))

# ─── Render bar plots ───
for (mname in names(metric_configs)) {
  cfg <- metric_configs[[mname]]
  is_absolute <- mname %in% c("asyn_nuc", "asyn_cyto")
  is_coloc    <- cfg$type == "coloc"
  
  # Module 1
  if (!is_coloc) {
    if (is_absolute) {
      for (ab in abs_present) {
        key <- paste0(mname, "__", ab)
        r <- results_M1[[key]]; if (is.null(r)) next
        p <- render_bar_plot(
          cells %>% filter(Antibody == ab),
          cfg$col, "condition", paste0(cfg$label, " [", ab, "]"),
          paste("Module 1 [EXPLORATORY] |", EXPLORATORY_CAVEAT),
          r$pairwise, group_colors, r$anova)
        if (!is.null(p)) save_plot_pdf_png(p, paste0("M1__", mname, "__", ab),
                                           width = 8.5, height = 6.5)
      }
    } else {
      r <- results_M1[[mname]]
      if (!is.null(r)) {
        p <- render_bar_plot(cells, cfg$col, "condition", cfg$label,
                             paste("Module 1 [EXPLORATORY] |", EXPLORATORY_CAVEAT),
                             r$pairwise, group_colors, r$anova)
        if (!is.null(p)) save_plot_pdf_png(p, paste0("M1__", mname),
                                           width = 8.5, height = 6.5)
      }
    }
  }
  
  # Module 2 — V1 and V2 separately
  if (!is_coloc && HAS_INFECTED) {
    for (vir in c("V1", "V2")) {
      res_list <- if (vir == "V1") results_M2_V1 else results_M2_V2
      if (is_absolute) {
        for (ab in abs_present) {
          key <- paste0(mname, "__", ab)
          r <- res_list[[key]]; if (is.null(r)) next
          p <- render_bar_plot(
            cells %>% filter(Antibody == ab, condition == vir),
            cfg$col, "IAV_state",
            paste0(cfg$label, " [", ab, "]"),
            sprintf("Module 2 [PRIMARY] — %s | Within-FOV (Infected vs Bystander)",
                    vir),
            r$pairwise, infection_colors, r$anova)
          if (!is.null(p)) save_plot_pdf_png(
            p, paste0("M2_", vir, "__", mname, "__", ab),
            width = 8.5, height = 6.5)
        }
      } else {
        r <- res_list[[mname]]; if (is.null(r)) next
        p <- render_bar_plot(
          cells %>% filter(condition == vir),
          cfg$col, "IAV_state", cfg$label,
          sprintf("Module 2 [PRIMARY] — %s | Within-FOV (Infected vs Bystander)",
                  vir),
          r$pairwise, infection_colors, r$anova)
        if (!is.null(p)) save_plot_pdf_png(p, paste0("M2_", vir, "__", mname),
                                           width = 8.5, height = 6.5)
      }
    }
  }
  
  # Module 3A — Paracrine
  if (!is_coloc && HAS_INFECTED) {
    para_df <- cells %>%
      mutate(para_group = case_when(
        condition == "Control" ~ "Control",
        condition == "V1" & IAV_state == "Bystander" ~ "V1.Bystander",
        condition == "V2" & IAV_state == "Bystander" ~ "V2.Bystander",
        TRUE ~ NA_character_)) %>%
      filter(!is.na(para_group)) %>%
      mutate(para_group = factor(para_group,
                                 levels = c("Control", "V1.Bystander",
                                            "V2.Bystander")))
    if (is_absolute) {
      for (ab in abs_present) {
        key <- paste0(mname, "__", ab)
        r <- results_M3A[[key]]; if (is.null(r)) next
        p <- render_bar_plot(
          para_df %>% filter(Antibody == ab),
          cfg$col, "para_group", paste0(cfg$label, " [", ab, "]"),
          "Module 3A [Paracrine Effect] | Control vs V1.Bystander vs V2.Bystander",
          r$pairwise, paracrine_colors, r$anova)
        if (!is.null(p)) save_plot_pdf_png(
          p, paste0("M3A__", mname, "__", ab),
          width = 8.5, height = 6.5)
      }
    } else {
      r <- results_M3A[[mname]]
      if (!is.null(r)) {
        p <- render_bar_plot(
          para_df, cfg$col, "para_group", cfg$label,
          "Module 3A [Paracrine Effect] | Control vs V1.Bystander vs V2.Bystander",
          r$pairwise, paracrine_colors, r$anova)
        if (!is.null(p)) save_plot_pdf_png(p, paste0("M3A__", mname),
                                           width = 8.5, height = 6.5)
      }
    }
  }
  
  # Module 3B — Strain interaction (cond_state combined)
  if (!is_coloc && HAS_INFECTED) {
    strain_df <- cells %>%
      filter(condition %in% c("V1", "V2"), !is.na(IAV_state)) %>%
      mutate(cond_state = factor(paste(condition, IAV_state, sep = "."),
                                 levels = c("V1.Bystander", "V1.Infected",
                                            "V2.Bystander", "V2.Infected")))
    strain_colors <- c("V1.Bystander" = "#A8C7F0", "V1.Infected" = "#1F61B0",
                       "V2.Bystander" = "#F2B4B4", "V2.Infected" = "#A02A2A")
    if (is_absolute) {
      for (ab in abs_present) {
        key <- paste0(mname, "__", ab)
        r <- results_M3B[[key]]; if (is.null(r)) next
        p <- render_bar_plot(
          strain_df %>% filter(Antibody == ab),
          cfg$col, "cond_state", paste0(cfg$label, " [", ab, "]"),
          "Module 3B [Strain Interaction] | V1 vs V2 across IAV state",
          r$pairwise, strain_colors, r$anova)
        if (!is.null(p)) save_plot_pdf_png(
          p, paste0("M3B__", mname, "__", ab),
          width = 9.5, height = 6.5)
      }
    } else {
      r <- results_M3B[[mname]]
      if (!is.null(r)) {
        p <- render_bar_plot(
          strain_df, cfg$col, "cond_state", cfg$label,
          "Module 3B [Strain Interaction] | V1 vs V2 across IAV state",
          r$pairwise, strain_colors, r$anova)
        if (!is.null(p)) save_plot_pdf_png(p, paste0("M3B__", mname),
                                           width = 9.5, height = 6.5)
      }
    }
  }
  
  # Module 3C — colocalisation profiling (colocalisation metrics only)
  if (is_coloc && HAS_INFECTED) {
    r <- results_M3C[[mname]]
    if (!is.null(r) && !is.null(r$filtered_data)) {
      p <- render_bar_plot(
        r$filtered_data, cfg$col, "condition", cfg$label,
        sprintf("Module 3C [Coloc, Costes-filtered p<%.2f] | V1 vs V2",
                COSTES_P_THRESHOLD),
        r$pairwise, group_colors, r$anova)
      if (!is.null(p)) save_plot_pdf_png(p, paste0("M3C__", mname),
                                         width = 8.5, height = 6.5)
    }
  }
}

# ─── Module 4 renderer: scatter plots against log10(iav_mean) ───
if (HAS_INFECTED && "iav_mean" %in% colnames(cells)) {
  for (mname in names(metric_configs)) {
    cfg <- metric_configs[[mname]]
    if (cfg$type == "coloc") next   # viral-load models for abundance metrics only
    is_absolute <- mname %in% c("asyn_nuc", "asyn_cyto")
    
    if (is_absolute) {
      ab_iter <- abs_present
    } else {
      ab_iter <- c(NA_character_)
    }
    
    for (ab in ab_iter) {
      key <- if (is.na(ab)) mname else paste0(mname, "__", ab)
      r_v1 <- results_M4A_V1[[key]]
      r_v2 <- results_M4A_V2[[key]]
      r_b  <- results_M4B[[key]]
      
      # At least one successful model is needed to draw the scatter
      if (is.null(r_v1) && is.null(r_v2) && is.null(r_b)) next
      
      # Combined data for the scatter (V1 + V2)
      plot_df <- cells %>%
        filter(condition %in% c("V1", "V2"),
               IAV_state == "Infected",
               !is.na(.data[[cfg$col]]), is.finite(.data[[cfg$col]]),
               !is.na(iav_mean), iav_mean > 0)
      if (!is.na(ab)) plot_df <- plot_df %>% filter(Antibody == ab)
      if (nrow(plot_df) < 10) next
      plot_df$condition <- droplevels(plot_df$condition)
      
      # Subtitle with slope information
      subtitle_parts <- c()
      for (vir in c("V1", "V2")) {
        rr <- if (vir == "V1") r_v1 else r_v2
        if (!is.null(rr) && !is.null(rr$slope) && nrow(rr$slope) > 0) {
          slope_val <- rr$slope$Estimate[1]
          p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)"),
                             colnames(rr$slope))[1]
          p_val <- if (!is.na(p_col)) rr$slope[[p_col]][1] else NA
          subtitle_parts <- c(subtitle_parts,
                              sprintf("%s slope=%.3f (p=%.3g %s)", vir, slope_val, p_val,
                                      p_to_stars(p_val)))
        }
      }
      # Slope diff test (V1 vs V2)
      if (!is.null(r_b) && !is.null(r_b$slope_diff) &&
          nrow(r_b$slope_diff) > 0) {
        sd_row <- r_b$slope_diff[1, ]
        subtitle_parts <- c(subtitle_parts,
                            sprintf("Δslope p=%.3g %s",
                                    sd_row$P_adjusted, sd_row$Significance))
      }
      subtitle_str <- if (length(subtitle_parts) > 0)
        paste(subtitle_parts, collapse = "  |  ")
      else ""
      
      ylim_clip <- compute_clip_ylim(plot_df, cfg$col, "condition")
      
      p4 <- ggplot(plot_df, aes(x = iav_mean, y = .data[[cfg$col]],
                                color = condition)) +
        geom_point(aes(shape = IAV_state), alpha = 0.45, size = 1.5) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
        scale_color_manual(values = group_colors[c("V1", "V2")]) +
        scale_x_log10() +
        labs(x = "Continuous IAV intensity (mean, log10)",
             y = if (is.na(ab)) cfg$label else paste0(cfg$label, " [", ab, "]"),
             title = paste(cfg$label,
                           if (is.na(ab)) "" else paste0("[", ab, "]"),
                           " — Dose-Response (Module 4)"),
             subtitle = paste("Continuous LME |", subtitle_str)) +
        theme_classic(base_size = 11) +
        theme(legend.position = "bottom",
              plot.subtitle = element_text(size = 9, color = "grey25",
                                           face = "italic"))
      
      if (!is.na(ab)) {
        # Already filtered to one antibody, no facet needed
      } else if (n_distinct(plot_df$Antibody) >= 2) {
        p4 <- p4 + facet_wrap(~ Antibody, scales = "free_y")
      }
      if (!is.null(ylim_clip) && all(is.finite(ylim_clip))) {
        p4 <- p4 + coord_cartesian(ylim = ylim_clip)
      }
      
      fname <- if (is.na(ab))
        paste0("M4_DoseResponse__", mname)
      else
        paste0("M4_DoseResponse__", mname, "__", ab)
      save_plot_pdf_png(p4, fname, width = 9.5, height = 6.5)
    }
  }
}

# ─── Module 3C — Frequency bar chart ───
if (!is.null(m3c_freq) && nrow(m3c_freq$frequency) > 0) {
  freq_df <- m3c_freq$frequency
  p_freq <- ggplot(freq_df,
                   aes(x = condition, y = Pct_True_Coloc, fill = condition)) +
    geom_col(width = 0.7, alpha = 0.7) +
    geom_text(aes(label = sprintf("%d/%d\n(%.1f%%)",
                                  True_Coloc_Cells, Total_Infected_Cells,
                                  Pct_True_Coloc)),
              vjust = -0.3, size = 3) +
    scale_fill_manual(values = group_colors) +
    facet_wrap(~ Antibody) +
    labs(x = NULL, y = "% Infected cells with true coloc",
         title = sprintf("Module 3C — Frequency of true colocalization (Costes p<%.2f)",
                         COSTES_P_THRESHOLD),
         subtitle = "Chi-square / Fisher exact test results in the M3C_Freq_Tests sheet") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 9, color = "grey40")) +
    coord_cartesian(ylim = c(0, max(freq_df$Pct_True_Coloc) * 1.25))
  save_plot_pdf_png(p_freq, "M3C_frequency_summary",
                    width = 8.5, height = 6)
}

# ==============================================================================
# DONE
# ==============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("DONE — confocal single-cell statistics complete")
message(paste(rep("=", 70), collapse = ""))
message("Stats log : ", stats_path)
message("PDF figs  : ", out_dir_pdf)
message("PNG figs  : ", out_dir_png)
message("\nReading guide:")
message("  M1 — Treatment (Control/V1/V2) bulk — EXPLORATORY")
message("  M2 — Within-FOV Infected vs Bystander, V1 and V2 separately — PRIMARY")
message("  M3A — Paracrine: Control vs V1.Bys vs V2.Bys")
message("  M3B — Strain interaction: V1 vs V2 (condition * IAV_state, no Control)")
message("  M3C — Colocalization profiling:")
message("        (a) Frequency test (chi-square/Fisher) — coloc rate")
message("        (b) LME on Costes-positive cells (Fisher z for r)")
message("  M4 — Continuous dose-response: slope of metric ~ log10(iav_mean)")
message("       (a) Per-strain slope (V1, V2 separately)")
message("       (b) Strain dynamics: emtrends comparison of V1 vs V2 slopes")
if (n_reps_global == 1) {
  message("\n*** REMINDER: n=1 replicate. M1 EXPLORATORY only; M2/M3A/M3B/M3C/M4 ",
          "leverage within-condition cell-level contrasts.")
}
