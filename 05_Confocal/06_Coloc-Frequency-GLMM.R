# ==============================================================================
# 06_Coloc-Frequency-GLMM.R
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(emmeans)
  library(car)
  library(openxlsx)
})

# ─── COMMAND-LINE ARGS ───
args <- commandArgs(trailingOnly = TRUE)
SEG_DIR           <- file.path("Confocal", "results_segmentation")
master_cells_csv  <- file.path(SEG_DIR, "master_cells.csv")
master_coloc_csv  <- file.path(SEG_DIR, "master_coloc.csv")
exclude_csv       <- NULL
out_xlsx          <- file.path("Confocal", "results_confocal_stats", "data",
                               "M3C_Freq_GLMM_supplement.xlsx")

for (i in seq_along(args)) {
  if (args[i] == "--master-cells" && i < length(args)) master_cells_csv <- args[i + 1]
  if (args[i] == "--master-coloc" && i < length(args)) master_coloc_csv <- args[i + 1]
  if (args[i] == "--exclude"      && i < length(args)) exclude_csv      <- args[i + 1]
  if (args[i] == "--out-xlsx"     && i < length(args)) out_xlsx         <- args[i + 1]
}

# ─── LOAD & PREPROCESS DATA ───
message("Loading master cells...")
if (!file.exists(master_cells_csv)) {
  stop("master_cells.csv not found: ", master_cells_csv)
}
cells_raw <- read_csv(master_cells_csv, show_col_types = FALSE)
message("Loaded ", nrow(cells_raw), " cells from ", n_distinct(cells_raw$image_id), " images.")

cells <- cells_raw %>%
  mutate(
    antibody = case_when(
      grepl("Syn1",              antibody, ignore.case = TRUE) ~ "Syn1",
      grepl("EPR|Abcam|C-?term", antibody, ignore.case = TRUE) ~ "EPR20535",
      TRUE ~ as.character(antibody)
    ),
    condition = factor(condition, levels = c("Control", "V1", "V2")),
    Antibody  = factor(antibody, levels = c("Syn1", "EPR20535"))
  ) %>%
  filter(!is.na(condition), !is.na(Antibody), !is.na(replicate))

if (any(grepl("PFF", as.character(cells$condition)))) {
  cells <- cells %>% filter(!grepl("PFF", as.character(condition)))
}

if (!is.null(exclude_csv) && file.exists(exclude_csv)) {
  excl <- read_csv(exclude_csv, show_col_types = FALSE)
  n_before <- nrow(cells)
  cells <- cells %>% filter(!image_id %in% excl$image_id)
  message("Excluded ", n_before - nrow(cells), " cells from ", length(unique(excl$image_id)), " flagged FOVs")
}

# Merge the Costes p-values unless they are already in the cell table
if (!all(c("nuc_asyn_iav_costes_p", "nuc_asyn_iav_pearson_r") %in% colnames(cells))) {
  if (!file.exists(master_coloc_csv)) {
    stop("Colocalisation columns are missing and master_coloc.csv was not found: ",
         master_coloc_csv, " (run 05_Confocal/03_Colocalization-Analysis.py or pass --master-coloc).")
  }
  coloc <- read_csv(master_coloc_csv, show_col_types = FALSE)
  message("Loaded ", nrow(coloc), " coloc rows; merging into cells")
  coloc_cols <- setdiff(colnames(coloc), c(colnames(cells), "nucleus_id", "image_id"))
  cells <- cells %>%
    left_join(coloc %>% select(nucleus_id, image_id, all_of(coloc_cols)),
              by = c("nucleus_id", "image_id"))
}

# ─── INFECTED CELLS IN V1 AND V2 ───
if (!"infected" %in% colnames(cells)) {
  stop("Column 'infected' not found in cells data.")
}

d_analysis <- cells %>%
  filter(condition %in% c("V1", "V2")) %>%
  filter(as.logical(infected) == TRUE) %>%
  mutate(
    condition = factor(condition, levels = c("V1", "V2")),
    antibody  = factor(antibody, levels = c("Syn1", "EPR20535")),
    # Binary outcome: Costes p < 0.05 and a defined Pearson r
    true_coloc = (nuc_asyn_iav_costes_p < 0.05) & !is.na(nuc_asyn_iav_pearson_r)
  )

message(sprintf("Running analysis on %d infected cells across %d images.",
                nrow(d_analysis), n_distinct(d_analysis$image_id)))


old_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))

# ─── FIT GLMM BINOMIAL ───
message("Fitting GLMM (binomial)...")
fit <- glmer(true_coloc ~ condition * antibody + (1 | image_id), data = d_analysis, family = binomial)
print(summary(fit))

# ─── ANOVA (TYPE III) ───
message("Computing Type III ANOVA...")
anova_df <- as.data.frame(car::Anova(fit, type = "III")) %>%
  rownames_to_column(var = "Term") %>%
  rename(P_value = `Pr(>Chisq)`)

# ─── PAIRWISE CONTRASTS VIA EMMEANS ───
message("Computing pairwise contrasts (V1 vs V2 stratified by antibody)...")
em <- emmeans(fit, ~ condition | antibody)
pairs_res <- contrast(em, method = "pairwise")
# Holm adjustment across the two contrasts (one per antibody)
pairs_summary <- as.data.frame(summary(pairs_res, by = NULL, adjust = "holm"))

if ("antibody" %in% colnames(pairs_summary)) {
  pairs_summary <- pairs_summary %>% rename(Antibody = antibody)
}

# ─── CLUSTER BOOTSTRAP (5000 RESAMPLES OF FIELDS OF VIEW) ───
BOOTSTRAP_SEED <- 20260524
set.seed(BOOTSTRAP_SEED)
n_boot <- 5000

boot_estimates <- matrix(NA_real_, nrow = n_boot, ncol = 2)
colnames(boot_estimates) <- c("Syn1", "EPR20535")

image_ids <- unique(d_analysis$image_id)

message(sprintf("Running cluster bootstrap (n = %d iterations)...", n_boot))
for (i in 1:n_boot) {
  if (i %% 1000 == 0) {
    message("  Completed ", i, " / ", n_boot, " iterations...")
  }

  # Resample fields of view with replacement; each draw becomes its own cluster
  sampled_ids <- sample(image_ids, size = length(image_ids), replace = TRUE)
  boot_df <- lapply(seq_along(sampled_ids), function(j) {
    id <- sampled_ids[j]
    clus_df <- d_analysis %>% filter(image_id == id)
    clus_df$image_id <- paste0(id, "_boot_", j)
    clus_df
  }) %>% bind_rows()

  ests <- tryCatch({
    suppressMessages(suppressWarnings({
      fit_boot <- glmer(true_coloc ~ condition * antibody + (1 | image_id), data = boot_df, family = binomial)
      em_boot <- emmeans(fit_boot, ~ condition | antibody)
      pairs_boot <- as.data.frame(contrast(em_boot, method = "pairwise"))

      if ("antibody" %in% colnames(pairs_boot)) {
        pairs_boot <- pairs_boot %>% rename(Antibody = antibody)
      }

      syn1_est <- pairs_boot$estimate[pairs_boot$Antibody == "Syn1"]
      epr_est  <- pairs_boot$estimate[pairs_boot$Antibody == "EPR20535"]

      res <- c(NA_real_, NA_real_)
      if (length(syn1_est) == 1) res[1] <- syn1_est
      if (length(epr_est) == 1) res[2] <- epr_est
      res
    }))
  }, error = function(e) NULL)

  if (!is.null(ests)) {
    boot_estimates[i, ] <- ests
  }
}

options(old_contrasts)

# Percentile CI and two-sided bootstrap p-value for one antibody
boot_summary <- function(est) {
  est <- est[!is.na(est)]
  if (length(est) == 0) return(list(lo = NA_real_, hi = NA_real_, p = NA_real_, n = 0L))
  p <- min(1, 2 * min(mean(est <= 0), mean(est >= 0)))
  list(lo = exp(quantile(est, 0.025)), hi = exp(quantile(est, 0.975)), p = p, n = length(est))
}
bs_syn1 <- boot_summary(boot_estimates[, "Syn1"])
bs_epr  <- boot_summary(boot_estimates[, "EPR20535"])

message(sprintf("Bootstrap success rates: Syn1: %d/%d (%.1f%%), EPR20535: %d/%d (%.1f%%)",
                bs_syn1$n, n_boot, 100 * bs_syn1$n / n_boot,
                bs_epr$n,  n_boot, 100 * bs_epr$n / n_boot))

# ─── COMBINE RESULTS FOR EXPORT ───
pairwise_df <- pairs_summary %>%
  mutate(
    Comparison = as.character(contrast),
    Boot_CI_Low_OR  = if_else(Antibody == "Syn1", bs_syn1$lo, bs_epr$lo),
    Boot_CI_High_OR = if_else(Antibody == "Syn1", bs_syn1$hi, bs_epr$hi),
    P_boot          = if_else(Antibody == "Syn1", bs_syn1$p,  bs_epr$p),
    OR = exp(estimate)
  ) %>%
  select(
    Comparison,
    Antibody,
    Estimate_logOR = estimate,
    OR,
    SE_logOR = SE,
    Z_ratio = z.ratio,
    P_GLMM = p.value,
    Boot_CI_Low_OR,
    Boot_CI_High_OR,
    P_boot
  )

# ─── EXPORT TO EXCEL ───
dir.create(dirname(out_xlsx), showWarnings = FALSE, recursive = TRUE)

wb <- createWorkbook()
addWorksheet(wb, "M3C_Freq_GLMM_ANOVA")
writeData(wb, "M3C_Freq_GLMM_ANOVA", anova_df)
addWorksheet(wb, "M3C_Freq_GLMM_Pairwise")
writeData(wb, "M3C_Freq_GLMM_Pairwise", pairwise_df)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Excel results exported to: ", out_xlsx)

# ─── CONSOLE SUMMARY ───
raw_summary <- d_analysis %>%
  group_by(antibody, condition) %>%
  summarize(
    Total_Cells = n(),
    Coloc_Cells = sum(true_coloc, na.rm = TRUE),
    Pct_Coloc = 100 * Coloc_Cells / Total_Cells,
    .groups = "drop"
  )

cat("\n======================================================================\n")
cat("                  M3C FREQUENCY GLMM ANALYSIS SUMMARY                 \n")
cat("======================================================================\n\n")

cat("--- Raw % True Colocalization (Costes p < 0.05) per Group ---\n")
print(as.data.frame(raw_summary))
cat("\n")

cat("--- GLMM Type III ANOVA (sum-to-zero contrasts) ---\n")
print(anova_df)
cat("\n")

cat("--- Pairwise Contrasts (V1 vs V2 stratified by Antibody, Holm-adjusted) ---\n")
for (i in seq_len(nrow(pairwise_df))) {
  row <- pairwise_df[i, ]
  cat(sprintf("Antibody: %s\n", row$Antibody))
  cat(sprintf("  Raw %% Coloc: V1 = %.1f%%, V2 = %.1f%%\n",
              raw_summary$Pct_Coloc[raw_summary$antibody == row$Antibody & raw_summary$condition == "V1"],
              raw_summary$Pct_Coloc[raw_summary$antibody == row$Antibody & raw_summary$condition == "V2"]))
  cat(sprintf("  GLMM Log Odds Ratio (V1 vs V2): %.4f (SE: %.4f, p-value: %.4f)\n",
              row$Estimate_logOR, row$SE_logOR, row$P_GLMM))
  cat(sprintf("  GLMM Odds Ratio: %.4f\n", row$OR))
  cat(sprintf("  Bootstrap OR 95%% CI: [%.4f, %.4f] (p-boot: %.4f, successful boot: %d/%d)\n",
              row$Boot_CI_Low_OR, row$Boot_CI_High_OR, row$P_boot,
              if (row$Antibody == "Syn1") bs_syn1$n else bs_epr$n, n_boot))
  cat("\n")
}
cat("======================================================================\n")
