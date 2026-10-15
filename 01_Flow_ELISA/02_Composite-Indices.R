################################################################################
# 02_Composite-Indices.R — Composite scores + Decoupling / Polarization indices
# Control-referenced effect size → composites → indices → lm (HC3) + cluster
# bootstrap (BCa) 95% CIs; leave-one-marker-out sensitivity.
#   e = (value_t − mean of Control at that timepoint) / SD of all samples at that timepoint
#   Decoupling Index   = Secreted M1 − Surface M1   (negative = surface-dominant)
#   Polarization Index = mean(all M1 markers) − M2
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(car)
  library(emmeans)
  library(sandwich)
  library(boot)
})

set.seed(2026)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir <- "Flow-ELISA"  # Flow/ELISA data folder (contains DATA/; receives R_analysis/)
out_dir <- file.path(base_dir, "R_analysis", "output")

# ── Load & transform ──────────────────────────────────────────────────────────
df_fact <- read.csv(file.path(out_dir, "tidy_factorial.csv"),
                    stringsAsFactors = FALSE) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c(4, 8, 24)),
    LPS       = factor(LPS,       levels = c("No", "Yes")),
    Virus     = factor(Virus,     levels = c("None", "V1", "V2")),
    Bio_Rep   = factor(Bio_Rep),
    Value_t   = ifelse(Platform == "flow",
                       log(Value / (100 - Value)),
                       log2(Value))
  )

# ── Marker sets (exclude CD68, NP, CCL2 from composites) ─────────────────────
surface_M1  <- c("CD40", "CD86", "HLA-I", "HLA-II")
secreted_M1 <- c("TNF-a", "IL-1b", "IL-6", "IL-12p70")
M2_markers  <- c("CD206", "IL-10")
all_M1      <- c(surface_M1, secreted_M1)
composite_markers <- c(all_M1, M2_markers)

cat("Surface M1:", paste(surface_M1, collapse = ", "), "\n")
cat("Secreted M1:", paste(secreted_M1, collapse = ", "), "\n")
cat("M2:", paste(M2_markers, collapse = ", "), "\n")
cat("Total composite markers:", length(composite_markers), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Control-referenced effect size
# ══════════════════════════════════════════════════════════════════════════════
# e_i = (value_t_i − mean(value_t[Control, tp])) / SD_all(tp)
# SD_all = SD of all samples of the 6 groups at that tp (one sd() over the
# whole timepoint, i.e. including between-group variance; not a pooled within-group SD)

# Step 1a: Get control means per marker × tp
ctrl_means <- df_fact %>%
  filter(Group == "C") %>%
  group_by(Marker, Timepoint) %>%
  summarise(ctrl_mean = mean(Value_t), .groups = "drop")

# Step 1b: SD of all samples per marker × tp (across all 6 groups)
all_sd <- df_fact %>%
  group_by(Marker, Timepoint) %>%
  summarise(sd_all = sd(Value_t), .groups = "drop")

# Step 1c: Compute effect size
df_eff <- df_fact %>%
  left_join(ctrl_means, by = c("Marker", "Timepoint")) %>%
  left_join(all_sd,     by = c("Marker", "Timepoint")) %>%
  mutate(e = (Value_t - ctrl_mean) / sd_all)

# Check
cat("── Effect size diagnostics ────────────────────────\n")
ctrl_check <- df_eff %>%
  filter(Group == "C") %>%
  group_by(Marker, Timepoint) %>%
  summarise(mean_e = mean(e), sd_e = sd(e), .groups = "drop")
cat("Control group mean effect (should be ~0):",
    round(range(ctrl_check$mean_e), 6), "\n")
cat("Control group SD effect:", round(range(ctrl_check$sd_e), 3), "\n\n")

# Save effect size table
write.csv(df_eff %>% select(Group, Timepoint, LPS, Virus, Bio_Rep, Platform,
                            Marker, Value, Value_t, e, Unit),
          file.path(out_dir, "L2_effect_sizes.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 2. Composites (per Bio_Rep × Group × Timepoint)
# ══════════════════════════════════════════════════════════════════════════════
# Keep only composite markers
df_comp <- df_eff %>% filter(Marker %in% composite_markers)

# Surface M1 composite
surface_comp <- df_comp %>%
  filter(Marker %in% surface_M1) %>%
  group_by(Group, Timepoint, LPS, Virus, Bio_Rep) %>%
  summarise(Surface_M1 = mean(e), .groups = "drop")

# Secreted M1 composite
secreted_comp <- df_comp %>%
  filter(Marker %in% secreted_M1) %>%
  group_by(Group, Timepoint, LPS, Virus, Bio_Rep) %>%
  summarise(Secreted_M1 = mean(e), .groups = "drop")

# M2 composite
m2_comp <- df_comp %>%
  filter(Marker %in% M2_markers) %>%
  group_by(Group, Timepoint, LPS, Virus, Bio_Rep) %>%
  summarise(M2 = mean(e), .groups = "drop")

# All M1
all_m1_comp <- df_comp %>%
  filter(Marker %in% all_M1) %>%
  group_by(Group, Timepoint, LPS, Virus, Bio_Rep) %>%
  summarise(all_M1 = mean(e), .groups = "drop")

# Merge
comp_df <- surface_comp %>%
  left_join(secreted_comp, by = c("Group","Timepoint","LPS","Virus","Bio_Rep")) %>%
  left_join(m2_comp,       by = c("Group","Timepoint","LPS","Virus","Bio_Rep")) %>%
  left_join(all_m1_comp,   by = c("Group","Timepoint","LPS","Virus","Bio_Rep"))

# ── Indices ────────────────────────────────────────────────────────────────────
comp_df <- comp_df %>%
  mutate(
    Decoupling_Index    = Secreted_M1 - Surface_M1,
    Polarization_Index  = all_M1 - M2
  )

cat("── Composites & Indices ───────────────────────────\n")
cat("Composite df dim:", nrow(comp_df), "x", ncol(comp_df), "\n")
cat("Groups:", paste(unique(comp_df$Group), collapse = ", "), "\n\n")

# Summary by group × tp
comp_summary <- comp_df %>%
  group_by(Group, Timepoint) %>%
  summarise(
    across(c(Surface_M1, Secreted_M1, M2, Decoupling_Index, Polarization_Index),
           list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

cat("── Composite means (selected) ────────────────────\n")
for (i in seq_len(nrow(comp_summary))) {
  r <- comp_summary[i,]
  cat(sprintf("  %-8s tp=%-3s  DI=%+.3f(%.3f)  PI=%+.3f(%.3f)  Surf=%+.3f  Sec=%+.3f  M2=%+.3f\n",
              r$Group, r$Timepoint,
              r$Decoupling_Index_mean, r$Decoupling_Index_sd,
              r$Polarization_Index_mean, r$Polarization_Index_sd,
              r$Surface_M1_mean, r$Secreted_M1_mean, r$M2_mean))
}
cat("\n")

# Save composites
write.csv(comp_df, file.path(out_dir, "L2_composites.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 3. Model each index (lm + HC3 + 7 contrasts)
# ══════════════════════════════════════════════════════════════════════════════
old_opts <- options(contrasts = c("contr.sum", "contr.poly"))

# Robust planned contrasts: build weight vectors from the ACTUAL emmeans grid.
# emmeans orders the LPS:Virus grid LPS-FASTEST (No.None, Yes.None, No.V1,
# Yes.V1, No.V2, Yes.V2), NOT Virus-fastest. Matching cells by (LPS,Virus)
# label avoids the silent mis-comparison caused by positional weight vectors.
build_contrasts_from_emm <- function(emm) {
  g   <- as.data.frame(summary(emm))
  tp1 <- as.character(g$Timepoint)[1]
  sub <- g[as.character(g$Timepoint) == tp1, c("LPS", "Virus")]
  key <- paste(as.character(sub$LPS), as.character(sub$Virus), sep = ".")
  pos <- function(lps, vir) {
    i <- which(key == paste(lps, vir, sep = "."))
    if (length(i) != 1)
      stop("emmeans cell not found / not unique: ", lps, ".", vir,
           " (grid = ", paste(key, collapse = ", "), ")")
    i
  }
  mk <- function(plus, minus) {
    v <- numeric(length(key)); v[plus] <- 1; v[minus] <- -1; v
  }
  list(
    "Ctrl-LPS"      = mk(pos("No",  "None"), pos("Yes", "None")),
    "Ctrl-V1"       = mk(pos("No",  "None"), pos("No",  "V1")),
    "Ctrl-V2"       = mk(pos("No",  "None"), pos("No",  "V2")),
    "V1-V2"         = mk(pos("No",  "V1"),   pos("No",  "V2")),
    "LPS-LPS.V1"    = mk(pos("Yes", "None"), pos("Yes", "V1")),
    "LPS-LPS.V2"    = mk(pos("Yes", "None"), pos("Yes", "V2")),
    "LPS.V1-LPS.V2" = mk(pos("Yes", "V1"),   pos("Yes", "V2"))
  )
}

analyze_index <- function(dat, index_col, index_name) {
  dat$y <- dat[[index_col]]

  mod <- lm(y ~ Timepoint * LPS * Virus + Timepoint:Bio_Rep, data = dat)

  # ANOVA
  anova_std <- car::Anova(mod, type = 3)
  anova_hc3 <- car::Anova(mod, type = 3, white.adjust = "hc3")

  all_terms  <- rownames(anova_std)
  terms_keep <- setdiff(all_terms, c("(Intercept)","Residuals","Timepoint:Bio_Rep"))
  ss_resid   <- anova_std["Residuals", "Sum Sq"]

  anova_out <- data.frame(
    Index        = index_name,
    Term         = terms_keep,
    SS           = anova_std[terms_keep, "Sum Sq"],
    partial_eta2 = anova_std[terms_keep, "Sum Sq"] /
                     (anova_std[terms_keep, "Sum Sq"] + ss_resid),
    F_hc3        = anova_hc3[terms_keep, "F"],
    p_hc3        = anova_hc3[terms_keep, "Pr(>F)"],
    row.names    = NULL,
    stringsAsFactors = FALSE
  )

  # Contrasts
  vcov_hc3 <- sandwich::vcovHC(mod, type = "HC3")
  emm      <- emmeans::emmeans(mod, ~ LPS:Virus | Timepoint, vcov. = vcov_hc3)
  con_s    <- as.data.frame(summary(
    emmeans::contrast(emm, method = build_contrasts_from_emm(emm)),
    infer = c(TRUE, TRUE), level = 0.95
  ))

  con_out <- data.frame(
    Index     = index_name,
    Contrast  = con_s$contrast,
    Timepoint = as.character(con_s$Timepoint),
    estimate  = con_s$estimate,
    SE        = con_s$SE,
    CI_low    = con_s$lower.CL,
    CI_high   = con_s$upper.CL,
    p_raw     = con_s$p.value,
    family    = paste0("index_", index_name),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  list(anova = anova_out, contrasts = con_out, model = mod)
}

# ── Worked example: Decoupling Index ─────────────────────────────────────────
cat("═══════════════════════════════════════════════════\n")
cat("   WORKED EXAMPLE — Decoupling Index (lm + HC3)\n")
cat("═══════════════════════════════════════════════════\n\n")

pilot_di <- analyze_index(comp_df, "Decoupling_Index", "Decoupling")
cat("── ANOVA ──\n")
print(pilot_di$anova, digits = 4, row.names = FALSE)
cat("\n── Contrasts (pre-Holm) ──\n")
print(pilot_di$contrasts, digits = 4, row.names = FALSE)
cat("\nResidual df:", df.residual(pilot_di$model), "\n")
cat("R²:", round(summary(pilot_di$model)$r.squared, 4), "\n\n")

cat("═══ Running all indices ═══\n\n")

# ── All indices + composites ─────────────────────────────────────────────────
index_names <- c("Surface_M1", "Secreted_M1", "M2",
                 "Decoupling_Index", "Polarization_Index")

all_anova <- list()
all_con   <- list()

for (idx in index_names) {
  cat("  Analyzing:", idx, "... ")
  res <- analyze_index(comp_df, idx, idx)
  all_anova[[idx]] <- res$anova
  all_con[[idx]]   <- res$contrasts
  cat("done\n")
}

idx_anova_df <- bind_rows(all_anova)
idx_con_df   <- bind_rows(all_con)

# Holm per index (each index = 1 family)
idx_con_df <- idx_con_df %>%
  group_by(Index) %>%
  mutate(p_holm = p.adjust(p_raw, method = "holm")) %>%
  ungroup()

# ══════════════════════════════════════════════════════════════════════════════
# 4. Cluster bootstrap (Bio_Rep) → 95% CI for indices by group × tp
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Cluster bootstrap (5000 reps) ──────────────────\n")

boot_index_ci <- function(comp_data, index_col, n_boot = 5000) {
  # cluster = Bio_Rep; resample Bio_Reps with replacement within each group×tp
  results <- list()

  for (grp in unique(comp_data$Group)) {
    for (tp in levels(comp_data$Timepoint)) {
      subset_data <- comp_data %>%
        filter(Group == grp, Timepoint == tp)

      bio_reps <- unique(subset_data$Bio_Rep)
      n_reps   <- length(bio_reps)

      boot_fn <- function(data, indices) {
        # Resample Bio_Reps
        sampled_reps <- bio_reps[indices]
        boot_data <- data %>%
          filter(Bio_Rep %in% sampled_reps)
        # Handle resampled duplicates by row binding
        boot_vals <- numeric(length(indices))
        for (j in seq_along(indices)) {
          boot_vals[j] <- data[[index_col]][data$Bio_Rep == sampled_reps[j]]
        }
        mean(boot_vals)
      }

      b <- boot::boot(subset_data, boot_fn, R = n_boot)
      ci <- boot::boot.ci(b, conf = 0.95, type = "bca")

      results[[length(results) + 1]] <- data.frame(
        Index     = index_col,
        Group     = grp,
        Timepoint = tp,
        mean      = mean(subset_data[[index_col]]),
        boot_mean = mean(b$t),
        CI_low    = ci$bca[4],
        CI_high   = ci$bca[5],
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(results)
}

boot_results <- list()
for (idx in c("Decoupling_Index", "Polarization_Index",
              "Surface_M1", "Secreted_M1", "M2")) {
  cat("  Bootstrap:", idx, "... ")
  boot_results[[idx]] <- boot_index_ci(comp_df, idx, n_boot = 5000)
  cat("done\n")
}
boot_df <- bind_rows(boot_results)

cat("\n── Bootstrap CIs for Decoupling Index ─────────────\n")
di_boot <- boot_df %>% filter(Index == "Decoupling_Index") %>% arrange(Group, Timepoint)
for (i in seq_len(nrow(di_boot))) {
  r <- di_boot[i,]
  cat(sprintf("  %-8s tp=%-3s  mean=%+.3f  95%%CI=[%+.3f, %+.3f]\n",
              r$Group, r$Timepoint, r$mean, r$CI_low, r$CI_high))
}
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Sensitivity: leave-one-marker-out for Decoupling Index (Virus effect)
# ══════════════════════════════════════════════════════════════════════════════
cat("── Sensitivity: leave-one-marker-out (Decoupling) ─\n")

all_di_markers <- c(surface_M1, secreted_M1)  # 8 markers in DI
sens_results <- list()

for (drop_marker in all_di_markers) {
  # Recompute composites without this marker
  surf_m <- setdiff(surface_M1, drop_marker)
  sec_m  <- setdiff(secreted_M1, drop_marker)

  loo_surf <- df_eff %>%
    filter(Marker %in% surf_m) %>%
    group_by(Group, Timepoint, LPS, Virus, Bio_Rep) %>%
    summarise(Surface_M1 = mean(e), .groups = "drop")

  loo_sec <- df_eff %>%
    filter(Marker %in% sec_m) %>%
    group_by(Group, Timepoint, LPS, Virus, Bio_Rep) %>%
    summarise(Secreted_M1 = mean(e), .groups = "drop")

  loo_df <- loo_surf %>%
    left_join(loo_sec, by = c("Group","Timepoint","LPS","Virus","Bio_Rep")) %>%
    mutate(DI_loo = Secreted_M1 - Surface_M1)

  # Test Virus effect on DI_loo
  mod_loo <- lm(DI_loo ~ Timepoint * LPS * Virus + Timepoint:Bio_Rep, data = loo_df)
  a_hc3   <- car::Anova(mod_loo, type = 3, white.adjust = "hc3")
  a_std   <- car::Anova(mod_loo, type = 3)
  ss_r    <- a_std["Residuals", "Sum Sq"]

  # Extract Virus term
  virus_f <- a_hc3["Virus", "F"]
  virus_p <- a_hc3["Virus", "Pr(>F)"]
  virus_eta2 <- a_std["Virus", "Sum Sq"] / (a_std["Virus", "Sum Sq"] + ss_r)

  sens_results[[length(sens_results) + 1]] <- data.frame(
    dropped_marker = drop_marker,
    Virus_F        = virus_f,
    Virus_p        = virus_p,
    Virus_eta2     = virus_eta2,
    stringsAsFactors = FALSE
  )
}

sens_df <- bind_rows(sens_results)

# Full model (no marker dropped)
mod_full <- lm(Decoupling_Index ~ Timepoint * LPS * Virus + Timepoint:Bio_Rep,
               data = comp_df)
a_full_hc3 <- car::Anova(mod_full, type = 3, white.adjust = "hc3")
a_full_std <- car::Anova(mod_full, type = 3)
ss_r_full  <- a_full_std["Residuals", "Sum Sq"]

full_row <- data.frame(
  dropped_marker = "(none — full)",
  Virus_F        = a_full_hc3["Virus", "F"],
  Virus_p        = a_full_hc3["Virus", "Pr(>F)"],
  Virus_eta2     = a_full_std["Virus", "Sum Sq"] /
                     (a_full_std["Virus", "Sum Sq"] + ss_r_full),
  stringsAsFactors = FALSE
)

sens_df <- bind_rows(full_row, sens_df)
cat("Sensitivity analysis — Virus effect on Decoupling Index:\n")
for (i in seq_len(nrow(sens_df))) {
  r <- sens_df[i,]
  sig_flag <- ifelse(r$Virus_p < 0.05, "***", ifelse(r$Virus_p < 0.1, ".", "ns"))
  cat(sprintf("  Drop %-12s → F=%.2f  p=%.4f  η²=%.4f  %s\n",
              r$dropped_marker, r$Virus_F, r$Virus_p, r$Virus_eta2, sig_flag))
}
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Save all outputs
# ══════════════════════════════════════════════════════════════════════════════
write.csv(idx_anova_df, file.path(out_dir, "L2_index_anova.csv"),      row.names = FALSE)
write.csv(idx_con_df,   file.path(out_dir, "L2_index_contrasts.csv"),  row.names = FALSE)
write.csv(boot_df,      file.path(out_dir, "L2_bootstrap_CI.csv"),     row.names = FALSE)
write.csv(sens_df,      file.path(out_dir, "L2_sensitivity_LOO.csv"),  row.names = FALSE)

options(old_opts)

# ── Final summary ─────────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════\n")
cat("         COMPOSITE INDICES — SUMMARY\n")
cat("═══════════════════════════════════════════════════\n\n")

cat("Index ANOVA:", nrow(idx_anova_df), "rows\n")
cat("Index contrasts:", nrow(idx_con_df), "rows\n")
cat("Bootstrap CIs:", nrow(boot_df), "rows\n")
cat("Sensitivity LOO:", nrow(sens_df), "rows\n\n")

cat("── Significant index contrasts (p_holm < 0.05) ───\n")
sig_idx <- idx_con_df %>%
  filter(p_holm < 0.05) %>%
  arrange(Index, Timepoint)
for (i in seq_len(nrow(sig_idx))) {
  r <- sig_idx[i,]
  cat(sprintf("  %-20s %-16s tp=%-3s est=%+.3f [%+.3f,%+.3f] p_holm=%.2e\n",
              r$Index, r$Contrast, r$Timepoint, r$estimate,
              r$CI_low, r$CI_high, r$p_holm))
}

cat("\n── Files saved ───────────────────────────────────\n")
cat("  L2_effect_sizes.csv\n")
cat("  L2_composites.csv\n")
cat("  L2_index_anova.csv\n")
cat("  L2_index_contrasts.csv\n")
cat("  L2_bootstrap_CI.csv\n")
cat("  L2_sensitivity_LOO.csv\n")
cat("\n═══ 02_Composite-Indices COMPLETE ═══\n")
