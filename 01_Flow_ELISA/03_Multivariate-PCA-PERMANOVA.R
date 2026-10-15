################################################################################
# 03_Multivariate-PCA-PERMANOVA.R — Correlation PCA + PERMANOVA
# Correlation PCA on control-referenced effects, PERMANOVA, betadisper
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(vegan)
})

set.seed(2026)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir <- "Flow-ELISA"  # Flow/ELISA data folder (contains DATA/; receives R_analysis/)
out_dir <- file.path(base_dir, "R_analysis", "output")

# ── Load effect sizes ─────────────────────────────────────────────────────────
df_eff <- read.csv(file.path(out_dir, "L2_effect_sizes.csv"),
                   stringsAsFactors = FALSE) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c(4, 8, 24)),
    LPS       = factor(LPS,       levels = c("No", "Yes")),
    Virus     = factor(Virus,     levels = c("None", "V1", "V2")),
    Bio_Rep   = factor(Bio_Rep)
  )

# ── Build wide matrix: rows=samples, cols=markers (effect e) ─────────────────
pca_markers <- c("CD40", "CD68", "CD86", "CD206", "HLA-I", "HLA-II",
                 "CCL2", "IL-1b", "IL-6", "IL-10", "TNF-a", "IL-12p70")

df_wide <- df_eff %>%
  filter(Marker %in% pca_markers) %>%
  select(Group, Timepoint, LPS, Virus, Bio_Rep, Marker, e) %>%
  pivot_wider(names_from = Marker, values_from = e)

cat("PCA matrix dim:", nrow(df_wide), "x", length(pca_markers), "markers\n")
cat("Markers:", paste(pca_markers, collapse = ", "), "\n\n")

# Extract numeric matrix
mat <- as.matrix(df_wide[, pca_markers])
rownames(mat) <- paste(df_wide$Group, df_wide$Timepoint, df_wide$Bio_Rep, sep = "_")

cat("── Check: CD86 effect in matrix ──────────────────\n")
cat("CD86 range:", round(range(mat[, "CD86"]), 3), "\n")
cat("CD86 mean:", round(mean(mat[, "CD86"]), 4), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 1. PCA (correlation PCA: scale = TRUE)
# ══════════════════════════════════════════════════════════════════════════════
pca_res <- prcomp(mat, scale. = TRUE, center = TRUE)

# Variance explained
var_exp <- summary(pca_res)$importance
cat("── PCA Variance Explained ────────────────────────\n")
cat(sprintf("  PC1: %.1f%% | PC2: %.1f%% | PC3: %.1f%% | PC4: %.1f%%\n",
            var_exp[2,1]*100, var_exp[2,2]*100, var_exp[2,3]*100, var_exp[2,4]*100))
cat(sprintf("  Cumulative PC1-4: %.1f%%\n\n", var_exp[3,4]*100))

# Loadings
loadings_df <- as.data.frame(pca_res$rotation[, 1:4])
loadings_df$Marker <- rownames(loadings_df)
loadings_df <- loadings_df %>% select(Marker, PC1, PC2, PC3, PC4)

cat("── PC1 Loadings (sorted) ─────────────────────────\n")
l1 <- loadings_df %>% arrange(desc(abs(PC1)))
for (i in seq_len(nrow(l1))) {
  cat(sprintf("  %-12s %+.3f\n", l1$Marker[i], l1$PC1[i]))
}

cat("\n── PC2 Loadings (sorted) ─────────────────────────\n")
l2 <- loadings_df %>% arrange(desc(abs(PC2)))
for (i in seq_len(nrow(l2))) {
  cat(sprintf("  %-12s %+.3f\n", l2$Marker[i], l2$PC2[i]))
}
cat("\n")

# Describe PC directions (DO NOT name axes)
cat("── PC Direction Description ──────────────────────\n")
cat("PC1: ", ifelse(sum(loadings_df$PC1 > 0) > sum(loadings_df$PC1 < 0),
    "Predominantly positive loadings across markers — general activation magnitude",
    "Mixed loading pattern"), "\n")
cat("PC2: Positive loaders vs negative loaders — see loadings above\n\n")

# Scores with metadata
scores_df <- as.data.frame(pca_res$x[, 1:4])
scores_df <- cbind(df_wide[, c("Group","Timepoint","LPS","Virus","Bio_Rep")],
                   scores_df)

# ══════════════════════════════════════════════════════════════════════════════
# 2. PERMANOVA (vegan::adonis2)
# ══════════════════════════════════════════════════════════════════════════════
cat("── PERMANOVA (Euclidean, 9999 permutations) ──────\n")

dist_mat <- dist(mat, method = "euclidean")

perm_res <- vegan::adonis2(
  dist_mat ~ Timepoint * LPS * Virus,
  data = df_wide,
  method = "euclidean",
  permutations = 9999,
  by = "terms"
)

cat("\n")
print(perm_res)
cat("\n")

# Extract results
perm_terms <- rownames(perm_res)
perm_terms <- setdiff(perm_terms, c("Residual", "Total"))

permanova_df <- data.frame(
  Term    = perm_terms,
  Df      = perm_res[perm_terms, "Df"],
  SumOfSqs = perm_res[perm_terms, "SumOfSqs"],
  R2      = perm_res[perm_terms, "R2"],
  F_value = perm_res[perm_terms, "F"],
  p_value = perm_res[perm_terms, "Pr(>F)"],
  row.names = NULL,
  stringsAsFactors = FALSE
)

cat("── R² by factor ──────────────────────────────────\n")
for (i in seq_len(nrow(permanova_df))) {
  r <- permanova_df[i,]
  sig <- ifelse(r$p_value < 0.001, "***",
         ifelse(r$p_value < 0.01, "**",
         ifelse(r$p_value < 0.05, "*", "ns")))
  cat(sprintf("  %-25s R²=%.4f  F=%.2f  p=%.4f  %s\n",
              r$Term, r$R2, r$F_value, r$p_value, sig))
}
cat(sprintf("  Residual                    R²=%.4f\n",
            perm_res["Residual", "R2"]))
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. betadisper + permutest (for each factor)
# ══════════════════════════════════════════════════════════════════════════════
cat("── betadisper + permutest ─────────────────────────\n")

betadisper_results <- list()
factors_to_test <- c("Timepoint", "LPS", "Virus")

for (fac in factors_to_test) {
  groups <- df_wide[[fac]]
  bd <- vegan::betadisper(dist_mat, groups)
  pt <- vegan::permutest(bd, permutations = 9999)
  p_val <- pt$tab[1, "Pr(>F)"]
  f_val <- pt$tab[1, "F"]

  flag <- ifelse(p_val < 0.05, "⚠️  DISPERSION CONFOUND", "OK — homogeneous")
  cat(sprintf("  %-12s F=%.2f  p=%.4f  → %s\n", fac, f_val, p_val, flag))

  betadisper_results[[length(betadisper_results) + 1]] <- data.frame(
    Factor           = fac,
    F_betadisper     = f_val,
    p_betadisper     = p_val,
    dispersion_flag  = ifelse(p_val < 0.05, "CONFOUND", "OK"),
    stringsAsFactors = FALSE
  )
}
cat("\n")

betadisper_df <- bind_rows(betadisper_results)

# ══════════════════════════════════════════════════════════════════════════════
# 4. Save outputs
# ══════════════════════════════════════════════════════════════════════════════
write.csv(loadings_df,  file.path(out_dir, "L3_pca_loadings.csv"),    row.names = FALSE)
write.csv(scores_df,    file.path(out_dir, "L3_pca_scores.csv"),      row.names = FALSE)
write.csv(permanova_df, file.path(out_dir, "L3_permanova.csv"),       row.names = FALSE)
write.csv(betadisper_df,file.path(out_dir, "L3_betadisper.csv"),      row.names = FALSE)

# Scree data
scree_df <- data.frame(
  PC      = paste0("PC", 1:ncol(pca_res$rotation)),
  Eigenvalue = pca_res$sdev^2,
  Var_pct = var_exp[2,] * 100,
  Cum_pct = var_exp[3,] * 100
)
write.csv(scree_df, file.path(out_dir, "L3_scree.csv"), row.names = FALSE)

cat("── Files saved ───────────────────────────────────\n")
cat("  L3_pca_loadings.csv\n")
cat("  L3_pca_scores.csv\n")
cat("  L3_permanova.csv\n")
cat("  L3_betadisper.csv\n")
cat("  L3_scree.csv\n")
cat("\n═══ 03_Multivariate-PCA-PERMANOVA COMPLETE ═══\n")
