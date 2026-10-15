################################################################################
# 01_Univariate-Factorial.R — Univariate factorial analysis per marker
# Transform (flow % → logit; ELISA pg/mL → log2)
#   → lm(Value_t ~ Timepoint * LPS * Virus + Timepoint:Bio_Rep)
#   → Type III ANOVA with HC3 F/p + planned emmeans contrasts (HC3), Holm per marker
# Poly(I:C) is compared with its reference groups by Welch t-tests.
# CD86 is printed in full first as a worked example, then all markers are run.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(car)
  library(emmeans)
  library(sandwich)
})

set.seed(2026)

# --Path setup -----
base_dir <- "Flow-ELISA"  # Flow/ELISA data folder (contains DATA/; receives R_analysis/)
out_dir <- file.path(base_dir, "R_analysis", "output")

# --Load data -----
df_fact <- read.csv(file.path(out_dir, "tidy_factorial.csv"),
                    stringsAsFactors = FALSE)
df_poly <- read.csv(file.path(out_dir, "tidy_polyic.csv"),
                    stringsAsFactors = FALSE)

df_fact <- df_fact %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c(4, 8, 24)),
    LPS       = factor(LPS,       levels = c("No", "Yes")),
    Virus     = factor(Virus,     levels = c("None", "V1", "V2")),
    Bio_Rep   = factor(Bio_Rep)
  )

df_poly <- df_poly %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c(4, 8, 24)),
    Bio_Rep   = factor(Bio_Rep)
  )

# --Transform ----
# flow % → logit (natural log);  ELISA pg/mL → log2
df_fact <- df_fact %>%
  mutate(Value_t = ifelse(Platform == "flow",
                          log(Value / (100 - Value)),   # logit
                          log2(Value)))                 # log2

df_poly <- df_poly %>%
  mutate(Value_t = ifelse(Platform == "flow",
                          log(Value / (100 - Value)),
                          log2(Value)))

old_opts <- options(contrasts = c("contr.sum", "contr.poly"))

# --Planed contrasts -----
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

# ── Reference groups for Poly(I:C) Welch tests ────────────────────────────────
ref_groups <- list(
  "PolyIC-Ctrl" = list(LPS = "No",  Virus = "None"),
  "PolyIC-LPS"  = list(LPS = "Yes", Virus = "None"),
  "PolyIC-V1"   = list(LPS = "No",  Virus = "V1"),
  "PolyIC-V2"   = list(LPS = "No",  Virus = "V2")
)

# ══════════════════════════════════════════════════════════════════════════════
# Helper: analyze one marker
# ══════════════════════════════════════════════════════════════════════════════
analyze_marker <- function(marker_name, dat_f, dat_p) {

  dat <- dat_f %>% filter(Marker == marker_name)

  # ---- Fit model ----
  mod <- lm(Value_t ~ Timepoint * LPS * Virus + Timepoint:Bio_Rep, data = dat)

  # ---- ANOVA ----
  # Standard Type III → SS for partial η²
  anova_std <- car::Anova(mod, type = 3)
  # HC3-adjusted → F, p for inference
  anova_hc3 <- car::Anova(mod, type = 3, white.adjust = "hc3")

  # Terms of interest (exclude intercept, residuals, block term)
  all_terms  <- rownames(anova_std)
  terms_keep <- setdiff(all_terms,
                        c("(Intercept)", "Residuals", "Timepoint:Bio_Rep"))

  ss_resid <- anova_std["Residuals", "Sum Sq"]

  anova_out <- data.frame(
    Marker       = marker_name,
    Term         = terms_keep,
    SS           = anova_std[terms_keep, "Sum Sq"],
    partial_eta2 = anova_std[terms_keep, "Sum Sq"] /
                     (anova_std[terms_keep, "Sum Sq"] + ss_resid),
    F_hc3        = anova_hc3[terms_keep, "F"],
    p_hc3        = anova_hc3[terms_keep, "Pr(>F)"],
    row.names    = NULL,
    stringsAsFactors = FALSE
  )

  # ---- emmeans contrasts (HC3 vcov) ----
  vcov_hc3 <- sandwich::vcovHC(mod, type = "HC3")
  emm      <- emmeans::emmeans(mod, ~ LPS:Virus | Timepoint, vcov. = vcov_hc3)
  con      <- emmeans::contrast(emm, method = build_contrasts_from_emm(emm))
  con_s    <- as.data.frame(summary(con, infer = c(TRUE, TRUE), level = 0.95))

  con_out <- data.frame(
    Marker    = marker_name,
    Contrast  = con_s$contrast,
    Timepoint = as.character(con_s$Timepoint),
    estimate  = con_s$estimate,
    SE        = con_s$SE,
    CI_low    = con_s$lower.CL,
    CI_high   = con_s$upper.CL,
    p_raw     = con_s$p.value,
    family    = "factorial",
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # ---- Poly(I:C) Welch t-tests ----
  dat_p_m  <- dat_p %>% filter(Marker == marker_name)
  poly_rows <- vector("list", 12)
  idx <- 1
  for (tp_val in levels(dat$Timepoint)) {
    poly_vals <- dat_p_m %>% filter(Timepoint == tp_val) %>% pull(Value_t)
    for (cname in names(ref_groups)) {
      ref <- ref_groups[[cname]]
      ref_vals <- dat %>%
        filter(Timepoint == tp_val,
               LPS   == ref$LPS,
               Virus  == ref$Virus) %>%
        pull(Value_t)
      tt <- t.test(poly_vals, ref_vals)  # Welch (default)
      poly_rows[[idx]] <- data.frame(
        Marker    = marker_name,
        Contrast  = cname,
        Timepoint = tp_val,
        estimate  = unname(tt$estimate[1] - tt$estimate[2]),
        SE        = tt$stderr,
        CI_low    = tt$conf.int[1],
        CI_high   = tt$conf.int[2],
        p_raw     = tt$p.value,
        family    = "polyic",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  poly_out <- bind_rows(poly_rows)

  list(anova = anova_out, contrasts = con_out, polyic = poly_out,
       model = mod, emm = emm)
}

# ══════════════════════════════════════════════════════════════════════════════
# WORKED EXAMPLE — CD86
# ══════════════════════════════════════════════════════════════════════════════
cat("═══════════════════════════════════════════════════\n")
cat("       WORKED EXAMPLE — CD86\n")
cat("═══════════════════════════════════════════════════\n\n")

pilot <- analyze_marker("CD86", df_fact, df_poly)

cat("── Model diagnostics ─────────────────────────────\n")
cat("Residual df:", df.residual(pilot$model), "\n")
cat("R²:", round(summary(pilot$model)$r.squared, 4), "\n\n")

cat("── emmeans grid (verify LPS×Virus ordering) ──────\n")
print(summary(pilot$emm))
cat("\n")

cat("── ANOVA table (Type III SS + HC3 F/p) ──────────\n")
print(pilot$anova, digits = 4, row.names = FALSE)
cat("\n")

cat("── Factorial contrasts (HC3, pre-Holm) ───────────\n")
print(pilot$contrasts, digits = 4, row.names = FALSE)
cat("\n")

cat("── Poly(I:C) Welch t-tests (pre-Holm) ────────────\n")
print(pilot$polyic, digits = 4, row.names = FALSE)
cat("\n")

cat("═══ Running all markers ═══════════════════════════\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# FULL LOOP — All markers
# ══════════════════════════════════════════════════════════════════════════════
markers <- sort(unique(df_fact$Marker))
cat("Markers:", paste(markers, collapse = ", "), "\n\n")

all_anova     <- vector("list", length(markers))
all_contrasts <- vector("list", length(markers))

for (i in seq_along(markers)) {
  m <- markers[i]
  cat("  [", i, "/", length(markers), "]", m, "... ")
  res <- analyze_marker(m, df_fact, df_poly)
  all_anova[[i]]     <- res$anova
  all_contrasts[[i]] <- bind_rows(res$contrasts, res$polyic)
  cat("done\n")
}

anova_df    <- bind_rows(all_anova)
contrast_df <- bind_rows(all_contrasts)

# ── Holm correction: per marker × per family ──────────────────────────────────
contrast_df <- contrast_df %>%
  group_by(Marker, family) %>%
  mutate(p_holm = p.adjust(p_raw, method = "holm")) %>%
  ungroup()

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(anova_df,    file.path(out_dir, "L1_anova.csv"),     row.names = FALSE)
write.csv(contrast_df, file.path(out_dir, "L1_contrasts.csv"), row.names = FALSE)

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════\n")
cat("         UNIVARIATE ANALYSIS — SUMMARY\n")
cat("═══════════════════════════════════════════════════\n\n")

cat("ANOVA table:", nrow(anova_df), "rows (",
    length(markers), "markers ×",
    nrow(anova_df) / length(markers), "terms)\n")
cat("Contrast table:", nrow(contrast_df), "rows\n")
cat("  Factorial:", sum(contrast_df$family == "factorial"), "\n")
cat("  Poly(I:C):", sum(contrast_df$family == "polyic"), "\n\n")

# Significant ANOVA terms
cat("── Significant ANOVA terms (p_hc3 < 0.05) ───────\n")
sig_anova <- anova_df %>%
  filter(p_hc3 < 0.05) %>%
  arrange(Marker, p_hc3) %>%
  select(Marker, Term, partial_eta2, F_hc3, p_hc3)
print(as.data.frame(sig_anova), digits = 4)
cat("\n")

# Significant contrasts after Holm
cat("── Significant contrasts (p_holm < 0.05) ─────────\n")
sig_con <- contrast_df %>%
  filter(p_holm < 0.05) %>%
  arrange(Marker, family, Timepoint, Contrast) %>%
  select(Marker, Contrast, Timepoint, estimate, CI_low, CI_high, p_raw, p_holm, family)
cat(nrow(sig_con), "significant contrasts (out of", nrow(contrast_df), "total)\n\n")
print(as.data.frame(sig_con))

# Restore contrasts
options(old_opts)

cat("\n── Files saved ───────────────────────────────────\n")
cat("  L1_anova.csv     :", nrow(anova_df), "rows\n")
cat("  L1_contrasts.csv :", nrow(contrast_df), "rows\n")
cat("\n═══ 01_Univariate-Factorial COMPLETE ═══\n")

