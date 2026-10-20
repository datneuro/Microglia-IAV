# ==============================================================================
# 01_Membrane-Diagnostic.R  
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(car)      # leveneTest
})
`%||%` <- function(a, b) if (is.null(a)) b else a
if (!exists("ingest_phase1")) source(file.path("04_Western_Blot", "00_Config-and-Ingestion.R"))

# ------------------------------------------------------------------ helpers ---
safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)

membrane_bg <- function(df, ycol, group_col, ctrl) {
  # spread of per-membrane control means + per-membrane SDs
  ctl <- df %>% filter(.data[[group_col]] == ctrl)
  cm <- ctl %>% group_by(membrane) %>%
    summarize(m = mean(.data[[ycol]], na.rm = TRUE), .groups = "drop")
  sds <- df %>% group_by(membrane) %>%
    summarize(s = sd(.data[[ycol]], na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(s), s > 0)
  list(
    bg_range_ctrl = if (nrow(cm) >= 2) max(cm$m) - min(cm$m) else NA_real_,
    n_membrane    = n_distinct(df$membrane),
    scale_ratio   = if (nrow(sds) >= 2) max(sds$s) / min(sds$s) else NA_real_
  )
}

var_het <- function(resid, membrane) {
  d <- tibble(r = resid, m = factor(membrane)) %>% filter(is.finite(r))
  if (n_distinct(d$m) < 2) return(NA_real_)
  safe(as.numeric(car::leveneTest(r ~ m, data = d)[1, "Pr(>F)"]))
}

message("\n=== DIAGNOSTIC: PHASE 1 (24h screen) ===")
p1 <- ingest_phase1()
targets1 <- unique(p1$target)

diag1 <- map_dfr(targets1, function(tg) {
  d <- filter(p1, target == tg)
  if (n_distinct(d$membrane) < 2 || nrow(d) < 8) {
    return(tibble(target = tg, n_membrane = n_distinct(d$membrane),
                  bg_range_ctrl = NA, scale_ratio = NA, icc_membrane = NA,
                  AB_estimate_r = NA, AB_mean_absdiff = NA, var_het_p = NA,
                  note = "single membrane / too few"))
  }
  bg <- membrane_bg(d, "y", "group", "C")

  # Branch A: lm with membrane fixed (2 levels -> fixed, not random)
  fitA <- safe(lm(y ~ LPS * Virus + membrane, data = d))
  residA <- if (inherits(fitA, "lm")) residuals(fitA) else rep(NA, nrow(d))
  vh <- var_het(residA, d$membrane[as.integer(names(residA))] %||% d$membrane)

  # Branch B: within-membrane control referencing (all 24h -> no match var)
  dB <- add_branchB(d, "y", "group", "C")
  fitB <- safe(lm(y_fc ~ LPS * Virus, data = dB))

  # Compare 6 cell means (LPS*Virus) effect vs grand: use coef of LPS/Virus terms
  cmp <- tibble(AB_estimate_r = NA_real_, AB_mean_absdiff = NA_real_)
  if (inherits(fitA, "lm") && inherits(fitB, "lm")) {
    # predicted group means (relative to control) for the 6 groups
    grid <- distinct(d, LPS, Virus)
    pa <- safe(predict(fitA, newdata = grid %>% mutate(membrane = d$membrane[1])))
    pb <- safe(predict(fitB, newdata = grid))
    if (length(pa) == nrow(grid) && length(pb) == nrow(grid)) {
      pa <- pa - pa[1]; pb <- pb - pb[1]      # reference to control cell
      cmp <- tibble(AB_estimate_r   = safe(cor(pa, pb)),
                    AB_mean_absdiff = mean(abs(pa - pb), na.rm = TRUE))
    }
  }

  # ICC proxy: membrane SS share from lm
  icc <- NA_real_
  if (inherits(fitA, "lm")) {
    av <- safe(as.data.frame(car::Anova(fitA, type = 2)))
    if (is.data.frame(av) && "membrane" %in% rownames(av)) {
      ss <- av[["Sum Sq"]]; names(ss) <- rownames(av)
      icc <- as.numeric(ss["membrane"] / sum(ss, na.rm = TRUE))
    }
  }

  tibble(target = tg, n_membrane = bg$n_membrane,
         bg_range_ctrl = bg$bg_range_ctrl, scale_ratio = bg$scale_ratio,
         icc_membrane = icc,
         AB_estimate_r = cmp$AB_estimate_r, AB_mean_absdiff = cmp$AB_mean_absdiff,
         var_het_p = vh, note = "")
})

write_csv(diag1, file.path(OUT_DIR, "diag_phase1_membrane.csv"))
message("  wrote diag_phase1_membrane.csv  (", nrow(diag1), " targets)")


p2 <- ingest_phase2()
tc <- phase2_timecourse(p2)

resp2 <- PHASE2_RESPONSES   # y_pocc excluded: psyn is a bogus SUM formula

diag2 <- map_dfr(names(resp2), function(yc) {
  d <- tc %>% filter(is.finite(.data[[yc]]))
  bg <- membrane_bg(d, yc, "Treatment", "C")

  # Branch A: random membrane intercept (4 reps)
  fitA <- safe(lmer(reformulate(c("timepoint*Treatment", "(1|membrane)"), response = yc),
                    data = d, REML = TRUE,
                    control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))))
  icc <- NA_real_; residA <- rep(NA, nrow(d)); sd_mem <- NA_real_; sd_res <- NA_real_
  if (inherits(fitA, "merMod")) {
    vc <- as.data.frame(VarCorr(fitA))
    sd_mem <- vc$sdcor[vc$grp == "membrane"]
    sd_res <- vc$sdcor[vc$grp == "Residual"]
    icc <- sd_mem^2 / (sd_mem^2 + sd_res^2)
    residA <- residuals(fitA)
  }
  vh <- var_het(residA, d$membrane)

  # Branch B: within-membrane, within-timepoint control referencing
  dB <- add_branchB(d, yc, "Treatment", "C", match_vars = "timepoint")
  ycb <- paste0(yc, "_fc")
  fitB <- safe(lm(reformulate("timepoint*Treatment", response = ycb), data = dB))

  # Compare Treatment-vs-Control effects per timepoint (A emmeans-style vs B means)
  cmp <- tibble(AB_estimate_r = NA_real_, AB_mean_absdiff = NA_real_)
  if (inherits(fitA, "merMod") && inherits(fitB, "lm")) {
    grid <- distinct(d, timepoint, Treatment)
    pa <- safe(predict(fitA, newdata = grid, re.form = NA))
    pb <- safe(predict(fitB, newdata = grid))
    gg <- grid %>% mutate(pa = pa, pb = pb) %>%
      group_by(timepoint) %>%
      mutate(pa = pa - pa[Treatment == "C"], pb = pb - pb[Treatment == "C"]) %>%
      ungroup() %>% filter(Treatment != "C")
    cmp <- tibble(AB_estimate_r   = safe(cor(gg$pa, gg$pb)),
                  AB_mean_absdiff = mean(abs(gg$pa - gg$pb), na.rm = TRUE))
  }

  tibble(response = yc, formula = unname(resp2[yc]),
         n_membrane = bg$n_membrane, bg_range_ctrl = bg$bg_range_ctrl,
         scale_ratio = bg$scale_ratio, sd_membrane = sd_mem, sd_resid = sd_res,
         icc_membrane = icc, AB_estimate_r = cmp$AB_estimate_r,
         AB_mean_absdiff = cmp$AB_mean_absdiff, var_het_p = vh)
})

write_csv(diag2, file.path(OUT_DIR, "diag_phase2_membrane.csv"))
message("  wrote diag_phase2_membrane.csv  (", nrow(diag2), " responses)")

# ==============================================================================
# §3  CONSOLE SUMMARY + HEURISTIC DECISION
# ==============================================================================
cat("\n================ PHASE 1 DIAGNOSTIC ================\n")
print(as.data.frame(diag1 %>% mutate(across(where(is.numeric), ~round(.,3)))), row.names = FALSE)
cat("\n================ PHASE 2 DIAGNOSTIC ================\n")
print(as.data.frame(diag2 %>% mutate(across(where(is.numeric), ~round(.,3)))), row.names = FALSE)

# NOTE: under additive log2 membrane + (near-)balanced design, Branch A and B give
# IDENTICAL point estimates (AB_estimate_r=1, AB_mean_absdiff=0). The real question
# is therefore about UNCERTAINTY: if residual variance is membrane-dependent
# (var_het_p<0.05) or membrane scales differ a lot (scale_ratio high), a plain
# homoscedastic fit (either branch) under-states SE -> need robust SE (HC3) /
# Bayesian per-observation variance. Branch B additionally protects against large
# baseline shifts (bg_range) when the additive assumption is shaky.
decide <- function(bg_range, ab_r, ab_diff, vh, scale_ratio) {
  heterosced <- isTRUE(vh < 0.05) || isTRUE(scale_ratio > 2.0)
  baseline_shift <- isTRUE(bg_range > 1.0)        # >1 log2 (2x) gap between membranes
  est_diverge <- isTRUE(ab_r < 0.95) || isTRUE(ab_diff > 0.3)
  if (est_diverge && baseline_shift) "PREFER_B"
  else if (heterosced) "ROBUST_SE"                # estimates fine; need HC3 / het-var
  else "A_OK"
}

cat("\n--- HEURISTIC (per target/response) ---\n")
d1 <- diag1 %>% mutate(decision = pmap_chr(list(bg_range_ctrl, AB_estimate_r, AB_mean_absdiff, var_het_p, scale_ratio), decide))
d2 <- diag2 %>% mutate(decision = pmap_chr(list(bg_range_ctrl, AB_estimate_r, AB_mean_absdiff, var_het_p, scale_ratio), decide))
cat("Phase1:\n"); print(d1 %>% count(decision))
cat("Phase1 targets needing attention:\n")
print(d1 %>% filter(decision != "A_OK") %>% select(target, decision, bg_range_ctrl, scale_ratio, var_het_p))
cat("Phase2:\n"); print(d2 %>% select(response, decision, bg_range_ctrl, scale_ratio, var_het_p))

write_csv(d1, file.path(OUT_DIR, "diag_phase1_membrane.csv"))
write_csv(d2, file.path(OUT_DIR, "diag_phase2_membrane.csv"))
message("\n[01_Membrane-Diagnostic] done.")
