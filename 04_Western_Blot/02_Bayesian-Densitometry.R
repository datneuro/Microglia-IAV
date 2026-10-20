# ==============================================================================
# 02_Bayesian-Densitometry.R  
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(emmeans)
  library(bayestestR)
  library(car)
  library(sandwich)
  library(openxlsx)
  library(posterior)
})
if (!exists("ingest_phase1")) source(file.path("04_Western_Blot", "00_Config-and-Ingestion.R"))
emm_options(rg.limit = 20000)


EFF_CMAT <- list(
  "LPS main"          = c(-1, 1, -1,  1, -1,  1) / 3,  # avg(+LPS) - avg(-LPS)
  "Virus main: PR8"   = c(-1,-1,  1,  1,  0,  0) / 2,  # avg(PR8) - avg(None)
  "Virus main: S33"   = c(-1,-1,  0,  0,  1,  1) / 2,
  "Virus×LPS: PR8"    = c( 1,-1, -1,  1,  0,  0),       # (LPS/PR8-PR8)-(LPS-C)
  "Virus×LPS: S33"    = c( 1,-1,  0,  0, -1,  1)
)
compute_effects <- function(fit) {
  emm <- emmeans(fit, ~ LPS * Virus)
  ef  <- bayes_describe(contrast(emm, method = EFF_CMAT))
  ef$effect <- ef$contrast
  ef
}

freq_hc3 <- function(fit_lm, emm_spec, cmat) {
  V <- sandwich::vcovHC(fit_lm, type = "HC3")
  emm <- emmeans(fit_lm, emm_spec, vcov. = V)
  ct  <- as.data.frame(summary(contrast(emm, method = cmat, adjust = "none"),
                               infer = c(TRUE, TRUE)))
  ct$p_holm <- p.adjust(ct$p.value, "holm")
  ct
}

# ==============================================================================
# §1  PHASE 1  — Bayesian factorial per target
# ==============================================================================
message("\n########## PHASE 1: Bayesian factorial (24h screen) ##########")
p1 <- ingest_phase1()
targets1 <- unique(p1$target)

# CMAT1 is defined in 00_Config-and-Ingestion.R

run_phase1_branch <- function(ycol, with_membrane, branch_tag) {
  fixed <- if (with_membrane) paste0(ycol, " ~ LPS*Virus + membrane")
           else               paste0(ycol, " ~ LPS*Virus")
  message("  [", branch_tag, "] compiling base model: ", fixed)
  base_dat <- filter(p1, target == targets1[1], is.finite(.data[[ycol]]))
  base_fit <- BRM(as.formula(fixed), base_dat, PRIORS_FIX)

  res <- list(); freq <- list(); effs <- list(); qc_list <- list()
  for (tg in targets1) {
    d <- filter(p1, target == tg, is.finite(.data[[ycol]]))
    fit <- tryCatch(update(base_fit, newdata = d, refresh = 0, silent = 2),
                    error = function(e) { message("    update fail ", tg, ": ", e$message); NULL })
    if (is.null(fit)) next
    
    # Convergence check (mcmc_qc)
    qc_list[[tg]] <- mcmc_qc(fit, paste0(tg, " (", branch_tag, ")"))
    
    emm <- emmeans(fit, ~ LPS * Virus)
    bd  <- bayes_describe(contrast(emm, method = CMAT1))
    bd$target <- tg; bd$branch <- branch_tag
    res[[tg]] <- bd
    if (branch_tag == "A_model_membrane") {
      ef <- tryCatch(compute_effects(fit), error = function(e) NULL)
      if (!is.null(ef) && nrow(ef)) { ef$target <- tg; effs[[tg]] <- ef }
    }
    # frequentist HC3
    lmf <- tryCatch(lm(as.formula(fixed), data = d), error = function(e) NULL)
    if (!is.null(lmf)) {
      fh <- tryCatch(freq_hc3(lmf, ~ LPS * Virus, CMAT1), error = function(e) NULL)
      if (!is.null(fh)) { fh$target <- tg; fh$branch <- branch_tag; freq[[tg]] <- fh }
    }
  }
  list(bayes = bind_rows(res), freq = bind_rows(freq), effects = bind_rows(effs), qc = bind_rows(qc_list))
}

phase1_out <- list()
if (RUN_BRANCH_A) phase1_out$A <- run_phase1_branch("y",    TRUE,  "A_model_membrane")
if (RUN_BRANCH_B) {
  p1 <- add_branchB(p1, "y", "group", "C", match_vars = "target")
  phase1_out$B <- run_phase1_branch("y_fc", FALSE, "B_within_membrane_FC")
}
phase1_bayes   <- bind_rows(lapply(phase1_out, `[[`, "bayes"))
phase1_freq    <- bind_rows(lapply(phase1_out, `[[`, "freq"))
phase1_effects <- bind_rows(lapply(phase1_out, `[[`, "effects"))
write_csv(phase1_bayes,   file.path(OUT_DIR, "WB_phase1_bayes_contrasts.csv"))
write_csv(phase1_freq,    file.path(OUT_DIR, "WB_phase1_freq_hc3_contrasts.csv"))
write_csv(phase1_effects, file.path(OUT_DIR, "WB_phase1_factor_effects.csv"))
message("  Phase1 done: ", nrow(phase1_bayes), " bayes rows.")

# Group summaries (raw log2 scale, branch A response)
phase1_summary <- p1 %>%
  group_by(target, group, LPS, Virus) %>%
  summarize(n = n(), mean_log2 = mean(y), sd_log2 = sd(y),
            sem_log2 = sd(y)/sqrt(n()), .groups = "drop")
write_csv(phase1_summary, file.path(OUT_DIR, "WB_phase1_group_summary.csv"))

# ==============================================================================
# §2  PHASE 2  — Bayesian fixed-block timecourse
# ==============================================================================
message("\n########## PHASE 2: Bayesian fixed-block (aSyn timecourse) ##########")
p2 <- ingest_phase2()
tc <- phase2_timecourse(p2)
pn <- phase2_panel24(p2)

run_phase2_timecourse <- function(ycol, branch_tag, dat) {
  message("  [", branch_tag, "] ", ycol, " ~ Timepoint*Treatment + membrane")
  d <- filter(dat, is.finite(.data[[ycol]]))
  fit <- BRM(as.formula(paste0(ycol, " ~ timepoint*Treatment + membrane")), d, PRIORS_FIX)

  # Convergence check (mcmc_qc)
  qc_check <- mcmc_qc(fit, paste0(ycol, " (", branch_tag, ")"))
  message("    [MCMC QC Check] Divergent: ", qc_check$divergent, 
          " | Max Rhat: ", round(qc_check$max_rhat, 3), 
          " | Min Tail ESS: ", round(qc_check$min_ess_tail, 0))

  # (a) Treatment vs Control within each timepoint
  emm_tt <- emmeans(fit, ~ Treatment | timepoint)
  vs_ctrl <- bayes_describe(contrast(emm_tt, method = "trt.vs.ctrl", ref = 1))
  vs_ctrl$kind <- "Treatment_vs_Control"

  # (b) polynomial time trend per treatment (linear + quadratic)
  emm_tp <- emmeans(fit, ~ timepoint | Treatment)
  poly <- bayes_describe(contrast(emm_tp, method = "poly", max.degree = 2))
  poly$kind <- "Polynomial_time"

  out <- bind_rows(vs_ctrl, poly)
  out$response <- ycol; out$branch <- branch_tag
  list(stats = out, qc = qc_check)
}

phase2_rows <- list()
for (yc in names(PHASE2_RESPONSES)) {
  if (RUN_BRANCH_A) phase2_rows[[paste0(yc,"_A")]] <- run_phase2_timecourse(yc, "A_model_membrane", tc)
  if (RUN_BRANCH_B) {
    tcB <- add_branchB(tc, yc, "Treatment", "C", match_vars = "timepoint")
    phase2_rows[[paste0(yc,"_B")]] <- run_phase2_timecourse(paste0(yc,"_fc"), "B_within_membrane_FC", tcB)
  }
}
phase2_temporal <- bind_rows(lapply(phase2_rows, `[[`, "stats"))
phase2_qc       <- bind_rows(lapply(phase2_rows, `[[`, "qc"))
write_csv(phase2_temporal, file.path(OUT_DIR, "WB_phase2_temporal.csv"))
write_csv(phase2_qc, file.path(OUT_DIR, "WB_phase2_qc.csv"))
message("  Phase2 timecourse done: ", nrow(phase2_temporal), " rows.")

# ---- 24h transfection-sensing panel {C,T,PO,TPO} ----
message("  [panel24] aSyn responses ~ Treatment + membrane, 24h only")
CMAT_PANEL <- list(
  "TPO vs T (intracellular dsRNA)" = c( 0, -1, 0, 1),   # order C,T,PO,TPO
  "TPO vs PO (intra vs extra)"     = c( 0,  0,-1, 1),
  "PO vs C (extracellular PolyIC)" = c(-1,  0, 1, 0),
  "T vs C (transfection vehicle)"  = c(-1,  1, 0, 0)
)
panel_rows <- list(); panel_qc <- list()
for (yc in names(PHASE2_RESPONSES)) {
  d <- filter(pn, is.finite(.data[[yc]]))
  if (nrow(d) < 8) next
  fit <- tryCatch(BRM(as.formula(paste0(yc, " ~ Treatment + membrane")), d, PRIORS_FIX),
                  error = function(e) { message("    panel fit fail: ", e$message); NULL })
  if (is.null(fit)) next

  # Convergence check (mcmc_qc)
  panel_qc[[yc]] <- mcmc_qc(fit, paste0(yc, "_Panel24"))
    
  emm <- emmeans(fit, ~ Treatment)
  bd  <- bayes_describe(contrast(emm, method = CMAT_PANEL))
  bd$response <- yc; bd$branch <- "A_model_membrane"
  panel_rows[[yc]] <- bd
}
phase2_panel <- bind_rows(panel_rows)
write_csv(phase2_panel, file.path(OUT_DIR, "WB_phase2_panel24.csv"))

phase2_summary <- p2 %>%
  pivot_longer(all_of(names(PHASE2_RESPONSES)), names_to = "response", values_to = "val") %>%
  filter(is.finite(val)) %>%
  group_by(response, timepoint, group) %>%
  summarize(n = n(), mean_log2 = mean(val), sd_log2 = sd(val),
            sem_log2 = sd(val)/sqrt(n()), .groups = "drop")
write_csv(phase2_summary, file.path(OUT_DIR, "WB_phase2_group_summary.csv"))

# ==============================================================================
# §3  EXCEL EXPORT
# ==============================================================================
message("\n########## EXCEL EXPORT ##########")

wb_qc_master <- bind_rows(
  phase1_out$A$qc,   
  phase1_out$B$qc,   
  phase2_qc,
  bind_rows(panel_qc) 
)

# Convergence criteria for the green "pass" highlight
wb_qc_master$pass <- with(wb_qc_master,
  divergent == 0 & max_rhat <= 1.01 & min_ess_bulk >= 400 & min_ess_tail >= 400)

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2")
sig <- createStyle(fgFill = "#E2EFDA")   # Structural green pass style
add_sheet <- function(name, df, pd_col = "pd") {
  addWorksheet(wb, name); writeData(wb, name, df)
  addStyle(wb, name, hdr, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
  if (pd_col %in% colnames(df)) {
    r <- which(!is.na(df[[pd_col]]) & df[[pd_col]] >= 0.975) + 1L
    if (length(r)) addStyle(wb, name, sig, rows = r, cols = seq_len(ncol(df)),
                            gridExpand = TRUE, stack = TRUE)
  }
}

methods <- tibble(
  Section = c("Design P1","Design P2","Response","Membrane P1","Membrane P2",
              "Model P1","Model P2","Priors","Branches","Bayes metrics",
              "Freq cross-check","Temporal","24h panel","DATA WARNING","Software"),
  Description = c(
    "2(LPS:No/Yes) x 3(Virus:None/PR8/S33), n=3, 24h",
    "Treatment{C,PO,PR8,S33} x {4h,8h,24h}, n=4; + 24h panel {C,T,PO,TPO}",
    "log2(target/actin); pathway ratios log2(pIkB/IkB),log2(pP65/P65),log2(LC3-II/LC3-I); aSyn: syn1,syn1.25,psyn all /actin",
    "mem column (2 membranes); LC3/CTSD override: rep1 | rep2+rep3",
    "rep = membrane (4 membranes)",
    "brm(y ~ LPS*Virus (+membrane), family=student)",
    "brm(y ~ Timepoint*Treatment + membrane, family=student)",
    "b~N(0,1.5) [log2]; sigma~student_t(3,0,2.5)  (weakly-informative fixed block)",
    paste0("A=model membrane (run=", RUN_BRANCH_A, "); B=within-membrane FC (run=", RUN_BRANCH_B, ")"),
    "posterior median, 95% HDI, pd, ROPE% (|0.5| log2); ROPE = 0.5 log2 (~1.4-fold) is the smallest biologically meaningful WB difference; BF10 is not computed.",
    "lm + HC3 robust SE + emmeans contrasts; Holm per target family",
    "emmeans poly contrasts (linear+quadratic) over 3 timepoints",
    "TPO-T (intracellular dsRNA), TPO-PO (intra vs extra), PO-C, T-C",
    "Phase2 'psyn' corrected to real densitometry (2026-06-14); y_pocc INCLUDED",
    paste0("R ", getRversion(), "; brms/rstan, emmeans, bayestestR, sandwich, openxlsx")),
)

addWorksheet(wb, "0_Methods"); writeData(wb, "0_Methods", methods)
addStyle(wb, "0_Methods", hdr, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "0_Methods", cols = 1:2, widths = c(20, 120))

add_sheet("1_P1_Bayes_Contrasts", phase1_bayes)
add_sheet("2_P1_Factor_Effects",  phase1_effects)
add_sheet("3_P1_Freq_HC3",        phase1_freq, pd_col = "p_holm")
add_sheet("4_P2_Temporal",        phase2_temporal)
add_sheet("5_P2_Panel24",         phase2_panel)
add_sheet("6_P1_Group_Summary",   phase1_summary)
add_sheet("7_P2_Group_Summary",   phase2_summary)

# Sheet 8 conditional formatting: highlight model fits that pass all criteria
add_sheet("8_MCMC_QC_Summary",    wb_qc_master) # Run without pd_col matching rule
ok <- which(wb_qc_master$pass) + 1L
if (length(ok)) {
  addStyle(wb, "8_MCMC_QC_Summary", sig, rows = ok, cols = seq_len(ncol(wb_qc_master)),
           gridExpand = TRUE, stack = TRUE)
}

saveWorkbook(wb, file.path(OUT_DIR, "WB_Stats_Report.xlsx"), overwrite = TRUE)
message("  Saved WB_Stats_Report.xlsx")

# Save fitted objects digest for figure script
saveRDS(list(phase1_bayes = phase1_bayes, phase1_effects = phase1_effects,
             phase1_summary = phase1_summary, phase2_temporal = phase2_temporal,
             phase2_panel = phase2_panel, phase2_summary = phase2_summary),
        file.path(OUT_DIR, "WB_results_digest.rds"))
message("\n[02_Bayesian-Densitometry] DONE.")