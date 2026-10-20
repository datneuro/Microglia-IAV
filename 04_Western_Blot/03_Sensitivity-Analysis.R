# ==============================================================================
# 03_Sensitivity-Analysis.R — Prior and ROPE sensitivity analysis
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(emmeans)
  library(bayestestR)
  library(openxlsx)
})

source(file.path("04_Western_Blot", "00_Config-and-Ingestion.R"))

message("\n########## WB SENSITIVITY ANALYSIS ##########")

# 1. Define alternative priors (PRIORS_FIX is defined in 00_Config-and-Ingestion.R)
# PRIORS_FIX: class = b normal(0, 1.5), class = sigma student_t(3, 0, 2.5)
PRIORS_BROADER  <- c(prior(normal(0, 3), class = b),
                     prior(student_t(3, 0, 5), class = sigma))
PRIORS_NARROWER <- c(prior(normal(0, 1), class = b),
                     prior(student_t(3, 0, 1), class = sigma))

prior_sets <- list(
  default  = PRIORS_FIX,
  broader  = PRIORS_BROADER,
  narrower = PRIORS_NARROWER
)

# 2. Phase 1 Ingestion and Fitting
p1 <- ingest_phase1()
targets1 <- unique(p1$target)

p1_list <- list()
if (RUN_BRANCH_A) {
  p1_list[["A"]] <- list(branch = "A_model_membrane", data = p1, ycol = "y", with_membrane = TRUE)
}
if (RUN_BRANCH_B) {
  p1B <- add_branchB(p1, "y", "group", "C", match_vars = "target")
  p1_list[["B"]] <- list(branch = "B_within_membrane_FC", data = p1B, ycol = "y_fc", with_membrane = FALSE)
}

sens_results_p1 <- list()

for (name in names(p1_list)) {
  branch_tag <- p1_list[[name]]$branch
  ycol <- p1_list[[name]]$ycol
  with_mem <- p1_list[[name]]$with_membrane
  
  fixed_fmla <- if (with_mem) paste0(ycol, " ~ LPS*Virus + membrane") else paste0(ycol, " ~ LPS*Virus")
  message("  [Phase 1] Fitting ", branch_tag, " ...")
  
  for (pname in names(prior_sets)) {
    pr <- prior_sets[[pname]]
    message("    -> Prior: ", pname)
    
    # Compile ONE base model per branch/prior combination
    base_dat <- filter(p1_list[[name]]$data, target == targets1[1], is.finite(.data[[ycol]]))
    base_fit <- brm(as.formula(fixed_fmla),
                    data = base_dat, family = student(), prior = pr,
                    seed = BAYES_SEED, chains = MCMC_CHAINS, iter = MCMC_ITER_SENS,
                    warmup = MCMC_WARMUP, cores = MCMC_CORES, refresh = 0, silent = 2, backend = "rstan",
                    control = list(adapt_delta = 0.95))
    
    for (tg in targets1) {
      d <- filter(p1_list[[name]]$data, target == tg, is.finite(.data[[ycol]]))
      fit <- tryCatch(update(base_fit, newdata = d, refresh = 0, silent = 2),
                      error = function(e) { message("      update fail ", tg, ": ", e$message); NULL })
      if (is.null(fit)) next
      
      emm <- emmeans(fit, ~ LPS * Virus)
      vs_ctrl <- contrast(emm, method = CMAT1)
      bd_default <- bayes_describe(vs_ctrl, rope = ROPE_LOG2)
      
      bd_default$prior_set <- pname
      bd_default$response  <- tg
      bd_default$branch    <- branch_tag
      bd_default$phase     <- "P1"
      bd_default$timepoint <- NA_character_
      
      if (pname == "default") {
        bd_rope04 <- bayestestR::describe_posterior(vs_ctrl, ci = 0.95, ci_method = "hdi", test = "rope", rope_range = c(-0.4, 0.4), rope_ci = 1)
        bd_rope06 <- bayestestR::describe_posterior(vs_ctrl, ci = 0.95, ci_method = "hdi", test = "rope", rope_range = c(-0.6, 0.6), rope_ci = 1)
        
        bd_default$ROPE_pct_04 <- bd_rope04$ROPE_Percentage
        bd_default$ROPE_pct_06 <- bd_rope06$ROPE_Percentage
      } else {
        bd_default$ROPE_pct_04 <- NA_real_
        bd_default$ROPE_pct_06 <- NA_real_
      }
      
      sens_results_p1[[paste(branch_tag, tg, pname, sep="_")]] <- bd_default
    }
  }
}
res_p1 <- bind_rows(sens_results_p1)

# 3. Phase 2 Ingestion and Fitting
p2 <- ingest_phase2()
tc <- phase2_timecourse(p2)

# Prepare responses list for Branch A and B
yc_list <- list()
for (yc in names(PHASE2_RESPONSES)) {
  if (RUN_BRANCH_A) yc_list[[paste0(yc, "_A")]] <- list(ycol = yc, branch = "A_model_membrane", data = tc)
  if (RUN_BRANCH_B) {
    tcB <- add_branchB(tc, yc, "Treatment", "C", match_vars = "timepoint")
    yc_list[[paste0(yc, "_B")]] <- list(ycol = paste0(yc, "_fc"), branch = "B_within_membrane_FC", data = tcB)
  }
}

# Refit under each prior and extract contrasts
sens_results_p2 <- list()

for (name in names(yc_list)) {
  yc <- yc_list[[name]]$ycol
  branch_tag <- yc_list[[name]]$branch
  d <- filter(yc_list[[name]]$data, is.finite(.data[[yc]]))
  
  message("  [Phase 2] Fitting ", name, " ...")
  for (pname in names(prior_sets)) {
    pr <- prior_sets[[pname]]
    message("    -> Prior: ", pname)
    
    # We use MCMC_ITER_SENS for this analysis and run in parallel
    fit <- brm(as.formula(paste0(yc, " ~ timepoint*Treatment + membrane")),
               data = d, family = student(), prior = pr,
               seed = BAYES_SEED, chains = MCMC_CHAINS, iter = MCMC_ITER_SENS,
               warmup = MCMC_WARMUP, cores = MCMC_CORES, refresh = 0, silent = 2, backend = "rstan",
               control = list(adapt_delta = 0.95))
    
    # Extract Treatment vs Control within timepoint
    emm_tt <- emmeans(fit, ~ Treatment | timepoint)
    vs_ctrl <- contrast(emm_tt, method = "trt.vs.ctrl", ref = 1)
    
    # Describe posterior at default ROPE
    bd_default <- bayes_describe(vs_ctrl, rope = ROPE_LOG2)
    
    bd_default$prior_set <- pname
    bd_default$response  <- yc
    bd_default$branch    <- branch_tag
    bd_default$phase     <- "P2"
    
    # If default prior, also run ROPE=0.4 and ROPE=0.6 for stability check
    if (pname == "default") {
      bd_rope04 <- bayestestR::describe_posterior(vs_ctrl, ci = 0.95, ci_method = "hdi", test = "rope", rope_range = c(-0.4, 0.4), rope_ci = 1)
      bd_rope06 <- bayestestR::describe_posterior(vs_ctrl, ci = 0.95, ci_method = "hdi", test = "rope", rope_range = c(-0.6, 0.6), rope_ci = 1)
      
      bd_default$ROPE_pct_04 <- bd_rope04$ROPE_Percentage
      bd_default$ROPE_pct_06 <- bd_rope06$ROPE_Percentage
    } else {
      bd_default$ROPE_pct_04 <- NA_real_
      bd_default$ROPE_pct_06 <- NA_real_
    }
    
    sens_results_p2[[paste(name, pname, sep="_")]] <- bd_default
  }
}

res_p2 <- bind_rows(sens_results_p2)
res_df <- bind_rows(res_p1, res_p2)

# 4. Classify accept/reject/undecided and compute stability flags
classify_rope <- function(rope_pct, hdi_low, hdi_high, rope_val) {
  # Normalize rope_pct to [0, 1] if it happens to be 0-100
  rp <- ifelse(!is.na(rope_pct) & rope_pct > 1, rope_pct / 100, rope_pct)
  reject <- (hdi_low > rope_val) | (hdi_high < -rope_val)
  accept <- (!is.na(rp) & rp >= 0.95)
  case_when(
    reject ~ "reject",
    accept ~ "accept",
    TRUE ~ "undecided"
  )
}

res_df <- res_df %>%
  mutate(
    class_05 = classify_rope(ROPE_pct, HDI_low, HDI_high, 0.5),
    class_04 = if_else(prior_set == "default", classify_rope(ROPE_pct_04, HDI_low, HDI_high, 0.4), NA_character_),
    class_06 = if_else(prior_set == "default", classify_rope(ROPE_pct_06, HDI_low, HDI_high, 0.6), NA_character_)
  )

# Compute rope_stable: identical class across all three widths
rope_stability <- res_df %>%
  filter(prior_set == "default") %>%
  mutate(rope_stable = (class_04 == class_05) & (class_05 == class_06)) %>%
  select(phase, response, branch, timepoint, contrast, rope_stable)

# Compute prior_stable: BOTH median sign and class_05 unchanged across default/broader/narrower
prior_stability <- res_df %>%
  group_by(phase, response, branch, timepoint, contrast) %>%
  summarize(
    prior_stable = (n_distinct(sign(median)) == 1) & (n_distinct(class_05) == 1),
    .groups = "drop"
  )

res_final <- res_df %>%
  left_join(rope_stability, by = c("phase", "response", "branch", "timepoint", "contrast")) %>%
  left_join(prior_stability, by = c("phase", "response", "branch", "timepoint", "contrast")) %>%
  arrange(phase, response, branch, timepoint, contrast, factor(prior_set, levels=c("default", "broader", "narrower"))) %>%
  select(phase, response, branch, timepoint, contrast, prior_set, median, HDI_low, HDI_high, pd, ROPE_pct, class_05, rope_stable, prior_stable, everything())

write_csv(res_final, file.path(OUT_DIR, "WB_sensitivity.csv"))
message("  Saved WB_sensitivity.csv")

# 5. Excel Export
excel_path <- file.path(OUT_DIR, "WB_Stats_Report.xlsx")
if (file.exists(excel_path)) {
  wb <- loadWorkbook(excel_path)
  sheet_name <- "9_Sensitivity"
  if (sheet_name %in% names(wb)) {
    removeWorksheet(wb, sheet_name)
  }
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, res_final)
  
  hdr <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2")
  addStyle(wb, sheet_name, hdr, rows = 1, cols = seq_len(ncol(res_final)), gridExpand = TRUE)
  
  sig_both <- createStyle(fgFill = "#E2EFDA") # both stable, Structural green
  
  # Highlight rows stable on BOTH flags
  stable_idx <- which(res_final$rope_stable == TRUE & res_final$prior_stable == TRUE)
  if (length(stable_idx) > 0) {
    addStyle(wb, sheet_name, sig_both, rows = stable_idx + 1L, cols = seq_len(ncol(res_final)), gridExpand = TRUE, stack = TRUE)
  }
  
  saveWorkbook(wb, excel_path, overwrite = TRUE)
  message("  Appended 9_Sensitivity to WB_Stats_Report.xlsx")
} else {
  message("  WB_Stats_Report.xlsx not found, run 02_Bayesian-Densitometry.R first.")
}

message("\n[03_Sensitivity-Analysis] DONE.")
