# ==============================================================================
# 00_Config-and-Ingestion.R 
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(brms)
  library(bayestestR)
  library(posterior)
})

# ==============================================================================
# §1  GLOBAL CONFIGURATION
# ==============================================================================

PROJ_DIR <- "."
WB_DIR   <- file.path(PROJ_DIR, "Western blot")
DATA_XLSX<- file.path(WB_DIR, "DATA", "WB_FORR.xlsx")
R_DIR    <- file.path(WB_DIR, "R_analysis")
OUT_DIR  <- file.path(R_DIR, "output")
FIG_DIR  <- file.path(R_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

RUN_BRANCH_A <- TRUE   
RUN_BRANCH_B <- TRUE   

# ---- Bayesian / MCMC settings ----
BAYES_SEED   <- 20260614
MCMC_ITER    <- 25000
MCMC_WARMUP  <- 2000
MCMC_CHAINS  <- 8
MCMC_CORES   <- 8
MCMC_ITER_SENS <- 10000
ROPE_LOG2    <- 0.5    
BAYES_DO_BF  <- FALSE 

STRAIN_V1 <- "A/PR/8/34"   # PR8
STRAIN_V2 <- "A/NWS/33"    # S33

# ==============================================================================
# §2  COLOR PALETTE 
# ==============================================================================

PAL_TREAT <- c(
  "C"       = "#D5D8DC",   # control light gray
  "PR8"     = "#A8C7F0",   # V1 light blue
  "S33"     = "#F2B4B4",   # V2 light red
  "LPS"     = "#6E7376",   # dark gray
  "LPS/PR8" = "#4285F4",   # dark blue
  "LPS/S33" = "#C0392B",   # dark red
  "PO"      = "#2ECC71",   # PolyIC treatment green
  "TPO"     = "#16A085",   # PolyIC transfected teal-green
  "T"       = "#95A5A6"    # transfection vehicle gray
)

p_to_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", "ns"))))
}
fmt_p  <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# ==============================================================================
# §2b  OFFICIAL TARGET LABELS 
# ==============================================================================
TARGET_LABELS <- c(
  "syn1.mono"       = "α-Synuclein monomer (Syn1)",
  "syn1.multi"      = "α-Synuclein multimer (Syn1)",
  "syn1.25"         = "α-Synuclein 25 kDa (Syn1)",
  "syn211"          = "α-Synuclein (Syn211)",
  "psyn"            = "phospho-α-Synuclein",
  "psyn_25kda"      = "phospho-α-Synuclein 25 kDa",
  "psmb5"           = "PSMB5",
  "psmb5.pre"       = "PSMB5 precursor",
  "k48"             = "K48-polyubiquitin",
  "psmb8"           = "PSMB8",
  "asc"             = "ASC",
  "pikb"            = "phospho-IκBα",
  "ikb"             = "IκBα (total)",
  "pP65"            = "phospho-NF-κB p65",
  "P65"             = "NF-κB p65 (total)",
  "lc3.1"           = "LC3-I",
  "lc3.2"           = "LC3-II",
  "pre.CTSD"        = "Cathepsin D precursor",
  "CTSD"            = "Cathepsin D (mature)",
  # derived phospho/total + LC3 ratios
  "ratio.pIkB.IkB"  = "p-IκBα / IκBα ratio",
  "ratio.pP65.P65"  = "p-p65 / p65 ratio",
  "ratio.LC3"       = "LC3-II / LC3-I ratio"
)
lab_target <- function(x) ifelse(x %in% names(TARGET_LABELS), TARGET_LABELS[x], x)

# ==============================================================================
# §3  PHASE 1 INGESTION  (24h target screen)
# ==============================================================================

#   pre.CTSD = procathepsin-D (>35 kDa), CTSD = mature cathepsin-D (>25 kDa)
PHASE1_TARGETS_OVERRIDE <- c("lc3.1", "lc3.2", "pre.CTSD", "CTSD")

ingest_phase1 <- function(xlsx = DATA_XLSX) {
  raw <- read_excel(xlsx, sheet = "Phase1_24h-Only")
  meta_cols   <- c("id", "group", "timepoint", "mem", "actin")
  target_cols <- setdiff(names(raw), meta_cols)

  long <- raw %>%
    pivot_longer(all_of(target_cols), names_to = "target", values_to = "signal") %>%
    mutate(
      rep      = as.integer(sub(".*-(\\d+)$", "\\1", id)),
      actin    = as.numeric(actin),
      signal   = as.numeric(signal),
      y        = log2(signal / actin),                       # branch-A response
      # factorial decode
      LPS      = factor(if_else(grepl("LPS", group), "LPS", "No LPS"),
                        levels = c("No LPS", "LPS")),
      Virus    = factor(case_when(grepl("PR8", group) ~ "PR8",
                                  grepl("S33", group) ~ "S33",
                                  TRUE ~ "None"),
                        levels = c("None", "PR8", "S33")),
      group    = factor(group,
                        levels = c("C","PR8","S33","LPS","LPS/PR8","LPS/S33")),
      # per-target membrane id
      membrane = if_else(
        target %in% PHASE1_TARGETS_OVERRIDE,
        if_else(rep == 1L, "M_r1", "M_r23"),     # LC3 / CTSD layout
        paste0("M", mem)                          # default layout
      )
    ) %>%
    filter(is.finite(y)) %>%
    mutate(membrane = factor(membrane))

  # ---- derived ratio targets (phospho/total, LC3-II/LC3-I) -------------------
  # y = log2(numerator/denominator) per well; actin cancels. Membrane per scheme.
  meta <- raw %>%
    mutate(rep = as.integer(sub(".*-(\\d+)$", "\\1", id)),
           LPS = factor(if_else(grepl("LPS", group), "LPS", "No LPS"),
                        levels = c("No LPS","LPS")),
           Virus = factor(case_when(grepl("PR8",group)~"PR8", grepl("S33",group)~"S33",
                                    TRUE~"None"), levels=c("None","PR8","S33")),
           group = factor(group, levels=c("C","PR8","S33","LPS","LPS/PR8","LPS/S33")))
  mk_ratio <- function(num, den, name, override = FALSE) {
    meta %>%
      transmute(id, group, LPS, Virus, rep,
                target = name, signal = NA_real_, actin = NA_real_,
                y = log2(.data[[num]] / .data[[den]]),
                membrane = if (override) factor(if_else(rep==1L,"M_r1","M_r23"))
                           else factor(paste0("M", mem))) %>%
      filter(is.finite(y))
  }
  ratios <- bind_rows(
    mk_ratio("pikb", "ikb",   "ratio.pIkB.IkB", override = FALSE),
    mk_ratio("pP65", "P65",   "ratio.pP65.P65", override = FALSE),
    mk_ratio("lc3.2","lc3.1", "ratio.LC3",      override = TRUE)
  )
  bind_rows(long, ratios) %>% mutate(membrane = factor(membrane))
}

# ==============================================================================
# §4  PHASE 2 INGESTION  (aSyn timecourse + 24h transfection panel)
# ==============================================================================

PHASE2_PSYN_VALID <- TRUE
PHASE2_RESPONSES <- c(
  y_main = "α-Synuclein (Syn1) / actin",
  y_aggr = "α-Synuclein 25 kDa (Syn1) / actin",
  y_pocc = "phospho-α-Synuclein / actin",
  y_pocc25 = "phospho-α-Synuclein 25 kDa / actin")

ingest_phase2 <- function(xlsx = DATA_XLSX) {
  raw <- read_excel(xlsx, sheet = "Phase2_Timecourse")
  long <- raw %>%
    mutate(
      across(c(syn1, syn1.25, actin, psyn, psyn_25kda), as.numeric),
      rep       = as.integer(rep),
      membrane  = factor(paste0("R", rep)),
      y_main    = log2(syn1    / actin),
      y_aggr    = log2(syn1.25 / actin),
      y_pocc    = log2(psyn    / actin),  # phospho-α-Syn level (/actin) — matches Phase-1 & cross-platform
      y_pocc25  = log2(psyn_25kda / actin),  # phospho-α-Syn 25 kDa species (/actin)
      timepoint = factor(timepoint, levels = c("4h", "8h", "24h")),
      Treatment = factor(group)
    )
  long
}

# Timecourse subset (balanced 4 treatments x 3 timepoints)
phase2_timecourse <- function(df) {
  df %>%
    filter(group %in% c("C", "PO", "PR8", "S33")) %>%
    mutate(Treatment = factor(group, levels = c("C", "PO", "PR8", "S33")),
           timepoint = factor(as.character(timepoint), levels = c("4h","8h","24h")))
}

# 24h transfection-sensing panel
phase2_panel24 <- function(df) {
  df %>%
    filter(timepoint == "24h", group %in% c("C", "PO", "T", "TPO")) %>%
    mutate(Treatment = factor(group, levels = c("C", "T", "PO", "TPO")))
}



add_branchB <- function(df, ycol, group_col, control_level,
                        match_vars = character(0)) {
  key <- c("membrane", match_vars)
  ctrl <- df %>%
    filter(.data[[group_col]] == control_level) %>%
    group_by(across(all_of(key))) %>%
    summarize(.ctrl_ref = mean(.data[[ycol]], na.rm = TRUE), .groups = "drop")
  df %>%
    left_join(ctrl, by = key) %>%
    mutate("{ycol}_fc" := .data[[ycol]] - .ctrl_ref) %>%
    select(-.ctrl_ref)
}

CMAT1 <- list(
  "LPS baseline: C vs LPS"             = c(-1, 1,  0,  0,  0,  0),
  "Infection(-LPS): C vs PR8"          = c(-1, 0,  1,  0,  0,  0),
  "Infection(-LPS): C vs S33"          = c(-1, 0,  0,  0,  1,  0),
  "Infection(-LPS): PR8 vs S33"        = c( 0, 0, -1,  0,  1,  0),
  "Co-stim(+LPS): LPS vs LPS/PR8"      = c( 0,-1,  0,  1,  0,  0),
  "Co-stim(+LPS): LPS vs LPS/S33"      = c( 0,-1,  0,  0,  0,  1),
  "Co-stim(+LPS): LPS/PR8 vs LPS/S33"  = c( 0, 0,  0, -1,  0,  1)
)

PRIORS_FIX <- c(prior(normal(0, 1.5), class = b),
                prior(student_t(3, 0, 2.5), class = sigma))
PRIORS_RE  <- c(prior(normal(0, 1.5), class = b),
                prior(student_t(3, 0, 2.5), class = sigma),
                prior(student_t(3, 0, 2.5), class = sd))

BRM <- function(formula, data, prior, ...) {
  brm(formula, data = data, family = student(), prior = prior,
      seed = BAYES_SEED, chains = MCMC_CHAINS, iter = MCMC_ITER,
      warmup = MCMC_WARMUP, cores = MCMC_CORES, refresh = 0, silent = 2, backend = "rstan",
      control = list(adapt_delta = 0.95), ...)
}

# ------------------------------------------------------------- bayes helper ---
# Describe an emmeans contrast() object -> tidy tibble (median, HDI, pd, ROPE, BF)
bayes_describe <- function(cobj, rope = ROPE_LOG2) {
  base <- as.data.frame(summary(cobj))           # gives contrast labels + medians
  dp <- tryCatch(
    bayestestR::describe_posterior(
      cobj, ci = 0.95, ci_method = "hdi",
      test = c("pd", "rope"), rope_range = c(-rope, rope), rope_ci = 1),
    error = function(e) NULL)
  # BF via Savage-Dickey refits the prior model per call (slow) -> gated off by
  # default; pd + ROPE% are the primary decision metrics. Enable with BAYES_DO_BF.
  bf <- if (isTRUE(BAYES_DO_BF)) tryCatch({
    b <- bayestestR::bayesfactor_parameters(cobj, null = 0, verbose = FALSE)
    as.data.frame(b)$log_BF
  }, error = function(e) NULL) else NULL

  out <- tibble(
    contrast = base$contrast,
    median   = if (!is.null(dp)) dp$Median   else base[[2]],
    HDI_low  = if (!is.null(dp)) dp$CI_low   else NA_real_,
    HDI_high = if (!is.null(dp)) dp$CI_high  else NA_real_,
    pd       = if (!is.null(dp)) dp$pd       else NA_real_,
    ROPE_pct = if (!is.null(dp)) dp$ROPE_Percentage else NA_real_,
    log_BF10 = if (!is.null(bf) && length(bf) == nrow(base)) bf else NA_real_
  )
  out$BF10 <- exp(out$log_BF10)
  # extra grouping cols (e.g., timepoint / Treatment) if present
  extra <- setdiff(colnames(base), c("contrast", colnames(out)))
  if (length(extra)) out <- bind_cols(base[extra], out)
  out
}

mcmc_qc <- function(fit, label) {
  conv <- posterior::summarise_draws(
            posterior::as_draws(fit),
            posterior::default_convergence_measures())
  div  <- sum(subset(brms::nuts_params(fit), Parameter == "divergent__")$Value)
  data.frame(
    model        = label,
    divergent    = div,
    max_rhat     = max(conv$rhat,     na.rm = TRUE),
    min_ess_bulk = min(conv$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(conv$ess_tail, na.rm = TRUE)
  )
}

message("[00_Config-and-Ingestion] loaded. branch A=", RUN_BRANCH_A, " B=", RUN_BRANCH_B)
