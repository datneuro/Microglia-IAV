# ==============================================================================
# 04_Supplementary-Tables.R 
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

if (!exists("OUT_DIR")) source(file.path("04_Western_Blot", "00_Config-and-Ingestion.R"))

message("\n########## COMPILING PHASE-SPLIT SUPPLEMENTARY workBOOKS ##########")

# ------------------------------------------------------------------------------
# §1  HUMAN-READABLE LABEL DICTIONARIES
# ------------------------------------------------------------------------------
humanize_target <- function(x) {
  case_when(
    x == "y_main"        ~ "α-Synuclein monomer (Syn1)",
    x == "y_aggr"        ~ "α-Synuclein 25 kDa species (Syn1)",
    x == "y_pocc"        ~ "phospho-α-Synuclein (p-α-Syn)",
    x == "y_pocc25"      ~ "phospho-α-Synuclein 25 kDa species",
    x == "syn1.mono"     ~ "α-Synuclein monomer (Syn1)",
    x == "syn1.multi"    ~ "α-Synuclein multimer (Syn1)",
    x == "syn1.25"       ~ "α-Synuclein 25 kDa species (Syn1)",
    x == "syn211"        ~ "α-Synuclein (Syn211)",
    x == "psyn"          ~ "phospho-α-Synuclein (p-α-Syn)",
    x == "psmb5"         ~ "PSMB5 (Mature)",
    x == "psmb5.pre"     ~ "PSMB5 precursor",
    x == "k48"           ~ "K48-polyubiquitin chains",
    x == "psmb8"         ~ "PSMB8",
    x == "asc"           ~ "ASC (Inflammasome adaptor)",
    x == "pikb"          ~ "phospho-IκBα",
    x == "ikb"           ~ "IκBα (Total)",
    x == "pP65"          ~ "phospho-NF-κB p65",
    x == "P65"           ~ "NF-κB p65 (Total)",
    x == "lc3.1"         ~ "LC3-I",
    x == "lc3.2"         ~ "LC3-II",
    x == "pre.CTSD"      ~ "Cathepsin D precursor (pro-CTSD)",
    x == "CTSD"          ~ "Cathepsin D (Mature lysosomal)",
    x == "ratio.pIkB.IkB" ~ "p-IκBα / IκBα Activation Ratio",
    x == "ratio.pP65.P65" ~ "p-p65 / p65 Activation Ratio",
    x == "ratio.LC3"      ~ "LC3-II / LC3-I Autophagy Ratio",
    TRUE ~ x
  )
}

humanize_contrast <- function(x) {
  case_when(
    x == "LPS baseline: C vs LPS"            ~ "LPS Baseline Effect (Control vs LPS)",
    x == "Infection(-LPS): C vs PR8"         ~ "PR8 Infection Effect without LPS (Control vs PR8)",
    x == "Infection(-LPS): C vs S33"         ~ "S33 Infection Effect without LPS (Control vs S33)",
    x == "Infection(-LPS): PR8 vs S33"       ~ "Differential Infection Profile (PR8 vs S33)",
    x == "Co-stim(+LPS): LPS vs LPS/PR8"     ~ "PR8 Co-stimulation Effect (LPS vs LPS/PR8)",
    x == "Co-stim(+LPS): LPS vs LPS/S33"     ~ "S33 Co-stimulation Effect (LPS vs LPS/S33)",
    x == "Co-stim(+LPS): LPS/PR8 vs LPS/S33" ~ "Differential Co-stimulation Profile (LPS/PR8 vs LPS/S33)",
    x == "LPS main"                          ~ "LPS Main Factorial Effect",
    x == "Virus main: PR8"                   ~ "PR8 Infection Main Factorial Effect",
    x == "Virus main: S33"                   ~ "S33 Infection Main Factorial Effect",
    x == "Virus×LPS: PR8"                    ~ "PR8 × LPS Factorial Synergy/Interaction",
    x == "Virus×LPS: S33"                    ~ "S33 × LPS Factorial Synergy/Interaction",
    x == "PO - C"                            ~ "Poly:IC Treatment vs Control",
    x == "PR8 - C"                           ~ "PR8 Infection vs Control",
    x == "S33 - C"                           ~ "S33 Infection vs Control",
    x == "T - C"                             ~ "Transfection Vehicle vs Control",
    x == "TPO - T"                           ~ "Intracellular dsRNA Sensing (TPO vs Transfection Vehicle)",
    x == "TPO - PO"                          ~ "Intra- vs Extra-cellular dsRNA Sensing (TPO vs Poly:IC)",
    x == "PO - C"                            ~ "Extracellular dsRNA Sensing (Poly:IC vs Control)",
    x == "linear"                            ~ "Linear Dynamic Time Trend",
    x == "quadratic"                         ~ "Quadratic Dynamic Time Trend",
    TRUE ~ x
  )
}

# ------------------------------------------------------------------------------
# §2  WORKBOOK STYLES
# ------------------------------------------------------------------------------
options("openxlsx.borderStyle" = "thin")
hdr_style  <- createStyle(textDecoration = "bold", fgFill = "#2F3E46", fontColour = "white", halign = "center", border = "Bottom", fontSize = 11)
data_style <- createStyle(halign = "center", fontSize = 10)
text_style <- createStyle(halign = "left", fontSize = 10)
pass_style <- createStyle(fgFill = "#D5E8D4") # pastel green: prior- and ROPE-stable rows

add_polished_sheet <- function(wb, name, df, is_sens = FALSE) {
  addWorksheet(wb, name)
  
  # Round numeric columns to 3 decimals
  if (any(sapply(df, is.numeric))) {
    df <- df %>% mutate(across(where(is.numeric), ~ round(., 3)))
  }
  writeData(wb, name, df)
  
  for (col in seq_len(ncol(df))) {
    max_len <- max(c(nchar(colnames(df)[col]), nchar(as.character(df[[col]]))), na.rm = TRUE)
    setColWidths(wb, name, cols = col, widths = min(max(max_len + 4, 12), 45))
  }
  
  addStyle(wb, name, hdr_style, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
  addStyle(wb, name, data_style, rows = seq_len(nrow(df)) + 1L, cols = seq_len(ncol(df)), gridExpand = TRUE)
  addStyle(wb, name, text_style, rows = seq_len(nrow(df)) + 1L, cols = 1:2, gridExpand = TRUE)
  
  if (is_sens) {
    ok_rows <- which(df$`Prior Stable?` == "Stable" & df$`ROPE Window Stable?` == "Stable") + 1L
    if (length(ok_rows) > 0) {
      addStyle(wb, name, pass_style, rows = ok_rows, cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
    }
  }
}

# MCMC QC summary from the workbook written by 02_Bayesian-Densitometry.R
master_report_path <- file.path(OUT_DIR, "WB_Stats_Report.xlsx")
master_qc_table <- if (file.exists(master_report_path)) {
  openxlsx::readWorkbook(master_report_path, sheet = "8_MCMC_QC_Summary")
} else {
  NULL
}

# Prior/ROPE sensitivity results from 03_Sensitivity-Analysis.R
sens_raw <- read_csv(file.path(OUT_DIR, "WB_sensitivity.csv"), show_col_types = FALSE) %>% filter(branch == "A_model_membrane")

# ------------------------------------------------------------------------------
# §3  COMPILE FILE 1: PHASE 1 SUPPLEMENTARY WORKBOOK
# ------------------------------------------------------------------------------
message("  Building: WB_Supplementary_Data_Phase1.xlsx...")
wb1 <- createWorkbook()

# Sheet 0: Methods
p1_methods <- tibble(
  `Statistical Dimension` = c("Experimental Design Phase 1", "Response Transformation", "Likelihood Specification", "Prior Structures", "ROPE Decision Rationale"),
  `Academic Specification` = c(
    "2 (LPS: No/Yes) x 3 (Virus: None/PR8/S33) Factorial Screen, n=3 independent biological replicates. Sampled at 24h post-treatment.",
    "Log2-transformed target-to-actin densitometric ratio. Pathway activation quantified via explicit internal log2 ratios (e.g., phospho-p65 / total p65) which naturally cancels loading variance.",
    "Bayesian Fixed-Block Factorial Model structured as: response ~ LPS * Virus + membrane. Heavy-tailed Student-t family applied to insulate posteriors from experimental gel outliers.",
    "Weakly-informative priors isolating parameter scale. Directional coefficients: b ~ Normal(0, 1.5). Residual scale variance: sigma ~ Student-t(3, 0, 2.5).",
    "Decision framework utilizes the Region of Practical Equivalence (ROPE) on continuous marginal posteriors, bounded at [-0.5, +0.5] log2 scale (equivalent to an intuitive ~1.4-fold minimal biological variation boundary)."
  )
)
addWorksheet(wb1, "0_Methods_Summary")
writeData(wb1, "0_Methods_Summary", p1_methods)
addStyle(wb1, "0_Methods_Summary", hdr_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb1, "0_Methods_Summary", cols = 1:2, widths = c(30, 100))
addStyle(wb1, "0_Methods_Summary", text_style, rows = 2:6, cols = 1:2, gridExpand = TRUE)

# Sheet 1: Results
p1_bayes <- read_csv(file.path(OUT_DIR, "WB_phase1_bayes_contrasts.csv"), show_col_types = FALSE) %>% filter(branch == "A_model_membrane")
p1_effs  <- read_csv(file.path(OUT_DIR, "WB_phase1_factor_effects.csv"), show_col_types = FALSE)
p1_clean <- bind_rows(p1_bayes %>% mutate(Block = "Planned Comparison"), p1_effs %>% mutate(Block = "Main Effect / Interaction", contrast = effect)) %>%
  mutate(`Target Protein / Pathway Ratio` = humanize_target(target), `Experimental Comparison` = humanize_contrast(contrast)) %>%
  select(`Target Protein / Pathway Ratio`, `Experimental Comparison`, Block, `Posterior Median` = median, `95% HDI Lower` = HDI_low, `95% HDI Upper` = HDI_high, `Directional Prob (pd)` = pd, `ROPE Mass (%)` = ROPE_pct) %>%
  arrange(`Target Protein / Pathway Ratio`, desc(Block))
add_polished_sheet(wb1, "1_Phase1_Bayes_Factorial", p1_clean)

# Sheet 2: MCMC QC
if (!is.null(master_qc_table)) {
  p1_qc <- master_qc_table %>% 
    filter(grepl("A_model_membrane", model) & !grepl("y_main|y_aggr|y_pocc", model)) %>%
    mutate(`Target Model Block` = humanize_target(sub(" *\\(.*", "", model))) %>%
    select(`Target Model Block`, `Divergent Transitions` = divergent, `Max R-hat (Convergence)` = max_rhat, `Min Bulk ESS` = min_ess_bulk, `Min Tail ESS` = min_ess_tail)
  add_polished_sheet(wb1, "2_P1_MCMC_Convergence_QC", p1_qc)
}

# Sheet 3: Sensitivity — all contrasts classified "reject" (outside ROPE), plus
# the NF-κB ratio contrasts reported in the Results (p-IκBα/IκBα: C vs PR8 and
# C vs S33; p-p65/p65: C vs S33 and LPS vs LPS/S33), which are listed even when
# classified "undecided".
p1_sens_filtered <- sens_raw %>%
  filter(
    phase == "P1" & (
      class_05 == "reject" | 
      (response == "ratio.pIkB.IkB" & contrast %in% c("Infection(-LPS): C vs PR8", "Infection(-LPS): C vs S33")) |
      (response == "ratio.pP65.P65" & contrast %in% c("Infection(-LPS): C vs S33", "Co-stim(+LPS): LPS vs LPS/S33"))
    )
  ) %>%
  mutate(
    `Target Protein Readout`  = humanize_target(response),
    `Experimental Comparison` = humanize_contrast(contrast),
    `Prior Configuration`     = case_when(
       prior_set == "default"  ~ "Standard Weakly-Informative",
       prior_set == "broader"  ~ "Broader Variance",
       prior_set == "narrower" ~ "Narrower Variance", TRUE ~ prior_set
    ),
    `Prior Stable?`           = if_else(prior_stable == TRUE, "Stable", "Sensitive"),
    `ROPE Window Stable?`     = if_else(rope_stable == TRUE, "Stable", "Sensitive")
  ) %>%
  select(`Target Protein Readout`, `Experimental Comparison`, `Prior Configuration`, `Posterior Median` = median, `95% HDI Lower` = HDI_low, `95% HDI Upper` = HDI_high, 
         `Directional Prob (pd)` = pd, `ROPE Mass (%)` = ROPE_pct, `Decision (ROPE)` = class_05, `Prior Stable?`, `ROPE Window Stable?`) %>%
  arrange(`Target Protein Readout`, `Experimental Comparison`)
add_polished_sheet(wb1, "3_P1_Prior_ROPE_Sensitivity", p1_sens_filtered, is_sens = TRUE)

saveWorkbook(wb1, file.path(OUT_DIR, "WB_Supplementary_Data_Phase1.xlsx"), overwrite = TRUE)

# ------------------------------------------------------------------------------
# §4  COMPILE FILE 2: PHASE 2 SUPPLEMENTARY WORKBOOK (TIMECOURSE)
# ------------------------------------------------------------------------------
message("  Building: WB_Supplementary_Data_Phase2.xlsx...")
wb2 <- createWorkbook()

# Sheet 0: Methods
p2_methods <- tibble(
  `Statistical Dimension` = c("Experimental Design Phase 2", "Response Transformation", "Likelihood Specification", "Prior Structures", "ROPE Decision Rationale"),
  `Academic Specification` = c(
    "Balanced Treatment grid (Control, Poly:IC, PR8, S33) across 3 dynamic timepoints (4h, 8h, 24h), n=4. Supported by a 24h Transfection sensing grid (Control, Vehicle, PO, TPO).",
    "Log2-transformed target-to-actin densitometric density ratio.",
    "Bayesian Fixed-Block Randomized Complete Block Design (RCBD) structured as: response ~ Timepoint * Treatment + membrane. Heavy-tailed Student-t family applied to penalize gel outliers.",
    "Weakly-informative priors isolating parameter scale. Directional coefficients: b ~ Normal(0, 1.5). Residual scale variance: sigma ~ Student-t(3, 0, 2.5).",
    "Decision framework utilizes the Region of Practical Equivalence (ROPE) on continuous marginal posteriors, bounded at [-0.5, +0.5] log2 scale (equivalent to an intuitive ~1.4-fold minimal biological variation boundary)."
  )
)
addWorksheet(wb2, "0_Methods_Summary")
writeData(wb2, "0_Methods_Summary", p2_methods)
addStyle(wb2, "0_Methods_Summary", hdr_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb2, "0_Methods_Summary", cols = 1:2, widths = c(30, 100))
addStyle(wb2, "0_Methods_Summary", text_style, rows = 2:6, cols = 1:2, gridExpand = TRUE)

# Sheet 1: Results
p2_temp  <- read_csv(file.path(OUT_DIR, "WB_phase2_temporal.csv"), show_col_types = FALSE) %>% filter(branch == "A_model_membrane")
p2_panel <- read_csv(file.path(OUT_DIR, "WB_phase2_panel24.csv"), show_col_types = FALSE)
p2_clean <- bind_rows(p2_temp %>% mutate(Block = "Timecourse Dynamics (4h/8h/24h)"), p2_panel %>% mutate(Block = "Transfection Panel (24h)", timepoint = "24h")) %>%
  mutate(`Target Readout` = humanize_target(response), `Comparison Metric` = humanize_contrast(contrast), `Timepoint` = timepoint) %>%
  select(`Target Readout`, `Comparison Metric`, Block, `Timepoint`, `Posterior Median` = median, `95% HDI Lower` = HDI_low, `95% HDI Upper` = HDI_high, `Directional Prob (pd)` = pd, `ROPE Mass (%)` = ROPE_pct) %>%
  arrange(`Target Readout`, desc(Block), `Timepoint`)
add_polished_sheet(wb2, "1_Phase2_aSyn_Dynamics", p2_clean)

# Sheet 2: MCMC QC
if (!is.null(master_qc_table)) {
  p2_qc <- master_qc_table %>% 
    filter(!grepl("A_model_membrane", model) | grepl("y_main|y_aggr|y_pocc", model)) %>%
    mutate(`Analyzed Model Block` = case_when(
      grepl("y_main", model) ~ "α-Synuclein monomer (Syn1 Timecourse)",
      grepl("y_aggr", model) ~ "α-Synuclein 25 kDa species (Syn1 Timecourse)",
      grepl("y_pocc25", model) ~ "phospho-α-Synuclein 25 kDa species (Timecourse)",
      grepl("y_pocc", model) ~ "phospho-α-Synuclein (p-α-Syn Timecourse)", TRUE ~ model
    )) %>%
    select(`Analyzed Model Block`, `Divergent Transitions` = divergent, `Max R-hat (Convergence)` = max_rhat, `Min Bulk ESS` = min_ess_bulk, `Min Tail ESS` = min_ess_tail)
  add_polished_sheet(wb2, "2_P2_MCMC_Convergence_QC", p2_qc)
}

# Sheet 3: Sensitivity
p2_sens_filtered <- sens_raw %>%
  filter(phase == "P2" & class_05 == "reject") %>%
  mutate(
    `Target Protein Readout`  = humanize_target(response),
    `Experimental Comparison` = humanize_contrast(contrast),
    `Timepoint`               = as.character(timepoint),
    `Prior Configuration`     = case_when(
       prior_set == "default"  ~ "Standard Weakly-Informative",
       prior_set == "broader"  ~ "Broader Variance",
       prior_set == "narrower" ~ "Narrower Variance", TRUE ~ prior_set
    ),
    `Prior Stable?`           = if_else(prior_stable == TRUE, "Stable", "Sensitive"),
    `ROPE Window Stable?`     = if_else(rope_stable == TRUE, "Stable", "Sensitive")
  ) %>%
  select(`Target Protein Readout`, `Experimental Comparison`, `Timepoint`, `Prior Configuration`, `Posterior Median` = median, `95% HDI Lower` = HDI_low, `95% HDI Upper` = HDI_high, 
         `Directional Prob (pd)` = pd, `ROPE Mass (%)` = ROPE_pct, `Decision (ROPE)` = class_05, `Prior Stable?`, `ROPE Window Stable?`) %>%
  arrange(`Target Protein Readout`, `Timepoint`, `Experimental Comparison`)
add_polished_sheet(wb2, "3_P2_Prior_ROPE_Sensitivity", p2_sens_filtered, is_sens = TRUE)

saveWorkbook(wb2, file.path(OUT_DIR, "WB_Supplementary_Data_Phase2.xlsx"), overwrite = TRUE)

message("\n====================================================================")
message("  SUCCESS: Two independent publication workbooks generated!")
message("  - File 1: WB_Supplementary_Data_Phase1.xlsx")
message("  - File 2: WB_Supplementary_Data_Phase2.xlsx")
message("====================================================================\n")