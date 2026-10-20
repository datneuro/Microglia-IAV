#!/usr/bin/env Rscript
suppressMessages({
  library(openxlsx)
  library(readxl)
})

GD  <- Sys.getenv("ANALYSIS_OUT", unset = ".")
OUT <- Sys.getenv("SUPP_TABLE_OUT", unset = file.path("Figures", "output", "supplementary_tables"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- helpers ---------------------------------------------------------------
star <- function(p) ifelse(is.na(p), "",
                    ifelse(p < .001, "***",
                    ifelse(p < .01,  "**",
                    ifelse(p < .05,  "*",
                    ifelse(p < .10,  ".", "ns")))))

fmt_p <- function(p) ifelse(is.na(p), NA_character_, formatC(p, format = "g", digits = 3))

rnd <- function(x, d = 3) if (is.numeric(x)) round(x, d) else x

vmap <- function(x) {
  x <- as.character(x)
  x[x == "V1"]   <- "A/PR/8/34"
  x[x == "V2"]   <- "A/NWS/33"
  x[x == "None"] <- "None"
  x
}

HS <- createStyle(textDecoration = "bold", fgFill = "#1F4E78",
                  fontColour = "white", halign = "center",
                  border = "TopBottomLeftRight", borderColour = "#BFBFBF",
                  wrapText = TRUE)
CAPS <- createStyle(textDecoration = "italic", fontColour = "#404040", wrapText = TRUE)

add_sheet <- function(wb, name, df, caption = NULL) {
  name <- substr(name, 1, 31)
  addWorksheet(wb, name)
  start <- 1
  if (!is.null(caption)) {
    writeData(wb, name, caption, startRow = 1, startCol = 1)
    addStyle(wb, name, CAPS, rows = 1, cols = 1)
    start <- 3
  }
  # round numeric cols (leave p-value display cols as character)
  for (cn in names(df)) if (is.numeric(df[[cn]])) df[[cn]] <- rnd(df[[cn]], 4)
  writeData(wb, name, df, startRow = start, headerStyle = HS, withFilter = TRUE)
  freezePane(wb, name, firstActiveRow = start + 1)
  setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
  invisible(wb)
}

methods_df <- function(...) {
  kv <- list(...)
  data.frame(`Statistical Dimension` = names(kv),
             `Academic Specification` = unlist(kv, use.names = FALSE),
             check.names = FALSE)
}

save_wb <- function(wb, fname) {
  fp <- file.path(OUT, fname)
  saveWorkbook(wb, fp, overwrite = TRUE)
  cat("  ✔ wrote", fname, "\n")
}

rd_csv  <- function(...) read.csv(file.path(GD, ...), check.names = FALSE, stringsAsFactors = FALSE)
rd_xlsx <- function(p, sh) as.data.frame(read_excel(file.path(GD, p), sheet = sh))

cat("SF1  Flow/ELISA factorial ...\n")
{
  an <- rd_csv("Flow-ELISA/R_analysis/output/L1_anova.csv")
  ct <- rd_csv("Flow-ELISA/R_analysis/output/L1_contrasts.csv")

  an2 <- data.frame(
    Marker        = an$Marker,
    Effect        = an$Term,
    SS            = an$SS,
    Partial_eta2  = an$partial_eta2,
    F_HC3         = an$F_hc3,
    p_HC3         = fmt_p(an$p_hc3),
    Sig           = star(an$p_hc3), check.names = FALSE)

  ct2 <- data.frame(
    Marker      = ct$Marker,
    Timepoint   = paste0(ct$Timepoint, "h"),
    Contrast    = ct$Contrast,
    Estimate    = ct$estimate,
    SE          = ct$SE,
    CI_low      = ct$CI_low,
    CI_high     = ct$CI_high,
    p_raw       = fmt_p(ct$p_raw),
    p_Holm      = fmt_p(ct$p_holm),
    Sig         = star(ct$p_holm), check.names = FALSE)

  m <- methods_df(
    "Design"               = "3 (Timepoint: 4/8/24 h) x 2 (LPS: No/Yes) x 3 (Virus: None / A/PR/8/34 / A/NWS/33), n = 5 biological replicates",
    "Strains"              = "A/PR/8/34 (PR8, H1N1); A/NWS/33 (NWS, H1N1)",
    "Markers"              = "Surface (flow, % live): CD40, CD68, CD86, CD206, HLA-I, HLA-II. Secreted (ELISA): TNF-a, IL-1b, IL-6, IL-12p70, CCL2",
    "Response transform"   = "Surface % -> logit  log(Value/(100-Value)); ELISA pg/mL -> log2",
    "Model"                = "lm(transform(Value) ~ Timepoint * LPS * Virus + Timepoint:Bio_Rep)  per marker",
    "Omnibus inference"    = "Type III SS (car::Anova, contr.sum) with HC3 robust SE (white.adjust='hc3') -> F, p; effect size = partial eta-squared",
    "Planned contrasts"    = "emmeans ~ LPS:Virus | Timepoint on HC3 vcov; pre-specified family of contrasts per timepoint",
    "Multiple testing"     = "Holm FWER per marker x contrast family",
    "Software"             = "R 4.4.2; car, emmeans, sandwich, effectsize",
    "Source"               = "01_Flow_ELISA/01_Univariate-Factorial.R -> L1_anova.csv, L1_contrasts.csv")

  wb <- createWorkbook()
  add_sheet(wb, "0_Methods", m)
  add_sheet(wb, "1_Omnibus_ANOVA", an2,
            "Type III ANOVA (HC3 robust SE). Partial eta-squared = effect size; Sig from HC3 p.")
  add_sheet(wb, "2_Planned_Contrasts", ct2,
            "Planned contrasts on the model scale (logit for %, log2 for ELISA). Estimate = log-odds / log2 difference; Sig from Holm-adjusted p.")
  save_wb(wb, "SuppTable_S3_Flow_Factorial.xlsx")
}

# =============================================================================
# SF2  — Composite / Decoupling & Polarization indices (Supplementary Table S4)
# =============================================================================
cat("SF2  Composite indices ...\n")
{
  ia <- rd_csv("Flow-ELISA/R_analysis/output/L2_index_anova.csv")
  ic <- rd_csv("Flow-ELISA/R_analysis/output/L2_index_contrasts.csv")
  lo <- rd_csv("Flow-ELISA/R_analysis/output/L2_sensitivity_LOO.csv")
  bc <- rd_csv("Flow-ELISA/R_analysis/output/L2_bootstrap_CI.csv")

  ia2 <- data.frame(Index = ia$Index, Effect = ia$Term, SS = ia$SS,
                    Partial_eta2 = ia$partial_eta2, F_HC3 = ia$F_hc3,
                    p_HC3 = fmt_p(ia$p_hc3), Sig = star(ia$p_hc3), check.names = FALSE)

  ic2 <- data.frame(Index = ic$Index, Timepoint = paste0(ic$Timepoint, "h"),
                    Contrast = ic$Contrast, Estimate = ic$estimate, SE = ic$SE,
                    CI_low = ic$CI_low, CI_high = ic$CI_high,
                    p_Holm = fmt_p(ic$p_holm), Sig = star(ic$p_holm), check.names = FALSE)

  m <- methods_df(
    "Composite construction" = "Per marker x timepoint, control-referenced effect size e = (Value_t - control_mean) / SD_all, where SD_all is the SD of all samples of the 6 groups at that timepoint (includes between-group variance). Composites = mean(e) over marker sets.",
    "Surface composite"      = "CD40, CD86, HLA-I, HLA-II",
    "Secreted composite"     = "TNF-a, IL-1b, IL-6, IL-12p70",
    "Decoupling Index"       = "Secreted composite - Surface composite (negative = surface-dominant)",
    "Polarization Index"     = "M1 composite - M2 composite (classical activation axis)",
    "Inference"              = "lm(Index ~ Timepoint*LPS*Virus + ...) ; Type III HC3 ANOVA + emmeans planned contrasts (Holm) -- IDENTICAL engine to SF1",
    "Robustness"             = "Leave-one-marker-out re-computation of the virus omnibus effect (sheet 3)",
    "Bootstrap CI"           = "Cluster bootstrap CI (sheet 4) is DESCRIPTIVE ONLY; all reported decisions are on the lm+HC3 emmeans+Holm contrasts, not on bootstrap",
    "Software"               = "R 4.4.2; car, emmeans, sandwich, boot",
    "Source"                 = "01_Flow_ELISA/02_Composite-Indices.R -> L2_index_anova.csv, L2_index_contrasts.csv, L2_sensitivity_LOO.csv, L2_bootstrap_CI.csv")

  wb <- createWorkbook()
  add_sheet(wb, "0_Methods", m)
  add_sheet(wb, "1_Index_ANOVA", ia2,
            "Type III HC3 ANOVA on each composite index (Surface, Secreted, Decoupling, Polarization).")
  add_sheet(wb, "2_Index_Contrasts", ic2,
            "Planned contrasts on indices (SD units). Estimate in control-referenced SD; Sig from Holm.")
  add_sheet(wb, "3_LeaveOneMarkerOut", lo,
            "Sensitivity: virus omnibus effect after dropping each marker. Confirms no single marker drives the composite.")
  add_sheet(wb, "4_Composite_BootCI", bc,
            "DESCRIPTIVE bootstrap CIs for composite point estimates. NOT used for any reported inferential decision.")
  save_wb(wb, "SuppTable_S4_Composite_Indices.xlsx")
}

# =============================================================================
# SF3  — Multivariate structure: PCA + PERMANOVA  (Supplementary Table S5)
# =============================================================================
cat("SF3  PCA / PERMANOVA ...\n")
{
  pm <- rd_csv("Flow-ELISA/R_analysis/output/L3_permanova.csv")
  bd <- rd_csv("Flow-ELISA/R_analysis/output/L3_betadisper.csv")
  sc <- rd_csv("Flow-ELISA/R_analysis/output/L3_scree.csv")
  ld <- rd_csv("Flow-ELISA/R_analysis/output/L3_pca_loadings.csv")

  pm2 <- data.frame(Term = pm$Term, Df = pm$Df, SumOfSqs = pm$SumOfSqs,
                    R2 = pm$R2, F = pm$F_value, p = fmt_p(pm$p_value),
                    Sig = star(pm$p_value), check.names = FALSE)
  bd2 <- data.frame(Factor = bd$Factor, F = bd$F_betadisper,
                    p = fmt_p(bd$p_betadisper), Dispersion_flag = bd$dispersion_flag,
                    check.names = FALSE)

  m <- methods_df(
    "Input"        = "Control-referenced effect-size matrix e (markers x samples), scaled & centred",
    "PCA"          = "prcomp(scale.=TRUE, center=TRUE) -> correlation PCA; eigenvalues, % variance, loadings",
    "PERMANOVA"    = "vegan::adonis2 on Euclidean distance, 9999 permutations (min attainable p = 1e-4)",
    "Dispersion"   = "vegan::betadisper + permutest (9999 perm) per factor to check homogeneity of multivariate dispersion",
    "Software"     = "R 4.4.2; vegan, stats",
    "Source"       = "01_Flow_ELISA/03_Multivariate-PCA-PERMANOVA.R -> L3_permanova.csv, L3_betadisper.csv, L3_scree.csv, L3_pca_loadings.csv")

  wb <- createWorkbook()
  add_sheet(wb, "0_Methods", m)
  add_sheet(wb, "1_PERMANOVA", pm2,
            "PERMANOVA (Euclidean, 9999 perm). R2 = variance explained; p = permutational.")
  add_sheet(wb, "2_Betadisper", bd2,
            "Multivariate dispersion homogeneity check (permutest). 'OK' = no significant dispersion confound.")
  add_sheet(wb, "3_PCA_Scree", sc, "Eigenvalues and cumulative variance per principal component.")
  add_sheet(wb, "4_PCA_Loadings", ld, "Marker loadings on PC1-PC4.")
  save_wb(wb, "SuppTable_S5_Multivariate_PCA_PERMANOVA.xlsx")
}

# =============================================================================
# SF4  — Single-cell infected vs bystander (beta-binomial GLMM) (Supplementary Table S6)
# =============================================================================
cat("SF4  Single-cell infected-vs-bystander ...\n")
{
  m2 <- rd_csv("Flow-ELISA/R-analysis-figure_VS-bystander/output/M2_contrasts.csv")
  m3a<- rd_csv("Flow-ELISA/R-analysis-figure_VS-bystander/output/M3A_contrasts.csv")
  m3b<- rd_csv("Flow-ELISA/R-analysis-figure_VS-bystander/output/M3B_interaction.csv")
  # Two analyses of the module are not tabulated here (not reported in the paper):
  #   - M4 whole-well dose-response (r = 0.97): r confounded with timepoint (15 wells = 3 clusters)
  #   - infection kinetics (%NP+): infection rates are not reported

  # M2: infected (NP+) vs bystander (NP-)  -- log-odds -> OR
  m2x <- data.frame(
    Marker     = m2$Marker, Strain = vmap(m2$Virus), LPS = m2$LPS,
    Timepoint  = paste0(m2$Tp, "h"), Contrast = m2$contrast,
    logOR      = m2$estimate, SE = m2$SE,
    OR         = exp(m2$estimate), OR_CI_low = exp(m2$CI_low), OR_CI_high = exp(m2$CI_high),
    z = m2$z, p_Holm = fmt_p(m2$p_holm), Sig = star(m2$p_holm), check.names = FALSE)

  # M3A: activation vs uninfected control (V - None) -> OR
  m3ax <- data.frame(
    Marker = m3a$Marker, LPS = m3a$LPS, Timepoint = paste0(m3a$Tp, "h"),
    Contrast = m3a$contrast, logOR = m3a$estimate, SE = m3a$SE,
    OR = exp(m3a$estimate), OR_CI_low = exp(m3a$CI_low), OR_CI_high = exp(m3a$CI_high),
    z = m3a$z, p_Holm = fmt_p(m3a$p_holm), Sig = star(m3a$p_holm), check.names = FALSE)

  # M3B: strain x state interaction (difference of log-odds) -> ratio of ORs
  m3bx <- data.frame(
    Marker = m3b$Marker, LPS = m3b$LPS, Timepoint = paste0(m3b$Tp, "h"),
    Contrast = m3b$contrast, logRoR = m3b$estimate, SE = m3b$SE,
    Ratio_of_ORs = exp(m3b$estimate), RoR_CI_low = exp(m3b$CI_low), RoR_CI_high = exp(m3b$CI_high),
    z = m3b$z, p_Holm = fmt_p(m3b$p_holm), Sig = star(m3b$p_holm), check.names = FALSE)

  # friendly contrast labels (V1/V2/None -> strain names)
  relab <- function(x) { x <- gsub("V1","A/PR/8/34",x); x <- gsub("V2","A/NWS/33",x)
                         x <- gsub("None","naive",x); x }
  m3ax$Contrast <- gsub(" - naive$"," bystander - naive", relab(m3ax$Contrast))
  m3bx$Contrast <- "A/NWS/33 vs A/PR/8/34 (autonomy difference)"

  m <- methods_df(
    "Design"            = "Per infected well, each cell gated NP+ (infected) or NP- (bystander). Timepoints 4/8/24 h; LPS No/Yes; strains A/PR/8/34, A/NWS/33; n = 5 biological replicates",
    "Model"             = "Beta-binomial GLMM (glmmTMB): cbind(positive, negative) ~ State * Virus * Timepoint + (1 | block), block = interaction(Timepoint, Bio_Rep)",
    "Effect metric"     = "Contrasts on logit scale; reported here also as Odds Ratio OR = exp(estimate). M3B interaction reported as Ratio-of-ORs.",
    "Sheet 1 (M2)"      = "CELL-AUTONOMOUS: infected (NP+) vs bystander (NP-) WITHIN each strain",
    "Sheet 2 (M3A)"     = "PARACRINE (direct): bystander (NP-) cells in infected wells vs naive control (total live). This is the genuine paracrine contrast.",
    "Sheet 3 (M3B)"     = "Strain x infection-state interaction (Ratio-of-ORs) = strain-graded autonomy",
    "Not tabulated (M4)" = "Whole-well dose-response (activation vs %NP+, r~0.97): infection rate confounded with timepoint (15 wells = 3 clusters); not reported.",
    "Not tabulated (kinetics)" = "Infection-rate (%NP+) contrasts: infection rates are not reported.",
    "Bayesian 4 h check"= "A targeted Bayesian re-estimate of the small 4 h sensing slice is in Supp. (bayes_4h_sensing) -- not duplicated here",
    "Multiple testing"  = "Holm per marker",
    "Software"          = "R 4.4.2; glmmTMB, emmeans",
    "Source"            = "01_Flow_ELISA/Infected_vs_Bystander/ -> M2_contrasts.csv, M3A_contrasts.csv, M3B_interaction.csv")

  wb <- createWorkbook()
  add_sheet(wb, "0_Methods", m)
  add_sheet(wb, "1_InfectedVsBystander", m2x,
            "(NP+) - (NP-) within each strain/timepoint. OR<1 = infected cells LOWER odds than bystanders (down-regulation/evasion); OR>1 = infected higher (cell-autonomous sensing).")
  add_sheet(wb, "2_Paracrine_BystanderVsNaive", m3ax,
            "PARACRINE: bystander (NP-) cells in infected wells vs naive control. OR>1 = paracrine up-regulation of bystanders.")
  add_sheet(wb, "3_Strain_Autonomy_Interaction", m3bx,
            "Strain x infection-state interaction. Ratio-of-ORs != 1 = strain-graded autonomy.")
  save_wb(wb, "SuppTable_S6_InfectedVsBystander.xlsx")
}

# =============================================================================
# SF5  — qPCR 3-way model  (Supplementary Table S7) -- restyles the module report
# =============================================================================
cat("SF5  qPCR ...\n")
{
  src <- "2606_PCR_VirusHMC3/analysis_master/04_master_qpcr_stats_report.xlsx"
  om  <- rd_xlsx(src, "1_Omnibus_ANOVA")
  ct  <- rd_xlsx(src, "2_All_Contrasts")
  if ("Pr(>F)" %in% names(om)) om$`Pr(>F)` <- fmt_p(suppressWarnings(as.numeric(om$`Pr(>F)`)))
  for (cc in c("P_raw","P_adjusted_Holm")) if (cc %in% names(ct)) ct[[cc]] <- fmt_p(suppressWarnings(as.numeric(ct[[cc]])))

  m <- methods_df(
    "Design"           = "2 (Timepoint: 4/8 h) x 2 (LPS: No/Yes) x 3 (Virus: None / A/PR/8/34 / A/NWS/33), randomized complete block",
    "Normalization"    = "geNORM geometric mean of GAPDH + RPS11 (NOT 2^-ddCt)",
    "Target genes"     = "IL1B, IL6, IRF7, RSAD2, MX1, OAS1, PSMB5, PSMB8",
    "Model"            = "lm(log2(FC) ~ Timepoint * LPS * Virus + Bio_Rep)  [fixed-block RCBD]",
    "Omnibus"          = "Type III SS (car::Anova, contr.sum); partial eta-squared",
    "Planned contrasts"= "7 planned contrasts per timepoint via emmeans (log2 scale)",
    "Multiple testing" = "Global Holm FWER per gene (14 factorial + 8 poly(I:C) = 22 p-values)",
    "Software"         = "R 4.4.2; car, emmeans",
    "Source"           = "02_RT-qPCR/03_Factorial-Statistics.R -> 04_master_qpcr_stats_report.xlsx (sheets 1, 2)")

  wb <- createWorkbook()
  add_sheet(wb, "0_Methods", m)
  add_sheet(wb, "1_Omnibus_ANOVA", om, "Type III ANOVA per gene on log2(FC). eta2p = partial eta-squared.")
  add_sheet(wb, "2_Planned_Contrasts", ct, "Planned contrasts on log2(FC). Holm-adjusted within gene.")
  save_wb(wb, "SuppTable_S7_qPCR.xlsx")
}

# =============================================================================
# SF6  — Proteasome activity -- restyles the module report
# =============================================================================
cat("SF6  Proteasome activity ...\n")
{
  src <- file.path("Proteasome", "processed", "Proteasome_Statistics_Log.xlsx")
  gs  <- rd_xlsx(src, "1_Group_Summary")
  ol  <- rd_xlsx(src, "2_Omnibus_LM")
  pc  <- rd_xlsx(src, "3_Planned_Contrasts")
  fa  <- rd_xlsx(src, "4_Factorial_ANOVA")
  mm  <- rd_xlsx(src, "5_Marginal_Means")
  for (cc in c("Omnibus_p","Robust_F_p")) if (cc %in% names(ol)) ol[[cc]] <- fmt_p(suppressWarnings(as.numeric(ol[[cc]])))
  for (cc in c("P_raw","P_adjusted"))     if (cc %in% names(pc)) pc[[cc]] <- fmt_p(suppressWarnings(as.numeric(pc[[cc]])))
  if ("P_Value" %in% names(fa)) fa$P_Value <- fmt_p(suppressWarnings(as.numeric(fa$P_Value)))

  m <- methods_df(
    "Design"            = "6 groups: Naive, A/PR/8/34, A/NWS/33, LPS, LPS+A/PR/8/34, LPS+A/NWS/33 (n=3)",
    "Activity readout"  = "Specific activity = total well rate - mean MG132-inhibited rate; rate = slope of lm(RFU ~ Time) over auto-detected linear window (max-R2)",
    "Immunoproteasome"  = "Ratio = Ac-ANW-AMC / Suc-LLVY-AMC, log-transformed before testing",
    "Inference"         = "lm(y ~ Group) + emmeans; HC3 robust SE",
    "Planned contrasts" = "11 pre-specified pairs via emmeans",
    "Factorial"         = "2x3 lm(y ~ LPS*Virus), Type III SS (contr.sum), HC3 robust SE; partial eta-squared",
    "Multiple testing"  = "BH within metric across 11 planned pairs",
    "Effect sizes"      = "Cohen's d (pairwise); eta-squared (omnibus); partial eta-squared (factorial)",
    "Software"          = "R 4.4.2; emmeans, car, effectsize, sandwich, lmtest",
    "Source"            = "03_Proteasome/01_Proteasome-Activity-Factorial.R -> Proteasome_Statistics_Log.xlsx")

  wb <- createWorkbook()
  add_sheet(wb, "0_Methods", m)
  add_sheet(wb, "1_Group_Summary", gs, "Per-group descriptive summary.")
  add_sheet(wb, "2_Omnibus_LM", ol, "Omnibus lm per metric with HC3 robust F and variance diagnostics.")
  add_sheet(wb, "3_Planned_Contrasts", pc, "11 planned pairwise contrasts (BH-adjusted). Cohen's d = effect size.")
  add_sheet(wb, "4_Factorial_ANOVA", fa, "2x3 factorial decomposition (LPS x Virus), Type III HC3, partial eta-squared.")
  add_sheet(wb, "5_Marginal_Means", mm, "Estimated marginal means.")
  save_wb(wb, "SuppTable_Proteasome_Activity.xlsx")
}

# =============================================================================
# SF7  — Confocal single-cell alpha-synuclein mixed models
# =============================================================================
cat("SF7  Confocal LMER ...\n")
{
  src  <- file.path("Confocal", "results_confocal_stats", "data", "Confocal_Statistics_LMER.xlsx")
  srcG <- file.path("Confocal", "results_confocal_stats", "data", "M3C_Freq_GLMM_supplement.xlsx")

  pmap <- function(df) { for (cc in intersect(c("P_raw","P_adjusted","Pr(>F)","Pr(>|t|)","P_boot"), names(df)))
                            df[[cc]] <- fmt_p(suppressWarnings(as.numeric(df[[cc]]))); df }

  pick <- list(
    c("1_Group_Summary",   "1_Group_Summary"),
    c("1b_Infection_Summary","1b_Infection_Summary"),
    c("M1_ANOVA",          "M1_Bulk_ANOVA"),
    c("M1_Pairwise",       "M1_Bulk_Pairwise"),
    c("M2_V1_ANOVA",       "M2_PR8_ANOVA"),
    c("M2_V1_Pairwise",    "M2_PR8_InfVsBystander"),
    c("M2_V2_ANOVA",       "M2_NWS_ANOVA"),
    c("M2_V2_Pairwise",    "M2_NWS_InfVsBystander"),
    c("M3A_Pairwise",      "M3A_Paracrine_Pairwise"),
    c("M3B_Pairwise",      "M3B_StrainInteraction"),
    c("M3C_Freq_Tests",    "M3C_Coloc_FreqTests"),
    c("M4A_V1_Slopes",     "M4_PR8_DoseSlopes"),
    c("M4A_V2_Slopes",     "M4_NWS_DoseSlopes"),
    c("M4B_SlopeDiff",     "M4_SlopeDifference"))

  m <- methods_df(
    "Design & context"   = "Conditions Control / A/PR/8/34 / A/NWS/33; antibodies Syn1, EPR20535. Cells classified IAV+ (infected) / IAV- (bystander) within each field",
    "Replicates"         = "n = 1 biological replicate (single-replicate study); 11 FOV, 470 cells. Module 1 (bulk) is EXPLORATORY",
    "Module 1 (bulk)"    = "lmer(metric ~ condition * Antibody + (1 | image_id))",
    "Module 2 (main)"    = "lmer(metric ~ IAV_state * Antibody + (1 | image_id)), per strain (PR8, NWS) -- infected vs bystander",
    "Module 3A"          = "lmer(metric ~ para_group * Antibody + (1 | image_id)) -- bystander vs naive (paracrine)",
    "Module 3B"          = "lmer(metric ~ condition * IAV_state + Antibody + (1 | image_id)) -- strain x state interaction",
    "Module 3C"          = "Coloc frequency (Costes p < 0.05): Fisher/chi-square tests + binomial GLMM glmer(true_coloc ~ condition * antibody + (1 | image_id)), Type III Wald tests with sum-to-zero contrasts, cluster bootstrap",
    "Module 4"           = "Dose-response: lmer(metric ~ log10(iav_mean) + Antibody + (1 | image_id)) per strain; slope difference between strains",
    "Correlation metric" = "Pearson r with Fisher-Z transform for inference",
    "Random effect"      = "FOV (image_id); lm fallback when single FOV. Cluster bootstrap CIs reported where available (Boot_* columns)",
    "Multiple testing"   = "Holm within metric",
    "Software"           = "R 4.4.2; lme4, lmerTest, emmeans, car",
    "Source"             = "05_Confocal/04_Single-Cell-Mixed-Models.R -> Confocal_Statistics_LMER.xlsx; 05_Confocal/06_Coloc-Frequency-GLMM.R -> M3C_Freq_GLMM_supplement.xlsx")

  wb <- createWorkbook()
  add_sheet(wb, "0_Methods", m)
  for (pp in pick) {
    df <- tryCatch(pmap(rd_xlsx(src, pp[1])), error = function(e) NULL)
    if (!is.null(df)) add_sheet(wb, pp[2], df)
  }
  # M3C GLMM supplement (separate file)
  for (sh in c("M3C_Freq_GLMM_ANOVA","M3C_Freq_GLMM_Pairwise")) {
    df <- tryCatch(pmap(rd_xlsx(srcG, sh)), error = function(e) NULL)
    if (!is.null(df)) add_sheet(wb, sh, df)
  }
  save_wb(wb, "SuppTable_Confocal_aSyn_MixedModels.xlsx")
}

cat("\nDONE. Output folder:\n  ", OUT, "\n")
