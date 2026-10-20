# ==============================================================================
# 06_Supplementary-Tables.R  —  Statistical supplement
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(openxlsx) })
PROJ <- "."
DIR  <- file.path(PROJ,"CrossPlatform_aSyn"); OUT <- file.path(DIR,"output")

options("openxlsx.borderStyle" = "thin")
hdr_style  <- createStyle(textDecoration="bold", fgFill="#2F3E46", fontColour="white", halign="center", border="Bottom", fontSize=11)
data_style <- createStyle(halign="center", fontSize=10)
text_style <- createStyle(halign="left",   fontSize=10)
pass_style <- createStyle(fgFill="#D5E8D4")     # academic pastel green = passed/stable

add_sheet <- function(wb, name, df, green_rows=NULL){
  addWorksheet(wb, name)
  df <- df %>% mutate(across(where(is.numeric), ~ round(., 3)))
  writeData(wb, name, df)
  for (col in seq_len(ncol(df))) {
    ml <- max(c(nchar(colnames(df)[col]), nchar(as.character(df[[col]]))), na.rm=TRUE)
    setColWidths(wb, name, cols=col, widths=min(max(ml+4,12),48))
  }
  addStyle(wb, name, hdr_style,  rows=1, cols=seq_len(ncol(df)), gridExpand=TRUE)
  addStyle(wb, name, data_style, rows=seq_len(nrow(df))+1L, cols=seq_len(ncol(df)), gridExpand=TRUE)
  addStyle(wb, name, text_style, rows=seq_len(nrow(df))+1L, cols=1:2, gridExpand=TRUE)
  if (length(green_rows)) addStyle(wb, name, pass_style, rows=green_rows+1L, cols=seq_len(ncol(df)), gridExpand=TRUE, stack=TRUE)
}

# ---- read raw outputs --------------------------------------------------------
tierA <- read_csv(file.path(OUT,"tierA_correlations.csv"), show_col_types=FALSE)
wbint <- read_csv(file.path(OUT,"wb_internal_corr.csv"),   show_col_types=FALSE)
qcA   <- read_csv(file.path(OUT,"tierA_qc.csv"),           show_col_types=FALSE)
qcW   <- read_csv(file.path(OUT,"wb_internal_qc.csv"),     show_col_types=FALSE)
sens  <- read_csv(file.path(OUT,"crossplatform_sensitivity.csv"), show_col_types=FALSE)

wb <- createWorkbook()

# ---- Sheet 0: Methods --------------------------------------------------------
methods <- tibble(
  `Statistical Dimension` = c("Experimental Pairing","Estimand","Likelihood & Priors",
    "MCMC Sampling","Convergence Gate","ROPE Decision Rationale","Robustness Variants","Sensitivity Analysis"),
  `Academic Specification` = c(
    "Each culture dish split across Western blot, RT-qPCR, ELISA and flow cytometry, yielding 48 timepoint-matched samples (between-subjects per timepoint). WB same-lysate panel uses strictly intra-lysate matching from identical 24 h extracts (n = 18).",
    "Within-condition partial correlation (rho) between the 25 kDa alpha-synuclein species and each marker, obtained as the residual correlation (rescor) of a single-stage Bayesian multivariate model adjusting for Timepoint x Treatment (cross-platform) or group (same-lysate). Associations are correlational, not causal.",
    "Gaussian multivariate model on z-scaled responses, set_rescor(TRUE). Weakly-informative priors: regression coefficients b ~ Normal(0, 1.5); residual correlation rescor ~ LKJ(2).",
    "rstan/brms, 4 chains, adapt_delta = 0.99. Cross-platform: 20,000 post-warmup draws. WB same-lysate (n=18): 28,000 post-warmup draws.",
    "Pre-specified gate on the estimand (rescor): 0 divergent transitions, Rhat <= 1.01, and effective sample size (bulk AND tail) > 10,000.",
    "Region of Practical Equivalence on the correlation scale, |rho| < 0.10 (negligible association). ROPE Mass (%) = posterior mass inside that band.",
    "Total = unadjusted Bayesian Pearson correlation; NP-adjusted = additionally adjusting for the infected-cell fraction (%NP+); 24 h = restricted to the 24 h block.",
    "Each partial correlation refit under default / broader / narrower priors and re-classified at ROPE half-widths 0.05 / 0.10 / 0.15. Prior-stable = sign and ROPE-class invariant across priors; ROPE-window-stable = ROPE-class invariant across the three widths."
  )
)
addWorksheet(wb, "0_Methods_Summary"); writeData(wb, "0_Methods_Summary", methods)
addStyle(wb, "0_Methods_Summary", hdr_style, rows=1, cols=1:2, gridExpand=TRUE)
setColWidths(wb, "0_Methods_Summary", cols=1:2, widths=c(28,115))
addStyle(wb, "0_Methods_Summary", text_style, rows=2:(nrow(methods)+1), cols=1:2, gridExpand=TRUE)

# ---- Sheet 1: Tier A cross-platform partial correlations ---------------------
rob <- function(df, mk, mod, col){ x<-df[df$marker==mk & df$model==mod,]; if(nrow(x)) x[[col]][1] else NA_real_ }
A_within <- tierA %>% filter(model=="within")
A_tab <- A_within %>% transmute(
  `Marker (display)`=label, Platform=platform, n=n,
  `Partial r (within-condition)`=rho, `95% HDI Lower`=CI_low, `95% HDI Upper`=CI_high,
  `Directional Prob (pd)`=pd, `ROPE Mass (%)`=100*ROPE_pct,
  `Total r (unadjusted)`=map_dbl(marker, rob, df=tierA, mod="total", col="rho"),
  `NP-adjusted r`=map_dbl(marker, rob, df=tierA, mod="np_adj", col="rho")) %>%
  arrange(desc(abs(`Partial r (within-condition)`)))
add_sheet(wb, "1_TierA_CrossPlatform", A_tab)

# ---- Sheet 2: WB same-lysate partial correlations ----------------------------
W_within <- wbint %>% filter(model=="within")
W_tab <- W_within %>% transmute(
  `Marker / Pathway (display)`=label, n=n,
  `Partial r (same lysate)`=rho, `95% HDI Lower`=CI_low, `95% HDI Upper`=CI_high,
  `Directional Prob (pd)`=pd, `ROPE Mass (%)`=100*ROPE_pct,
  `Total r (unadjusted)`=map_dbl(target, function(t){x<-wbint[wbint$target==t & wbint$model=="total",]; if(nrow(x)) x$rho[1] else NA_real_})) %>%
  arrange(desc(abs(`Partial r (same lysate)`)))
add_sheet(wb, "2_WB_SameLysate", W_tab)

# ---- Sheet 3: MCMC convergence QC --------------------------------------------
qc_tab <- bind_rows(
  qcA %>% mutate(Tier="Cross-platform (Tier A)", Marker=sub("tierA:","",model)),
  qcW %>% mutate(Tier="WB same-lysate",          Marker=sub("wb_internal:","",model))) %>%
  transmute(Tier, Marker, `Divergent Transitions`=divergent, `Max R-hat`=max_rhat,
    `rescor ESS (Bulk)`=rescor_ess_bulk, `rescor ESS (Tail)`=rescor_ess_tail,
    `Min ESS (all params)`=min_ess_bulk, `Gate Passed`=if_else(pass,"PASS","FAIL"))
green_qc <- which(qc_tab$`Gate Passed`=="PASS")
add_sheet(wb, "3_MCMC_Convergence_QC", qc_tab, green_rows=green_qc)

# ---- Sheet 4: Prior + ROPE-width sensitivity ---------------------------------
S_tab <- sens %>% transmute(
  Tier=if_else(tier=="tierA","Cross-platform","WB same-lysate"),
  `Marker (display)`=label,
  `Prior Configuration`=recode(prior_set, default="Standard Weakly-Informative",
                               broader="Broader Variance", narrower="Narrower Variance"),
  `Posterior Median r`=median, `95% HDI Lower`=HDI_low, `95% HDI Upper`=HDI_high,
  `Directional Prob (pd)`=pd, `ROPE Mass (%) [|0.1|]`=100*ROPE_pct_10,
  `Decision (ROPE |0.1|)`=class_10,
  `Prior Stable?`=if_else(prior_stable,"Stable","Sensitive"),
  `ROPE Window Stable?`=if_else(rope_stable,"Stable","Sensitive")) %>%
  arrange(Tier, desc(abs(`Posterior Median r`)), `Prior Configuration`)
green_s <- which(S_tab$`Prior Stable?`=="Stable" & S_tab$`ROPE Window Stable?`=="Stable")
add_sheet(wb, "4_Prior_ROPE_Sensitivity", S_tab, green_rows=green_s)

saveWorkbook(wb, file.path(OUT,"CrossPlatform_Supplementary_Data.xlsx"), overwrite=TRUE)
cat("\n[06_Supplementary-Tables] Saved CrossPlatform_Supplementary_Data.xlsx (5 sheets)\n")
cat(sprintf("  Tier A markers: %d | WB same-lysate markers: %d | QC pass: %d/%d | sensitivity rows: %d\n",
            nrow(A_tab), nrow(W_tab), sum(qc_tab$`Gate Passed`=="PASS"), nrow(qc_tab), nrow(S_tab)))
