# ==============================================================================
# _Labels.R  
# ==============================================================================

MK_LABEL <- c(
  # ELISA (secreted protein)
  ELISA_CCL2="CCL2", ELISA_IL_1b="IL-1β", ELISA_IL_6="IL-6", ELISA_IL_10="IL-10",
  ELISA_TNFa="TNF-α", ELISA_IL_12p70="IL-12p70",
  # Flow cytometry (surface protein)
  FACS_CD40="CD40", FACS_CD68="CD68", FACS_CD86="CD86", FACS_CD206="CD206",
  FACS_HLA_I="HLA-I", FACS_HLA_II="HLA-II", FACS_NP="NP (influenza)",
  # RT-qPCR (transcript)
  qPCR_IL1B="IL1B mRNA", qPCR_IL6="IL6 mRNA", qPCR_IRF7="IRF7 mRNA",
  qPCR_RSAD2="RSAD2 mRNA", qPCR_MX1="MX1 mRNA", qPCR_OAS1="OAS1 mRNA",
  qPCR_PSMB5="PSMB5 mRNA", qPCR_PSMB8="PSMB8 mRNA",
  # WB (protein) — names matched to Phase-1 WB targets
  WB_aSyn_25kDa="25-kDa α-Synuclein (Syn1)",
  WB_aSyn_multimer="α-Synuclein multimer (Syn1)",
  WB_aSyn_pSyn="phospho-α-Synuclein"
)
disp_marker <- function(x) ifelse(x %in% names(MK_LABEL), unname(MK_LABEL[x]), x)

PLAT_LABEL <- c(WB="Western blot", ELISA="ELISA", FACS="Flow cytometry", qPCR="RT-qPCR")
disp_platform <- function(x) ifelse(x %in% names(PLAT_LABEL), unname(PLAT_LABEL[x]), x)
PLAT_LEVELS <- c("ELISA","Flow cytometry","RT-qPCR","Western blot")
PLAT_COL    <- c("ELISA"="#E69F73","Flow cytometry"="#7BA6C9","RT-qPCR"="#7FB28A","Western blot"="#9B8BBF")

AXIS_COL <- c("Proteostasis"="#5AA786","NF-κB"="#C56B6B","Autophagy–lysosome"="#6A93C0",
              "α-Synuclein species"="#9B7BB5","Inflammasome"="#D9A066","Other"="grey55")

TRT_LABEL <- c(Control="Control", V1="A/PR/8/34", V2="A/NWS/33", PolyIC="Poly(I:C)")
disp_treat <- function(x) ifelse(x %in% names(TRT_LABEL), unname(TRT_LABEL[x]), x)
TRT_LEVELS <- c("Control","A/PR/8/34","A/NWS/33","Poly(I:C)")
TRT_COL    <- c("Control"="#7F8C8D","A/PR/8/34"="#4285F4","A/NWS/33"="#C0392B","Poly(I:C)"="#27AE60")

pd_stars    <- function(pd) ifelse(pd>=0.995,"***",ifelse(pd>=0.985,"**",ifelse(pd>=0.975,"*","")))
pd_stars_md <- function(pd) gsub("\\*","&ast;", pd_stars(pd))

fmt_bf <- function(bf) ifelse(is.na(bf),"—", ifelse(bf>=10, sprintf("%.0f",bf), sprintf("%.1f",bf)))
stat_md <- function(rho, lo, hi, pd, bf=NA) sprintf(
  "%s *ρ* = %+.2f [%.2f, %.2f], *p*<sub>d</sub> = %.2f, BF<sub>10</sub> = %s",
  pd_stars_md(pd), rho, lo, hi, pd, fmt_bf(bf))

stat_inline_md <- function(rho, lo, hi, pd) {
  s <- sprintf("%s *ρ* = %+.2f [%.2f, %.2f], *p*<sub>d</sub> = %.2f",
               pd_stars_md(pd), rho, lo, hi, pd)
  ifelse(pd >= 0.985, paste0("**", s, "**"), s)
}
