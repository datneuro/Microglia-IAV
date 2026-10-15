# ==============================================================================
# Run-All.R  —  Reproduce the full Infected-vs-Bystander analysis end to end.
#   Run from the repository root:  Rscript 01_Flow_ELISA/Infected_vs_Bystander/Run-All.R
#   Data/output locations are set by BASE / HERE at the top of each script.
# Stages run in order; each writes to output/ (CSV) and figure/ (PDF+PNG).
# ==============================================================================
CODE <- "01_Flow_ELISA/Infected_vs_Bystander"
steps <- c(file.path(CODE, "00_Count-Preprocessing.R"),
           file.path(CODE, "01_Infection-Kinetics.R"),
           file.path(CODE, "02_Infected-vs-Bystander-GLMM.R"),
           file.path(CODE, "03_Paracrine-Activation.R"),
           file.path(CODE, "04_Dose-Response.R"),
           file.path(CODE, "05_Bayesian-4h-Sensing.R"))   
for (s in steps) {
  cat("\n################  ", basename(s), "  ################\n")
  source(s)
}
cat("\n ALL STAGES COMPLETE — see output/.\n",
    "   Panels f-g of Figure 1 are drawn by Figures/Figure-01_Microglial-Activation.R\n")
