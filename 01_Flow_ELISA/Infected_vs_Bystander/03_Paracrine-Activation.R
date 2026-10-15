# ==============================================================================
# 03_Paracrine-Activation.R 
# ==============================================================================

BASE <- "./Flow-ELISA"                                     # data folder (contains DATA/)
HERE <- file.path(BASE, "R-analysis-figure_VS-bystander")  # module output folder
OUT  <- file.path(HERE, "output")
CODE <- "01_Flow_ELISA/Infected_vs_Bystander"              # module scripts (repo-relative)
source(file.path(CODE, "_Model-Helpers.R"))
source(file.path(CODE, "_Style-Bystander.R"))
suppressPackageStartupMessages({ library(purrr); library(tidyr) })

ov <- read.csv(file.path(OUT, "counts_overall_long.csv")) %>%
  mutate(Virus = factor(Virus, levels = c("None","V1","V2")),
         LPS   = factor(LPS, levels = c("No","Yes")),
         Tp    = factor(Timepoint, levels = c(4,8,24)),
         block = interaction(Timepoint, Bio_Rep, drop = TRUE))

emm_l <- ct_l <- list()
for (lps in c("No","Yes")) {
  cat("\n=========== LPS =", lps, "===========\n")
  for (mk in GATE_MARKERS) {
    d <- ov %>% filter(Marker == mk, LPS == lps) %>% mutate(Virus = droplevels(Virus))
    cat(sprintf("── %-7s fitting … ", mk))
    m <- fit_bb(pos ~ Virus * Tp + (1 | block), d); cat(attr(m, "engine"), "\n")
    emm_l[[paste(lps,mk)]] <- emm_means(m, ~ Virus | Tp) %>%
      rename(any_of(c(prob = "prob", prob = "response"))) %>%
      mutate(Marker = mk, LPS = lps, pct = 100 * prob)
    ct_l[[paste(lps,mk)]] <- emm_contrasts(m, ~ Virus | Tp, method = "trt.vs.ctrl") %>%
      mutate(Marker = mk, LPS = lps)
  }
}

holm_lps <- function(df) df %>% group_by(LPS, Marker) %>%
  mutate(p_holm = p.adjust(p, "holm")) %>% ungroup() %>%
  mutate(sig = sig_star(p_holm), contrast = gsub("Virus", "", contrast))

emm_df <- bind_rows(emm_l)
ct_df  <- bind_rows(ct_l) %>% holm_lps()
write.csv(emm_df, file.path(OUT, "M3A_emmeans.csv"),   row.names = FALSE)
write.csv(ct_df,  file.path(OUT, "M3A_contrasts.csv"), row.names = FALSE)

cat("\n══ M3A @24h — bystander vs Ctrl log-OR (−LPS vs +LPS) ══\n")
cmp <- ct_df %>% filter(Tp == 24) %>%
  transmute(Marker, contrast, LPS, logOR = round(estimate,2), sig) %>%
  pivot_wider(names_from = LPS, values_from = c(logOR, sig)) %>% arrange(contrast, Marker)
print(as.data.frame(cmp))
cat("\n═══ 03_Paracrine-Activation COMPLETE ═══\n")
