# ==============================================================================
# 01_Infection-Kinetics.R  —  %NP+ replication kinetics (Fig 1 backbone)
#   Beta-binomial GLMM on the NP+ gate:  cbind(NP, Live-NP) ~ Virus*LPS*Timepoint
#   Restricted to the four infected groups (C_PR8, C_S33, LPS_PR8, LPS_S33).
#   Question:  which strain replicates faster, and does LPS change burden?
# ==============================================================================

BASE <- "./Flow-ELISA"                                     # data folder (contains DATA/)
HERE <- file.path(BASE, "R-analysis-figure_VS-bystander")  # module output folder
OUT  <- file.path(HERE, "output")
CODE <- "01_Flow_ELISA/Infected_vs_Bystander"              # module scripts (repo-relative)
source(file.path(CODE, "_Model-Helpers.R"))
source(file.path(CODE, "_Style-Bystander.R"))

kin <- read.csv(file.path(OUT, "infection_kinetics.csv"))
inf <- kin %>%
  filter(Virus %in% c("V1","V2")) %>%
  mutate(pos = NP_pos, den = NP_den,
         Virus = factor(Virus, levels = c("V1","V2")),
         LPS   = factor(LPS, levels = c("No","Yes")),
         Tp    = factor(Timepoint, levels = c(4,8,24)),
         block = interaction(Timepoint, Bio_Rep, drop = TRUE))

cat("── Beta-binomial GLMM: %NP+ ~ Virus*LPS*Tp + (1|block) ──\n")
m_np <- fit_bb(pos ~ Virus * LPS * Tp + (1 | block), inf)
cat("engine:", attr(m_np, "engine"), "\n")

# Model-based %NP+ per Virus×LPS×Timepoint
emm_np <- emm_means(m_np, ~ Virus | LPS + Tp) %>%
  rename(any_of(c(pctNP = "prob", pctNP = "response",
                  CI_low = "asymp.LCL", CI_high = "asymp.UCL",
                  CI_low = "lower.CL",  CI_high = "upper.CL"))) %>%
  mutate(across(c(pctNP, CI_low, CI_high), ~ 100 * .x))
write.csv(emm_np, file.path(OUT, "kinetics_emmeans.csv"), row.names = FALSE)

# V1 vs V2 within each LPS × Timepoint (log-odds ratio of being infected)
ct_np <- emm_contrasts(m_np, ~ Virus | LPS + Tp) %>%
  mutate(Marker = "%NP+") %>% holm_by_marker()
write.csv(ct_np, file.path(OUT, "kinetics_contrasts.csv"), row.names = FALSE)

cat("\nModel %NP+ (Virus × LPS × Timepoint):\n")
print(emm_np %>% transmute(Virus, LPS, Tp,
        pctNP = round(pctNP,2), CI = sprintf("[%.2f, %.2f]", CI_low, CI_high)))
cat("\nV1 vs V2 contrast (log-OR, Holm per family):\n")
print(ct_np %>% transmute(LPS, Tp, logOR = round(estimate,2),
        p_holm = signif(p_holm,3), sig = sig_star(p_holm)))

# ── Viability (cell death) — secondary, all groups ────────────────────────────
via <- kin %>% filter(Virus != "PolyIC") %>%
  mutate(pos = Live_cell, den = Single_cell,
         Group = factor(Group, levels = GROUP_LV),
         Tp = factor(Timepoint, levels = c(4,8,24)),
         block = interaction(Timepoint, Bio_Rep, drop = TRUE))
via_sum <- via %>% group_by(Group, Timepoint) %>%
  summarise(viability = mean(viability),
            sd = sd(100*Live_cell/Single_cell), n = n(), .groups = "drop") %>%
  mutate(se = sd/sqrt(n))
write.csv(via_sum, file.path(OUT, "viability_summary.csv"), row.names = FALSE)

cat("\n═══ 01_Infection-Kinetics COMPLETE ═══\n")
