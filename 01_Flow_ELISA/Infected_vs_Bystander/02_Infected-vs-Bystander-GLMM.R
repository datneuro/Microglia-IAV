# ==============================================================================
# 02_Infected-vs-Bystander-GLMM.R 
# ==============================================================================

BASE <- "./Flow-ELISA"                                     # data folder (contains DATA/)
HERE <- file.path(BASE, "R-analysis-figure_VS-bystander")  # module output folder
OUT  <- file.path(HERE, "output")
CODE <- "01_Flow_ELISA/Infected_vs_Bystander"              # module scripts (repo-relative)
source(file.path(CODE, "_Model-Helpers.R"))
source(file.path(CODE, "_Style-Bystander.R"))
suppressPackageStartupMessages({ library(purrr); library(tidyr) })

gated <- read.csv(file.path(OUT, "counts_gated_long.csv")) %>%
  mutate(State = factor(State, levels = c("NP-","NP+")),
         Virus = factor(Virus, levels = c("V1","V2")),
         LPS   = factor(LPS, levels = c("No","Yes")),
         Tp    = factor(Timepoint, levels = c(4,8,24)),
         block = interaction(Timepoint, Bio_Rep, drop = TRUE))

markers <- GATE_MARKERS

# ── fit one marker within one LPS stratum ─────────────────────────────────────
fit_marker <- function(d, mk, lps) {
  m <- fit_bb(pos ~ State * Virus * Tp + (1 | block), d)
  eng <- attr(m, "engine")
  emm <- emm_means(m, ~ State | Virus + Tp) %>%
    rename(any_of(c(prob = "prob", prob = "response"))) %>%
    mutate(Marker = mk, LPS = lps, pct = 100 * prob)
  ct  <- emm_contrasts(m, ~ State | Virus + Tp) %>%       # "NP+ - NP-"
    mutate(Marker = mk, LPS = lps)
  emm_sv <- emmeans(m, ~ State * Virus | Tp)
  int <- contrast(emm_sv,
            method = list("autonomy_V2_minus_V1" =
              c(`NP-,V1` = 1, `NP+,V1` = -1, `NP-,V2` = -1, `NP+,V2` = 1)),
            adjust = "none")
  intd <- as.data.frame(summary(int, infer = c(TRUE, TRUE))) %>%
    rename(any_of(c(CI_low = "asymp.LCL", CI_high = "asymp.UCL",
                    CI_low = "lower.CL", CI_high = "upper.CL",
                    z = "z.ratio", p = "p.value"))) %>%
    mutate(Marker = mk, LPS = lps)
  list(emm = emm, ct = ct, int = intd, eng = data.frame(Marker = mk, LPS = lps, engine = eng))
}

emm_l <- ct_l <- int_l <- eng_l <- list()
for (lps in c("No","Yes")) {
  cat("\n=========== LPS =", lps, "===========\n")
  for (mk in markers) {
    d <- gated %>% filter(Marker == mk, LPS == lps) %>%
      mutate(Virus = droplevels(Virus))
    cat(sprintf("── %-7s fitting … ", mk))
    r <- fit_marker(d, mk, lps); cat(r$eng$engine, "\n")
    k <- paste(lps, mk)
    emm_l[[k]] <- r$emm; ct_l[[k]] <- r$ct; int_l[[k]] <- r$int; eng_l[[k]] <- r$eng
  }
}

# Holm per marker is applied WITHIN each LPS stratum (separate families).
holm_lps <- function(df) df %>% group_by(LPS, Marker) %>%
  mutate(p_holm = p.adjust(p, "holm")) %>% ungroup() %>% mutate(sig = sig_star(p_holm))

emm_df <- bind_rows(emm_l)
ct_df  <- bind_rows(ct_l)  %>% holm_lps()
int_df <- bind_rows(int_l) %>% holm_lps()
eng_df <- bind_rows(eng_l)

write.csv(emm_df, file.path(OUT, "M2_emmeans.csv"),      row.names = FALSE)
write.csv(ct_df,  file.path(OUT, "M2_contrasts.csv"),    row.names = FALSE)
write.csv(int_df, file.path(OUT, "M3B_interaction.csv"), row.names = FALSE)
write.csv(eng_df, file.path(OUT, "M2_model_engine.csv"), row.names = FALSE)

# ── side-by-side comparison of the 24h cross-over: −LPS vs +LPS ────────────────
cat("\n══ M2 @24h — NP+ vs NP- log-OR  (−LPS vs +LPS) ══\n")
cmp <- ct_df %>% filter(Tp == 24) %>%
  transmute(Marker, Virus, LPS, logOR = round(estimate,2), sig) %>%
  pivot_wider(names_from = LPS, values_from = c(logOR, sig)) %>%
  arrange(Virus, Marker)
print(as.data.frame(cmp))

cat("\n══ M3B @24h — strain×autonomy interaction (−LPS vs +LPS) ══\n")
cmpi <- int_df %>% filter(Tp == 24) %>%
  transmute(Marker, LPS, dlogOR = round(estimate,2), sig) %>%
  pivot_wider(names_from = LPS, values_from = c(dlogOR, sig))
print(as.data.frame(cmpi))

cat("\n═══ 02_Infected-vs-Bystander-GLMM COMPLETE ═══\n")
