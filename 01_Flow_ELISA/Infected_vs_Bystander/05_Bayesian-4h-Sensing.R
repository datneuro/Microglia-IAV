# ==============================================================================
# 05_Bayesian-4h-Sensing.R  —  Bayesian cross-check of the EARLY (4h) NP+ vs NP- contrast
#   At 4h the NP+ gate is small, so the asymptotic 95% CI of the frequentist
#   GLMM may be distorted. MCMC does not rely on large-sample asymptotics; the
#   posterior 95% HDI, pd and prior-declared Bayes factors are reported.
# ==============================================================================

BASE <- "./Flow-ELISA"                                     # data folder
HERE <- file.path(BASE, "R-analysis-figure_VS-bystander")  # module output folder
OUT  <- file.path(HERE, "output")
suppressPackageStartupMessages({
  library(dplyr); library(brms); library(bayestestR); library(posterior)
})
set.seed(2026); options(mc.cores = 4)

# Highest-density interval from a posterior draw vector (no extra dependency)
hdi95 <- function(x, p = 0.95) {
  x <- sort(x); n <- length(x); k <- floor(p * n)
  w <- x[(k + 1):n] - x[1:(n - k)]; i <- which.min(w)
  c(lo = x[i], hi = x[i + k])
}

gat <- read.csv(file.path(OUT, "counts_gated_long.csv")) %>%
  filter(Timepoint == 4, LPS == "No") %>%
  mutate(infected = ifelse(State == "NP+", 1L, 0L),
         Virus = factor(Virus), Bio_Rep = factor(Bio_Rep))

markers <- c("HLA-I","CD86","CD68","CD40","HLA-II","CD206")   # HLA-I/CD86 = headline
priors  <- c(set_prior("normal(0,1.5)", class = "b"),
             set_prior("normal(0,5)",   class = "Intercept"))

fit_one <- function(d, base = NULL) {
  if (is.null(base))
    brm(pos | trials(den) ~ infected + Virus + (1 | Bio_Rep),
        data = d, family = beta_binomial(), prior = priors,
        sample_prior = "yes", chains = 4, iter = 3000, warmup = 1500,
        seed = 2026, refresh = 0, silent = 2,
        control = list(adapt_delta = 0.99, max_treedepth = 12))
  else
    update(base, newdata = d, recompile = FALSE, refresh = 0, silent = 2,
           control = list(adapt_delta = 0.99, max_treedepth = 12))
}

# Savage-Dickey BF10 (NP+ effect ≠ 0 vs = 0) under a normal approximation of the
# posterior; prior is the analytic Normal(0, 1.5) we set. Prior-dependent by
# construction — reported with that caveat. Huge BFs are capped to ">100".
savage_dickey <- function(dr, prior_sd = 1.5) {
  prior_d0 <- dnorm(0, 0, prior_sd)
  post_d0  <- dnorm(0, mean(dr), sd(dr))     # normal-approx posterior density at 0
  BF01 <- post_d0 / prior_d0
  c(BF10 = 1 / BF01, BF01 = BF01)
}

res <- list(); base_fit <- NULL
for (mk in markers) {
  cat("── Bayesian fit:", mk, "…\n")
  d  <- gat %>% filter(Marker == mk)
  fit <- fit_one(d, base_fit); if (is.null(base_fit)) base_fit <- fit
  dr  <- as_draws_df(fit)$b_infected
  h   <- hdi95(dr)
  pd  <- max(mean(dr > 0), mean(dr < 0))
  bf  <- savage_dickey(dr)
  res[[mk]] <- data.frame(
    Marker = mk, n_wells = nrow(d)/2,
    med_logOR = median(dr), HDI_low = unname(h["lo"]), HDI_high = unname(h["hi"]),
    pd = pd, BF10 = unname(bf["BF10"]), BF01 = unname(bf["BF01"]))
}

out <- bind_rows(res) %>%
  mutate(across(c(med_logOR, HDI_low, HDI_high), ~ round(.x, 3)),
         pd = round(pd, 4),
         BF10_fmt = ifelse(is.na(BF10), "NA",
                    ifelse(BF10 > 100, ">100", sprintf("%.1f", BF10))),
         BF01_fmt = ifelse(is.na(BF01), "NA",
                    ifelse(BF01 < 0.01, "<0.01", sprintf("%.3f", BF01))))
write.csv(out, file.path(OUT, "bayes_4h_sensing.csv"), row.names = FALSE)

cat("\n══ Bayesian 4h sensing (NP+ vs NP-, −LPS, both strains) ══\n")
print(out %>% transmute(Marker, med_logOR,
        HDI = sprintf("[%.2f, %.2f]", HDI_low, HDI_high),
        pd, BF10 = BF10_fmt, BF01 = BF01_fmt))
cat("\nHDI excludes 0 ⇒ effect statistically robust to small NP+ gate (no asymptotic assumption).\n")
cat("Note: primary timepoint remains 8h; BF values are prior-dependent (Normal(0,1.5)).\n")
cat("\n═══ 05_Bayesian-4h-Sensing COMPLETE ═══\n")
