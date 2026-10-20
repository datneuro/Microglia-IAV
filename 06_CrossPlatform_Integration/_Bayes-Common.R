# ==============================================================================
# _Bayes-Common.R  
# ==============================================================================
suppressPackageStartupMessages({ library(brms); library(bayestestR); library(posterior) })

# ---- MCMC settings -----------------------------------------------------------
CP_SEED        <- 20260615
CP_CHAINS      <- 4
CP_ITER        <- 6000     
CP_WARMUP      <- 1000     
CP_CORES       <- 4
CP_ADAPT_DELTA <- 0.99    
CP_ITER_SENS   <- 4000     
ROPE_RHO       <- 0.1      
CP_DO_BF       <- FALSE    
CP_ESS_MIN     <- 10000    

# ---- weakly-informative priors (responses are z-scaled, coeffs in SD units) ---
# multivariate model: class=b priors must be set per response (resp = a / b)
CP_PRIORS <- c(prior(normal(0, 1.5), class = b, resp = a),
               prior(normal(0, 1.5), class = b, resp = b),
               prior(lkj(2),         class = rescor))

# prior sets for sensitivity (vary both the regression and the correlation prior)
CP_PRIOR_SETS <- list(
  default  = c(prior(normal(0, 1.5), class = b, resp = a),
               prior(normal(0, 1.5), class = b, resp = b),
               prior(lkj(2), class = rescor)),
  broader  = c(prior(normal(0, 3),   class = b, resp = a),
               prior(normal(0, 3),   class = b, resp = b),
               prior(lkj(1), class = rescor)),
  narrower = c(prior(normal(0, 1),   class = b, resp = a),
               prior(normal(0, 1),   class = b, resp = b),
               prior(lkj(4), class = rescor))
)

# ---- model wrapper -----------------------------------------------------------
cp_brm <- function(formula, data, prior = CP_PRIORS, iter = CP_ITER, ...) {
  brm(formula, data = data, prior = prior, seed = CP_SEED,
      chains = CP_CHAINS, iter = iter, warmup = CP_WARMUP, cores = CP_CORES,
      refresh = 0, silent = 2, backend = "rstan",
      control = list(adapt_delta = CP_ADAPT_DELTA), ...)
}

# ---- MCMC convergence QC (mirrors WB mcmc_qc) + rescor-specific ESS -----------
cp_mcmc_qc <- function(fit, label) {
  conv <- posterior::summarise_draws(posterior::as_draws(fit),
            posterior::default_convergence_measures())
  div  <- sum(subset(brms::nuts_params(fit), Parameter == "divergent__")$Value)
  rc   <- conv[conv$variable == "rescor__a__b", , drop = FALSE]
  data.frame(
    model           = label,
    divergent       = div,
    max_rhat        = max(conv$rhat,     na.rm = TRUE),
    min_ess_bulk    = min(conv$ess_bulk, na.rm = TRUE),
    min_ess_tail    = min(conv$ess_tail, na.rm = TRUE),
    rescor_ess_bulk = if (nrow(rc)) rc$ess_bulk else NA_real_,
    rescor_ess_tail = if (nrow(rc)) rc$ess_tail else NA_real_,
    # gate is on the ESTIMAND (rescor / partial correlation): ESS>10k bulk+tail,
    # 0 divergences, Rhat<=1.01. (min_ess over nuisance intercepts reported too.)
    pass = (div == 0 & max(conv$rhat, na.rm = TRUE) <= 1.01 &
            nrow(rc) > 0 &&
            rc$ess_bulk >= CP_ESS_MIN & rc$ess_tail >= CP_ESS_MIN)
  )
}

# ---- describe a rescor draw vector -> tidy row -------------------------------
# median, 95% HDI, pd, ROPE% (proportion of posterior within +/- rope), BF (NA unless enabled)
cp_describe_rho <- function(dr, rope = ROPE_RHO) {
  h  <- bayestestR::hdi(dr, ci = .95)
  rp <- bayestestR::rope(dr, range = c(-rope, rope), ci = 1)
  tibble::tibble(
    rho      = median(dr),
    CI_low   = as.numeric(h$CI_low),
    CI_high  = as.numeric(h$CI_high),
    pd       = as.numeric(bayestestR::p_direction(dr)),
    ROPE_pct = as.numeric(rp$ROPE_Percentage),
    BF       = NA_real_
  )
}

message("[_Bayes-Common] loaded. draws=", CP_CHAINS, "x(", CP_ITER, "-", CP_WARMUP,
        ")=", CP_CHAINS*(CP_ITER-CP_WARMUP), " | ROPE=|", ROPE_RHO, "| | ESS gate=", CP_ESS_MIN)
