# ==============================================================================
# _Model-Helpers.R  —  Beta-binomial GLMM fitting + emmeans contrast utilities
# Shared by 01_Infection-Kinetics.R, 02_Infected-vs-Bystander-GLMM.R and
# 03_Paracrine-Activation.R.
#   • Counts modelled at event level → gate size enters the likelihood, so small
#     4h NP+ gates are automatically down-weighted.
#   • family = betabinomial absorbs bio-rep overdispersion; block enters as a
#     random intercept (1 | Timepoint:Bio_Rep).
# ==============================================================================

suppressPackageStartupMessages({
  library(glmmTMB); library(emmeans); library(dplyr); library(lme4); library(broom.mixed)
})

# Fit one beta-binomial GLMM. `dat` must contain integer `pos` and `den`.
fit_bb <- function(formula, dat) {
  dat <- dat %>% mutate(.fail = den - pos)
  rhs <- paste(deparse(formula[[3]]), collapse = " ")
  f   <- as.formula(paste("cbind(pos, .fail) ~", rhs))
  m <- tryCatch(
    glmmTMB(f, family = betabinomial, data = dat,
            control = glmmTMBControl(parallel = 1)),
    error = function(e) NULL, warning = function(w) NULL)
  ok <- !is.null(m) && is.finite(logLik(m)) &&
        !is.null(m$sdr) && m$sdr$pdHess
  if (ok) { attr(m, "engine") <- "betabinomial-GLMM"; return(m) }

  # ── fallback: binomial GLMM + OLRE ──
  dat$.olre <- factor(seq_len(nrow(dat)))
  f2 <- as.formula(paste("cbind(pos, .fail) ~", rhs, "+ (1 | .olre)"))
  m2 <- glmer(f2, family = binomial, data = dat,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5)))
  attr(m2, "engine") <- "binomial-GLMM+OLRE"
  m2
}

# emmeans contrast → tidy data frame on the log-odds (link) scale.
#   `specs`   : emmeans rhs (e.g. ~ State | Virus + Timepoint)
#   `method`  : "revpairwise" (default) or a custom contrast list
emm_contrasts <- function(model, specs, method = "revpairwise", by = NULL) {
  emm <- emmeans(model, specs = specs)
  # NB: when `specs` uses the `|` form, the by-grouping is already stored in
  # `emm`; passing by = NULL to contrast() would wipe it, so only override when
  # `by` is explicitly supplied.
  ct  <- if (is.null(by)) contrast(emm, method = method, adjust = "none")
         else             contrast(emm, method = method, by = by, adjust = "none")
  as.data.frame(summary(ct, infer = c(TRUE, TRUE))) %>%
    rename(any_of(c(estimate = "estimate", SE = "SE",
                    CI_low = "asymp.LCL", CI_high = "asymp.UCL",
                    CI_low = "lower.CL",  CI_high = "upper.CL",
                    z = "z.ratio", p = "p.value")))
}

# emmeans marginal means → tidy (response %) for plotting.
emm_means <- function(model, specs) {
  as.data.frame(summary(emmeans(model, specs = specs, type = "response"),
                        infer = c(TRUE, FALSE)))
}

# Holm-adjust a p-value vector within a per-marker family.
holm_by_marker <- function(df, pcol = "p", group = "Marker") {
  df %>% group_by(.data[[group]]) %>%
    mutate(p_holm = p.adjust(.data[[pcol]], method = "holm")) %>%
    ungroup()
}
