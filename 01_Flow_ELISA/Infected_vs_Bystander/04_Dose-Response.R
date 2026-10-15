# ==============================================================================
# 04_Dose-Response.R 
# ==============================================================================

BASE <- "./Flow-ELISA"                                     # data folder (contains DATA/)
HERE <- file.path(BASE, "R-analysis-figure_VS-bystander")  # module output folder
OUT  <- file.path(HERE, "output")
CODE <- "01_Flow_ELISA/Infected_vs_Bystander"              # module scripts (repo-relative)
source(file.path(CODE, "_Style-Bystander.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr);
                                  library(broom); library(readxl); library(stringr) })

kin <- read.csv(file.path(OUT, "infection_kinetics.csv")) %>%
  filter(Virus %in% c("V1","V2"))

# ── (a) whole-well M1-surface activation vs %NP+ ──────────────────────────────
# logit-link quasibinomial GLM per strain: surfM1 ~ pctNP, weight = Live_cell
fit_dose <- function(d, yfrac, wcol) {
  d$.y <- d[[yfrac]]; d$.w <- d[[wcol]]
  m <- glm(.y ~ pctNP, family = quasibinomial("logit"), weights = .w, data = d)
  s <- summary(m)$coefficients
  ci <- suppressMessages(confint.default(m))
  data.frame(slope = s["pctNP","Estimate"], SE = s["pctNP","Std. Error"],
             CI_low = ci["pctNP",1], CI_high = ci["pctNP",2],
             p = s["pctNP","Pr(>|t|)"],
             # Pearson correlation on the raw % scale for reporting
             r = cor(d$pctNP, 100*d$.y))
}

# fitted curve (% scale) on a regular %NP+ grid, with 95% CI
predict_dose <- function(d, yfrac, wcol) {
  d$.y <- d[[yfrac]]; d$.w <- d[[wcol]]
  m <- glm(.y ~ pctNP, family = quasibinomial("logit"), weights = .w, data = d)
  g  <- data.frame(pctNP = seq(min(d$pctNP), max(d$pctNP), length.out = 80))
  pr <- predict(m, g, type = "link", se.fit = TRUE)
  ilink <- family(m)$linkinv
  transform(g, pct = 100*ilink(pr$fit),
               CI_low  = 100*ilink(pr$fit - 1.96*pr$se.fit),
               CI_high = 100*ilink(pr$fit + 1.96*pr$se.fit))
}

surf_slopes <- kin %>%
  mutate(yfrac = surfM1_pos / surfM1_den) %>%
  group_by(Virus) %>%
  group_modify(~ fit_dose(.x %>% mutate(Live_cell = surfM1_den/4),
                          "yfrac", "Live_cell")) %>%
  ungroup() %>% mutate(Marker = "M1-surface (mean)")

# ── (b) per-marker total-live % vs %NP+ ───────────────────────────────────────
cnt <- read_excel(file.path(BASE, "DATA", "FACS_Count_ForR.xlsx"), sheet = 1) %>%
  mutate(Virus = case_when(group %in% c("C_PR8","LPS_PR8") ~ "V1",
                           group %in% c("C_S33","LPS_S33") ~ "V2", TRUE ~ "other"),
         Timepoint = as.integer(str_extract(timepoint, "\\d+")),
         LPS = if_else(grepl("^LPS", group), "Yes", "No"),
         pctNP = 100 * NP / Live_cell) %>%
  filter(Virus %in% c("V1","V2"))

nice <- c(CD40="CD40", CD68="CD68", CD86="CD86", CD206="CD206",
          HLA_I="HLA-I", HLA_II="HLA-II")
per_marker <- map_dfr(names(nice), function(m) {
  cnt %>% transmute(Virus, Timepoint, LPS, pctNP,
                    yfrac = pmin(.data[[m]], Live_cell)/Live_cell,
                    Live_cell, Marker = nice[[m]])
})

marker_slopes <- per_marker %>%
  group_by(Marker, Virus, LPS) %>%
  group_modify(~ fit_dose(.x, "yfrac", "Live_cell")) %>% ungroup()

# surface composite per LPS too
surf_slopes_lps <- kin %>%
  mutate(yfrac = surfM1_pos / surfM1_den, Live_cell = surfM1_den/4) %>%
  group_by(Virus, LPS) %>%
  group_modify(~ fit_dose(.x, "yfrac", "Live_cell")) %>%
  ungroup() %>% mutate(Marker = "M1-surface (mean)")

slopes <- bind_rows(surf_slopes_lps, marker_slopes) %>%
  mutate(sig = sig_star(p))
write.csv(slopes, file.path(OUT, "M4_slopes.csv"), row.names = FALSE)

# points for plotting (per-marker + surface)
pts <- bind_rows(
  per_marker %>% transmute(Virus, Timepoint, LPS, pctNP,
                           pct = 100*yfrac, Marker),
  kin %>% transmute(Virus, Timepoint, LPS, pctNP,
                    pct = surfM1_pct, Marker = "M1-surface (mean)"))
write.csv(pts, file.path(OUT, "M4_points.csv"), row.names = FALSE)

# fitted curves (% scale) for plotting
fit_surf <- kin %>% mutate(yfrac = surfM1_pos/surfM1_den, Live_cell = surfM1_den/4) %>%
  group_by(Virus) %>% group_modify(~ predict_dose(.x, "yfrac", "Live_cell")) %>%
  ungroup() %>% mutate(Marker = "M1-surface (mean)")
fit_marker <- per_marker %>% group_by(Marker, Virus) %>%
  group_modify(~ predict_dose(.x, "yfrac", "Live_cell")) %>% ungroup()
fit_lines <- bind_rows(fit_surf, fit_marker)
write.csv(fit_lines, file.path(OUT, "M4_fit.csv"), row.names = FALSE)

cat("══ M4 — activation vs infection burden (logit slope per %NP+, by LPS) ══\n")
print(slopes %>% transmute(Marker, Virus, LPS, slope = round(slope,3),
        r = round(r,2), p = signif(p,2), sig) %>%
      arrange(Marker, Virus, LPS), n = Inf)
cat("\n═══ 04_Dose-Response COMPLETE ═══\n")
