# ==============================================================================
# 01_Proteasome-Activity-Factorial.R

# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
library(tidyverse)
library(readxl)
library(lubridate)
library(hms)
library(car)
library(emmeans)
library(effectsize)
library(sandwich) 
library(lmtest) 
library(ggsignif)
library(ggpubr)
library(grid)
library(patchwork)
library(openxlsx)

# CONFIG
PROT_DIR <- "Proteasome"  
data_file <- file.path(PROT_DIR, "260515_Virus-HMC3_Proteasome-activity.xlsx")
out_dir_pdf <- file.path(PROT_DIR, "results", "pdf")
out_dir_png <- file.path(PROT_DIR, "results", "png")
data_out_dir <- file.path(PROT_DIR, "processed")

if (!dir.exists(out_dir_pdf)) dir.create(out_dir_pdf, recursive = TRUE)
if (!dir.exists(out_dir_png)) dir.create(out_dir_png, recursive = TRUE)
if (!dir.exists(data_out_dir)) dir.create(data_out_dir, recursive = TRUE)

# GROUPS & PRE-SPECIFIED COMPARISONS
group_order <- c(
  "Naive", "A/PR/8/34", "A/NWS/33",
  "LPS", "LPS+A/PR/8/34", "LPS+A/NWS/33"
)

my_comparisons <- list(
  c("Naive", "A/PR/8/34"), c("Naive", "A/NWS/33"), c("A/PR/8/34", "A/NWS/33"),
  c("LPS", "LPS+A/PR/8/34"), c("LPS", "LPS+A/NWS/33"), c("LPS+A/PR/8/34", "LPS+A/NWS/33"),
  c("Naive", "LPS"), c("Naive", "LPS+A/PR/8/34"), c("Naive", "LPS+A/NWS/33"),
  c("A/PR/8/34", "LPS+A/PR/8/34"), c("A/NWS/33", "LPS+A/NWS/33")
)

# ==============================================================================
# RATE CALCULATION METHOD


RATE_METHOD <- "lm_window"

time_settings <- list(
  "Ac-ANW-AMC"   = c(start = 0, end = 240),
  "Suc-LLVY-AMC" = c(start = 0, end = 150)
)

# ==============================================================================
# COLOR PALETTES
# ==============================================================================
palette_blue <- c(
  "Naive"          = "#D5D8DC", # Platinum
  "A/PR/8/34"      = "#8AB4F8", # Cornflower Blue
  "A/NWS/33"       = "#4285F4", # Bright Blue
  "LPS"            = "#9AA0A6", # Darker Grey
  "LPS+A/PR/8/34"  = "#1967D2", # Cobalt
  "LPS+A/NWS/33"   = "#114B9A"  # Navy
)

palette_grbr <- c(
  "Naive"          = "#D5D8DC", # Light gray
  "A/PR/8/34"      = "#A8C7F0", # Pale blue
  "A/NWS/33"       = "#F2B4B4", # Pale red/pink
  "LPS"            = "#6E7376", # Dark gray
  "LPS+A/PR/8/34"  = "#4285F4", # Bright blue
  "LPS+A/NWS/33"   = "#C0392B"  # Deep red
)

palette_orange <- c(
  "Naive"          = "#FCE5D4", # Very pale orange
  "A/PR/8/34"      = "#F8BFA2", # Pale orange
  "A/NWS/33"       = "#F39872", # Soft orange
  "LPS"            = "#E07B3D", # Medium orange
  "LPS+A/PR/8/34"  = "#C25711", # Deep orange
  "LPS+A/NWS/33"   = "#8B3A07"  # Burnt orange
)

# >>> SWITCH HERE <<<
# Options: "blue", "grbr" (gray-blue-red), "orange"
COLOR_PALETTE <- "grbr"

group_colors <- switch(COLOR_PALETTE,
                       "blue"   = palette_blue,
                       "grbr"   = palette_grbr,
                       "orange" = palette_orange,
                       stop("Unknown COLOR_PALETTE: ", COLOR_PALETTE)
)

lps_pair_colors <- switch(COLOR_PALETTE,
                          "blue"   = c("No LPS" = "#9AA0A6", "+ LPS" = "#114B9A"), # dark grey vs navy
                          "grbr"   = c("No LPS" = "#6E7376", "+ LPS" = "#C0392B"), # dark grey vs deep red
                          "orange" = c("No LPS" = "#F8BFA2", "+ LPS" = "#8B3A07"), # pale vs burnt orange
                          c("No LPS" = group_colors[["Naive"]], "+ LPS" = group_colors[["LPS"]])
)

ADJUST_METHOD <- "BH"

# ==============================================================================
# ROBUST STANDARD ERRORS (HC3 sandwich estimator)
# ==============================================================================
USE_ROBUST_SE <- TRUE

# ==============================================================================
# 2. READ PLATE MAPS & DATA
# ==============================================================================
read_plate_map <- function(filepath, sheet_name, col_name) {
  raw_map <- read_excel(filepath, sheet = sheet_name, col_names = FALSE, range = "A1:L8")
  colnames(raw_map) <- as.character(1:12)
  raw_map$Row <- LETTERS[1:8]
  raw_map %>%
    pivot_longer(cols = -Row, names_to = "Col", values_to = "Value") %>%
    mutate(Well = paste0(Row, Col)) %>%
    select(Well, Value) %>%
    rename(!!sym(col_name) := Value)
}

map_groups     <- read_plate_map(data_file, "map-groups", "Group")
map_inhibit    <- read_plate_map(data_file, "map-inhibit", "Inhibitor")
map_substrates <- read_plate_map(data_file, "map-substrates", "Substrate")

master_map <- map_groups %>%
  left_join(map_inhibit, by = "Well") %>%
  left_join(map_substrates, by = "Well") %>%
  mutate(Group = str_trim(Group), Substrate = str_trim(Substrate)) %>%
  filter(!is.na(Group)) %>%
  mutate(Condition = if_else(!is.na(Inhibitor) & Inhibitor == "MG132", "Background", "Total"))

raw_data <- read_excel(data_file, sheet = "Sheet1")
tidy_data <- raw_data %>%
  mutate(Time_Min = hour(as_hms(Time)) * 60 + minute(as_hms(Time)) + second(as_hms(Time)) / 60) %>%
  select(Time_Min, matches("^[A-H][0-9]+")) %>%
  pivot_longer(-Time_Min, names_to = "Well", values_to = "RFU")

merged_data <- inner_join(tidy_data, master_map, by = "Well")

# ==============================================================================
# 3. RATE & SPECIFIC ACTIVITY
# ------------------------------------------------------------------------------
# ==============================================================================
calc_well_rate <- function(df, sub_type, method = RATE_METHOD) {
  if (is.na(sub_type) || is.null(time_settings[[sub_type]])) {
    return(list(rate = NA_real_, r2 = NA_real_, n_points = 0L))
  }
  limits <- time_settings[[sub_type]]
  df_win <- df %>%
    filter(Time_Min >= limits["start"] & Time_Min <= limits["end"]) %>%
    arrange(Time_Min)
  
  if (nrow(df_win) < 2) return(list(rate = NA_real_, r2 = NA_real_, n_points = nrow(df_win)))
  
  if (method == "delta") {
    dr <- df_win$RFU[nrow(df_win)] - df_win$RFU[1]
    dt <- df_win$Time_Min[nrow(df_win)] - df_win$Time_Min[1]
    rate <- if (dt == 0) NA_real_ else dr / dt
    return(list(rate = rate, r2 = NA_real_, n_points = nrow(df_win)))
  } else if (method == "lm_window") {
    if (nrow(df_win) < 3 || sd(df_win$RFU, na.rm = TRUE) == 0) {
      return(list(rate = NA_real_, r2 = NA_real_, n_points = nrow(df_win)))
    }
    fit <- tryCatch(lm(RFU ~ Time_Min, data = df_win), error = function(e) NULL)
    if (is.null(fit)) return(list(rate = NA_real_, r2 = NA_real_, n_points = nrow(df_win)))
    s <- summary(fit)
    return(list(rate = unname(coef(fit)["Time_Min"]),
                r2 = s$r.squared, n_points = nrow(df_win)))
  } else {
    stop("Unknown RATE_METHOD: ", method)
  }
}

endpoint_rates <- merged_data %>%
  filter(!is.na(Substrate)) %>%
  group_by(Group, Substrate, Condition, Well) %>%
  nest() %>%
  mutate(rate_info = map2(data, Substrate, calc_well_rate)) %>%
  mutate(Rate = map_dbl(rate_info, "rate"),
         R2_fit = map_dbl(rate_info, "r2"),
         N_points = map_int(rate_info, "n_points")) %>%
  select(-data, -rate_info) %>%
  ungroup()

bg_summary <- endpoint_rates %>%
  filter(Condition == "Background") %>%
  group_by(Group, Substrate) %>%
  summarise(Mean_BG = mean(Rate, na.rm = TRUE), .groups = "drop")

specific_activity_df <- endpoint_rates %>%
  filter(Condition == "Total") %>%
  left_join(bg_summary, by = c("Group", "Substrate")) %>%
  mutate(Mean_BG = replace_na(Mean_BG, 0),
         Specific_Activity = Rate - Mean_BG) %>%
  group_by(Group, Substrate) %>%
  arrange(Well) %>%
  mutate(Rep_ID = row_number()) %>%
  ungroup() %>%
  mutate(Group = factor(Group, levels = group_order))

ratio_df <- specific_activity_df %>%
  select(Group, Substrate, Rep_ID, Specific_Activity) %>%
  pivot_wider(names_from = Substrate, values_from = Specific_Activity) %>%
  mutate(Ratio_ANW_LLVY = `Ac-ANW-AMC` / `Suc-LLVY-AMC`) %>%
  drop_na(Ratio_ANW_LLVY) %>%
  mutate(Group = factor(Group, levels = group_order))

write.csv(specific_activity_df,
          file.path(data_out_dir, "Proteasome_SpecificActivity.csv"),
          row.names = FALSE)
write.csv(ratio_df,
          file.path(data_out_dir, "Proteasome_Ratio.csv"),
          row.names = FALSE)

# ==============================================================================
# 4. STATISTICS - PLANNED CONTRASTS
# ==============================================================================
p_to_stars <- function(p) {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", TRUE ~ "ns")
}

# Cohen's d for two independent groups (pooled SD)
cohens_d_pair <- function(data, vcol, gcol, g1, g2) {
  x1 <- data[[vcol]][data[[gcol]] == g1]
  x2 <- data[[vcol]][data[[gcol]] == g2]
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s1 <- sd(x1); s2 <- sd(x2)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (is.na(sp) || sp == 0) return(NA_real_)
  (mean(x2) - mean(x1)) / sp
}

run_stats_planned <- function(data, vcol, gcol, comparisons,
                              adjust_method = "holm",
                              log_transform = FALSE,
                              metric_name = "metric") {
  
  # Keep groups with n >= 2
  grp_n <- data %>%
    group_by(across(all_of(gcol))) %>%
    summarize(n = sum(!is.na(.data[[vcol]])), .groups = "drop")
  keep_groups <- grp_n %>% filter(n >= 2) %>% pull(1)
  data <- data %>%
    filter(.data[[gcol]] %in% keep_groups) %>%
    droplevels()
  
  empty_pw <- data.frame(
    Metric = character(), Group1 = character(), Group2 = character(),
    N1 = integer(), N2 = integer(),
    Mean1 = numeric(), Mean2 = numeric(),
    Estimate = numeric(), SE = numeric(), t_ratio = numeric(),
    df = numeric(), P_raw = numeric(), P_adjusted = numeric(),
    Cohens_d = numeric(), Significance = character()
  )
  
  if (length(unique(data[[gcol]])) < 2) {
    return(list(
      omnibus = data.frame(Metric = metric_name, Summary = "< 2 testable groups"),
      pairwise = empty_pw
    ))
  }
  
  # Optional log-transform
  data$.y <- if (log_transform) log(data[[vcol]]) else data[[vcol]]
  
  # Fit linear model
  formula_str <- as.formula(paste(".y ~", gcol))
  fit <- lm(formula_str, data = data)
  
  # Robust vcov (HC3) for heteroscedasticity-consistent inference.
  # NULL -> emmeans uses classical model-based vcov.
  robust_vcov <- if (USE_ROBUST_SE) {
    tryCatch(sandwich::vcovHC(fit, type = "HC3"),
             error = function(e) NULL)
  } else NULL
  se_type_label <- if (!is.null(robust_vcov)) "HC3 robust SE" else "Classical SE"
  

  anova_tab <- anova(fit)
  om_F  <- anova_tab$`F value`[1]
  om_df1 <- anova_tab$Df[1]
  om_df2 <- anova_tab$Df[2]
  om_p  <- anova_tab$`Pr(>F)`[1]
  
  # Robust omnibus F (Wald test on all group dummies, HC3-adjusted)
  robust_F <- NA_real_; robust_F_p <- NA_real_
  if (!is.null(robust_vcov)) {
    rhs <- names(coef(fit))[-1]  # all non-intercept coefs
    if (length(rhs) > 0) {
      rob_test <- tryCatch(
        car::linearHypothesis(fit, rhs, vcov. = robust_vcov, test = "F"),
        error = function(e) NULL
      )
      if (!is.null(rob_test)) {
        robust_F <- rob_test$F[2]
        robust_F_p <- rob_test$`Pr(>F)`[2]
      }
    }
  }
  

  eta_sq <- tryCatch(
    as.numeric(effectsize::eta_squared(fit, partial = FALSE)$Eta2[1]),
    error = function(e) NA_real_
  )
  
  # Estimated marginal means (with robust vcov if available)
  emm <- if (!is.null(robust_vcov)) {
    emmeans(fit, specs = gcol, vcov. = robust_vcov)
  } else {
    emmeans(fit, specs = gcol)
  }
  available_groups <- as.character(summary(emm)[[gcol]])
  
  # Group means on the original scale
  grp_means <- data %>%
    group_by(across(all_of(gcol))) %>%
    summarize(M = mean(.data[[vcol]], na.rm = TRUE),
              N = sum(!is.na(.data[[vcol]])),
              .groups = "drop")
  
  valid_comps <- comparisons[sapply(comparisons, function(p)
    all(p %in% available_groups))]
  
  if (length(valid_comps) == 0) {
    return(list(
      omnibus = data.frame(
        Metric = metric_name,
        Omnibus_test = paste("Linear model (lm) +", se_type_label),
        F_stat = om_F, df1 = om_df1, df2 = om_df2, Omnibus_p = om_p,
        Robust_F = robust_F, Robust_F_p = robust_F_p,
        Eta2 = eta_sq,
        Adjust_method = adjust_method,
        Log_transformed = log_transform,
        N_comparisons = 0,
        Summary = "no valid comparisons"
      ),
      pairwise = empty_pw
    ))
  }
  
  contrast_list <- lapply(valid_comps, function(pair) {
    vec <- rep(0, length(available_groups))
    vec[available_groups == pair[1]] <- -1
    vec[available_groups == pair[2]] <-  1
    vec
  })
  names(contrast_list) <- sapply(valid_comps, function(p)
    paste(p[2], "vs", p[1]))
  
  ctr_raw <- as.data.frame(summary(
    contrast(emm, method = contrast_list, adjust = "none")
  ))
  ctr_adj <- as.data.frame(summary(
    contrast(emm, method = contrast_list, adjust = adjust_method)
  ))
  
  # Build pairwise output
  pw_list <- lapply(seq_along(valid_comps), function(i) {
    g1 <- valid_comps[[i]][1]
    g2 <- valid_comps[[i]][2]
    n1 <- grp_means$N[grp_means[[gcol]] == g1]
    n2 <- grp_means$N[grp_means[[gcol]] == g2]
    m1 <- grp_means$M[grp_means[[gcol]] == g1]
    m2 <- grp_means$M[grp_means[[gcol]] == g2]
    d_val <- cohens_d_pair(data, vcol, gcol, g1, g2)
    data.frame(
      Metric = metric_name,
      Group1 = g1, Group2 = g2,
      N1 = n1, N2 = n2,
      Mean1 = m1, Mean2 = m2,
      Estimate = ctr_raw$estimate[i],
      SE = ctr_raw$SE[i],
      t_ratio = ctr_raw$t.ratio[i],
      df = ctr_raw$df[i],
      P_raw = ctr_raw$p.value[i],
      P_adjusted = ctr_adj$p.value[i],
      Cohens_d = d_val,
      Significance = p_to_stars(ctr_adj$p.value[i]),
      stringsAsFactors = FALSE
    )
  })
  pw <- bind_rows(pw_list)
  
  resid_vec <- residuals(fit)
  shap_resid <- tryCatch(
    if (length(resid_vec) >= 3 && length(resid_vec) <= 5000 && sd(resid_vec) > 0)
      shapiro.test(resid_vec)$p.value else NA_real_,
    error = function(e) NA_real_
  )
  lev_p <- tryCatch(
    car::leveneTest(fit)$`Pr(>F)`[1],
    error = function(e) NA_real_
  )
  
  # Variance homogeneity diagnostic (with n=3, Levene is unreliable;
  # SD ratio is a simple practical check. Conservatively flag if > 4.)
  grp_sds <- tapply(data[[vcol]], data[[gcol]], function(x) sd(x, na.rm = TRUE))
  grp_sds <- grp_sds[!is.na(grp_sds) & grp_sds > 0]
  sd_ratio <- if (length(grp_sds) >= 2) max(grp_sds) / min(grp_sds) else NA_real_
  min_n <- min(grp_n$n, na.rm = TRUE)
  
  omnibus <- data.frame(
    Metric = metric_name,
    Omnibus_test = paste("Linear model (lm) +", se_type_label, "+ planned contrasts"),
    F_stat = om_F, df1 = om_df1, df2 = om_df2,
    Omnibus_p = om_p,
    Robust_F = robust_F,
    Robust_F_p = robust_F_p,
    SE_type = se_type_label,
    Eta2 = eta_sq,
    Residual_Shapiro_p = shap_resid,
    Levene_p = lev_p,
    SD_ratio_max_min = sd_ratio,
    Variance_flag = if (is.na(sd_ratio)) "" else if (sd_ratio > 4) "HIGH (>4)" else "OK",
    Min_n_per_group = min_n,
    Adjust_method = adjust_method,
    Log_transformed = log_transform,
    N_comparisons = length(valid_comps),
    Summary = sprintf("LM F(%d,%d)=%.3f, p=%.3g; %s; %d planned contrasts adj=%s",
                      om_df1, om_df2, om_F, om_p,
                      se_type_label, length(valid_comps), adjust_method),
    stringsAsFactors = FALSE
  )
  
  list(omnibus = omnibus, pairwise = pw, fit = fit, emm = emm)
}

# ==============================================================================
# 5. FACTORIAL ANALYSIS 
add_factorial_cols <- function(df, group_col = "Group") {
  df %>%
    mutate(
      LPS = if_else(str_detect(.data[[group_col]], "LPS"), "Yes", "No"),
      LPS = factor(LPS, levels = c("No", "Yes")),
      Virus = case_when(
        str_detect(.data[[group_col]], "A/PR/8/34") ~ "PR8",
        str_detect(.data[[group_col]], "A/NWS/33") ~ "NWS",
        TRUE ~ "None"
      ),
      Virus = factor(Virus, levels = c("None", "PR8", "NWS"))
    )
}

run_factorial <- function(data, vcol, log_transform = FALSE, metric_name = "metric") {
  d <- add_factorial_cols(data, "Group")
  d$.y <- if (log_transform) log(d[[vcol]]) else d[[vcol]]
  
  # Type III SS (proper for unbalanced factorial)
  old_contrasts <- options("contrasts")
  options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old_contrasts))
  
  fit <- lm(.y ~ LPS * Virus, data = d)
  
  # Robust vcov (HC3) for heteroscedasticity-consistent inference
  robust_vcov <- if (USE_ROBUST_SE) {
    tryCatch(sandwich::vcovHC(fit, type = "HC3"),
             error = function(e) NULL)
  } else NULL
  se_type_label <- if (!is.null(robust_vcov)) "HC3 robust SE" else "Classical SE"
  
  # Anova with robust SE via white.adjust (car uses HC3 by default when white.adjust=TRUE)
  aov_tab <- if (!is.null(robust_vcov)) {
    car::Anova(fit, type = "III", white.adjust = "hc3")
  } else {
    car::Anova(fit, type = "III")
  }
  

  eta_partial <- tryCatch(
    effectsize::eta_squared(fit, partial = TRUE, ci = NULL),
    error = function(e) NULL
  )
  

  get_col <- function(tab, ...) {
    for (nm in c(...)) {
      if (nm %in% colnames(tab)) return(tab[[nm]])
    }
    rep(NA_real_, nrow(tab))
  }
  
  terms_of_interest <- c("LPS", "Virus", "LPS:Virus")
  rows <- lapply(terms_of_interest, function(term) {
    if (!term %in% rownames(aov_tab)) return(NULL)
    idx <- which(rownames(aov_tab) == term)
    eta_val <- NA_real_
    if (!is.null(eta_partial)) {
      eta_row <- eta_partial[eta_partial$Parameter == term, ]
      if (nrow(eta_row) > 0) eta_val <- eta_row$Eta2_partial[1]
    }
    data.frame(
      Metric = metric_name,
      Term = term,
      SS = get_col(aov_tab, "Sum Sq")[idx],   # NA if HC3 (no SS)
      df = get_col(aov_tab, "Df")[idx],
      F_stat = get_col(aov_tab, "F", "F value")[idx],
      P_Value = get_col(aov_tab, "Pr(>F)")[idx],
      Partial_Eta2 = eta_val,
      SE_type = se_type_label,
      Significance = p_to_stars(get_col(aov_tab, "Pr(>F)")[idx]),
      stringsAsFactors = FALSE
    )
  })
  factorial_tab <- bind_rows(rows)
  
  # Marginal means with robust vcov for SE & CIs
  emm_lps <- if (!is.null(robust_vcov)) emmeans(fit, ~ LPS, vcov. = robust_vcov)
  else emmeans(fit, ~ LPS)
  emm_vir <- if (!is.null(robust_vcov)) emmeans(fit, ~ Virus, vcov. = robust_vcov)
  else emmeans(fit, ~ Virus)
  emm_int <- if (!is.null(robust_vcov)) emmeans(fit, ~ LPS * Virus, vcov. = robust_vcov)
  else emmeans(fit, ~ LPS * Virus)
  
  emm_lps <- as.data.frame(summary(emm_lps))
  emm_vir <- as.data.frame(summary(emm_vir))
  emm_int <- as.data.frame(summary(emm_int))
  
  margin_tab <- bind_rows(
    emm_lps %>% transmute(Metric = metric_name, Factor = "LPS",
                          Level = as.character(LPS),
                          EMM = emmean, SE = SE,
                          LCL = lower.CL, UCL = upper.CL),
    emm_vir %>% transmute(Metric = metric_name, Factor = "Virus",
                          Level = as.character(Virus),
                          EMM = emmean, SE = SE,
                          LCL = lower.CL, UCL = upper.CL),
    emm_int %>% transmute(Metric = metric_name, Factor = "LPS:Virus",
                          Level = paste(LPS, Virus, sep = "_"),
                          EMM = emmean, SE = SE,
                          LCL = lower.CL, UCL = upper.CL)
  )
  
  list(factorial = factorial_tab, marginals = margin_tab, fit = fit)
}

# ==============================================================================
# 6. RUN STATISTICS
# ==============================================================================
df_anw  <- specific_activity_df %>% filter(Substrate == "Ac-ANW-AMC")
df_llvy <- specific_activity_df %>% filter(Substrate == "Suc-LLVY-AMC")

# Planned contrasts
stat_anw   <- run_stats_planned(df_anw, "Specific_Activity", "Group",
                                my_comparisons, adjust_method = ADJUST_METHOD,
                                log_transform = FALSE, metric_name = "Ac-ANW-AMC")
stat_llvy  <- run_stats_planned(df_llvy, "Specific_Activity", "Group",
                                my_comparisons, adjust_method = ADJUST_METHOD,
                                log_transform = FALSE, metric_name = "Suc-LLVY-AMC")
stat_ratio <- run_stats_planned(ratio_df, "Ratio_ANW_LLVY", "Group",
                                my_comparisons, adjust_method = ADJUST_METHOD,
                                log_transform = TRUE,
                                metric_name = "Ratio_ANW_LLVY (log)")

# Factorial
fac_anw   <- run_factorial(df_anw, "Specific_Activity", log_transform = FALSE,
                           metric_name = "Ac-ANW-AMC")
fac_llvy  <- run_factorial(df_llvy, "Specific_Activity", log_transform = FALSE,
                           metric_name = "Suc-LLVY-AMC")
fac_ratio <- run_factorial(ratio_df, "Ratio_ANW_LLVY", log_transform = TRUE,
                           metric_name = "Ratio_ANW_LLVY (log)")

# ==============================================================================
# 7. EXPORT FULL STATISTICAL LOG (multi-sheet Excel)
# ==============================================================================
omnibus_log <- bind_rows(stat_anw$omnibus, stat_llvy$omnibus, stat_ratio$omnibus)
pairwise_log <- bind_rows(stat_anw$pairwise, stat_llvy$pairwise, stat_ratio$pairwise)
factorial_log <- bind_rows(fac_anw$factorial, fac_llvy$factorial, fac_ratio$factorial)
marginals_log <- bind_rows(fac_anw$marginals, fac_llvy$marginals, fac_ratio$marginals)

# R^2 QC: per-well fit quality from lm(RFU ~ Time_Min) within linear window
r2_qc_log <- endpoint_rates %>%
  filter(Condition == "Total") %>%
  transmute(Group, Substrate, Well, Rate, R2 = R2_fit, N_points) %>%
  arrange(Substrate, Group, Well)
r2_qc_summary <- r2_qc_log %>%
  group_by(Substrate) %>%
  summarize(
    N_wells = n(),
    Median_R2 = median(R2, na.rm = TRUE),
    Min_R2 = min(R2, na.rm = TRUE),
    Wells_R2_below_0.95 = sum(R2 < 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Group n's
group_n_log <- specific_activity_df %>%
  group_by(Substrate, Group) %>%
  summarize(N = sum(!is.na(Specific_Activity)),
            Mean = mean(Specific_Activity, na.rm = TRUE),
            SD = sd(Specific_Activity, na.rm = TRUE),
            SEM = SD / sqrt(N),
            .groups = "drop") %>%
  bind_rows(
    ratio_df %>%
      group_by(Group) %>%
      summarize(Substrate = "Ratio_ANW_LLVY",
                N = sum(!is.na(Ratio_ANW_LLVY)),
                Mean = mean(Ratio_ANW_LLVY, na.rm = TRUE),
                SD = sd(Ratio_ANW_LLVY, na.rm = TRUE),
                SEM = SD / sqrt(N),
                .groups = "drop") %>%
      select(Substrate, Group, N, Mean, SD, SEM)
  )

# Write multi-sheet xlsx
wb <- createWorkbook()
addWorksheet(wb, "1_Group_Summary"); writeData(wb, "1_Group_Summary", group_n_log)
addWorksheet(wb, "2_Omnibus_LM"); writeData(wb, "2_Omnibus_LM", omnibus_log)
addWorksheet(wb, "3_Planned_Contrasts"); writeData(wb, "3_Planned_Contrasts", pairwise_log)
addWorksheet(wb, "4_Factorial_ANOVA"); writeData(wb, "4_Factorial_ANOVA", factorial_log)
addWorksheet(wb, "5_Marginal_Means"); writeData(wb, "5_Marginal_Means", marginals_log)
addWorksheet(wb, "6_R2_QC_PerWell"); writeData(wb, "6_R2_QC_PerWell", r2_qc_log)
addWorksheet(wb, "7_R2_QC_Summary"); writeData(wb, "7_R2_QC_Summary", r2_qc_summary)

# Methods sheet
methods_text <- data.frame(
  Section = c(
    "Design",
    "Rate calculation",
    "Linear windows",
    "Specific activity",
    "Ratio (immunoproteasome shift)",
    "Inference framework",
    "Robust SE",
    "Planned contrasts",
    "Multiple testing",
    "Factorial decomposition",
    "Effect sizes",
    "Variance diagnostics",
    "Software"
  ),
  Description = c(
    "6 groups (Naive, A/PR/8/34, A/NWS/33, LPS, LPS+A/PR/8/34, LPS+A/NWS/33), n=3/group, 2 substrates",
    sprintf("Method = '%s'. Rate per well = slope of lm(RFU ~ Time_Min) within substrate-specific linear window.",
            RATE_METHOD),
    sprintf("Linear windows (auto-detected by max-R^2 criterion, see 02_Rate-Estimation-QC.R): Ac-ANW-AMC [%g, %g] min; Suc-LLVY-AMC [%g, %g] min",
            time_settings[["Ac-ANW-AMC"]]["start"], time_settings[["Ac-ANW-AMC"]]["end"],
            time_settings[["Suc-LLVY-AMC"]]["start"], time_settings[["Suc-LLVY-AMC"]]["end"]),
    "Specific Activity = well rate (Total) minus mean rate of MG132-inhibited wells (Background) within the same Group and Substrate",
    "Ratio = Ac-ANW-AMC / Suc-LLVY-AMC; log-transformed before testing (symmetric, additive)",
    "Linear model lm(y ~ Group) + estimated marginal means (emmeans). Same lm for omnibus F and contrasts; SEs may be classical or robust (see below).",
    sprintf("USE_ROBUST_SE = %s. When TRUE, HC3 sandwich SE (MacKinnon-White 1985) is used via sandwich::vcovHC. Heteroscedasticity-consistent: valid under unequal group variances while preserving point estimates and df. Recommended for n=3/group with SD ratio > 4.",
            USE_ROBUST_SE),
    sprintf("%d pre-specified pairs (planned contrasts via emmeans).", length(my_comparisons)),
    sprintf("p-values adjusted within metric using '%s' across the %d planned pairs",
            ADJUST_METHOD, length(my_comparisons)),
    sprintf("2x3 factorial lm(y ~ LPS * Virus) with Type III SS (contr.sum). %s. Tests LPS main, Virus main, LPS:Virus interaction.",
            if (USE_ROBUST_SE) "HC3-adjusted F via car::Anova(white.adjust='hc3')" else "Classical F"),
    "Pairwise: Cohen's d (pooled SD, raw scale). Omnibus: eta^2. Factorial: partial eta^2.",
    "Per-group SD ratio (max/min) flagged 'HIGH' if > 4. Levene's test and Shapiro on residuals reported but not used to drive analysis (unreliable at n=3).",
    sprintf("R %s; emmeans, car, effectsize, sandwich, lmtest, openxlsx",
            paste(R.version$major, R.version$minor, sep = "."))
  )
)
addWorksheet(wb, "0_Methods"); writeData(wb, "0_Methods", methods_text)
worksheetOrder(wb) <- c(which(names(wb) == "0_Methods"),
                        which(names(wb) != "0_Methods"))

stats_path <- file.path(data_out_dir, "Proteasome_Statistics_Log.xlsx")
saveWorkbook(wb, stats_path, overwrite = TRUE)
message("Stats log saved: ", stats_path)

# ==============================================================================
# 8. PLOTTING - same style as the Western-blot and in vivo flow figures
# ==============================================================================
# Helper: significant pairs from a stats result (adjusted p)
get_sig <- function(res, p_col = "P_adjusted") {
  if (is.null(res) || nrow(res$pairwise) == 0) {
    return(list(pairs = list(), labels = character()))
  }
  sig <- res$pairwise %>% filter(.data[[p_col]] < 0.05)
  if (nrow(sig) == 0) {
    return(list(pairs = list(), labels = character()))
  }
  list(
    pairs = lapply(seq_len(nrow(sig)), function(i) c(sig$Group1[i], sig$Group2[i])),
    labels = sig$Significance
  )
}

# Bar style: col alpha 0.6, errorbar, jitter shape 21 white fill
create_bar <- function(df, xcol, ycol, cmap, ylab,
                       sig_p = list(), sig_l = character(),
                       title_text = NULL) {
  sdf <- df %>%
    group_by(across(all_of(xcol))) %>%
    summarise(
      m = mean(.data[[ycol]], na.rm = TRUE),
      s = sd(.data[[ycol]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[ycol]]))),
      .groups = "drop"
    )
  
  mx <- max(df[[ycol]], na.rm = TRUE)
  ns <- length(sig_p)
  hr <- if (ns > 0) 0.5 + 0.08 * ns else 1.15
  
  p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_col(data = sdf, aes(x = .data[[xcol]], y = m, fill = .data[[xcol]]),
             alpha = 0.6, width = 0.7, color = NA) +
    geom_errorbar(data = sdf,
                  aes(x = .data[[xcol]], y = m,
                      ymin = pmax(m - s, 0), ymax = m + s,
                      color = .data[[xcol]]),
                  width = 0.22, linewidth = 0.6) +
    geom_jitter(aes(color = .data[[xcol]]),
                shape = 21, fill = "white", size = 3, stroke = 0.8,
                alpha = 0.85, width = 0.15) +
    scale_fill_manual(values = cmap, drop = FALSE) +
    scale_color_manual(values = cmap, drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.5), drop = FALSE) +
    expand_limits(y = c(0, mx * hr)) +
    labs(y = ylab, x = NULL, title = title_text) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.margin = margin(10, 5, 0, 5),
      legend.position = "none"
    )
  
  if (ns > 0) {
    p <- p + ggsignif::geom_signif(
      comparisons = sig_p, annotations = sig_l,
      step_increase = 0.075, tip_length = 0.01,
      vjust = 0.4, textsize = 3.5
    )
  }
  p + coord_cartesian(clip = "off")
}

# Stratification axis: separates "No LPS" and "+LPS" with segments + labels
build_strat_axis <- function() {
  ggplot() +
    annotate("segment", x = 1, xend = 3, y = 0.8, yend = 0.8, linewidth = 1.2) +
    annotate("text", x = 2, y = 0.4, label = "No LPS",
             size = 4.5, fontface = "bold") +
    annotate("segment", x = 4, xend = 6, y = 0.8, yend = 0.8, linewidth = 1.2) +
    annotate("text", x = 5, y = 0.4, label = "+ LPS",
             size = 4.5, fontface = "bold") +
    scale_x_continuous(limits = c(0.5, 6.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.2), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 5, 15, 5))
}

# Title box in the draw_box style
draw_box <- function(label, fill_c = "#FBEFE3", x_pos = 0.5) {
  ggplot() +
    geom_label(aes(x = x_pos, y = 0.5, label = label),
               fill = fill_c, color = "black",
               size = 5.5, fontface = "bold",
               label.padding = unit(0.8, "lines"),
               label.size = 0.3) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1)
}

# Align the stratification axis with the plot panel
align_marks <- function(g_plot, g_marks) {
  pp <- g_plot$layout[g_plot$layout$name == "panel", ]
  pm <- g_marks$layout[g_marks$layout$name == "panel", ]
  lp <- pp$l; lm <- pm$l
  g_marks$widths[1:(lm - 1)] <- unit(0, "cm")
  g_marks$widths[1] <- sum(g_plot$widths[1:(lp - 1)])
  g_marks$widths[lm] <- g_plot$widths[lp]
  g_marks$widths[(lm + 1):length(g_marks$widths)] <- unit(0, "cm")
  g_marks$widths[length(g_marks$widths)] <- sum(g_plot$widths[(lp + 1):length(g_plot$widths)])
  g_marks
}

# Shared legend
shared_legend <- ggpubr::get_legend(
  ggplot(data.frame(group = factor(group_order, levels = group_order), val = 1),
         aes(x = group, y = val, fill = group)) +
    geom_point(shape = 22, size = 6, color = "transparent") +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
)

# ==============================================================================
# 9. BUILD 3-PANEL FIGURE
# ==============================================================================
# Build each plot
sig_anw   <- get_sig(stat_anw)
sig_llvy  <- get_sig(stat_llvy)
sig_ratio <- get_sig(stat_ratio)

p_anw <- create_bar(df_anw, "Group", "Specific_Activity", group_colors,
                    ylab = "Ac-ANW-AMC (RFU/min)",
                    sig_p = sig_anw$pairs, sig_l = sig_anw$labels,
                    title_text = NULL)
p_llvy <- create_bar(df_llvy, "Group", "Specific_Activity", group_colors,
                     ylab = "Suc-LLVY-AMC (RFU/min)",
                     sig_p = sig_llvy$pairs, sig_l = sig_llvy$labels,
                     title_text = NULL)
p_ratio <- create_bar(ratio_df, "Group", "Ratio_ANW_LLVY", group_colors,
                      ylab = "Ac-ANW-AMC / Suc-LLVY-AMC (ratio)",
                      sig_p = sig_ratio$pairs, sig_l = sig_ratio$labels,
                      title_text = NULL)

formal_title <- "Proteasome and Immunoproteasome Profiling"

# Render PDF (8.27 x 11.69, portrait A4, 3 plots stacked with stratification axis)
render_figure <- function(out_path, dev_fun) {
  dev_fun(out_path, width = 8.27, height = 11.69)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 8, ncol = 1,
    heights = c(0.8,   # title box
                3.0,   # plot ANW
                0.6,   # strat axis
                3.0,   # plot LLVY
                0.6,   # strat axis
                3.0,   # plot ratio
                0.6,   # strat axis
                0.8)   # legend
  )))
  
  # Title box
  print(draw_box(formal_title, fill_c = "#FBEFE3"),
        vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  
  # Plots + strat axes
  plots <- list(p_anw, p_llvy, p_ratio)
  plot_rows <- c(2, 4, 6)
  mark_rows <- c(3, 5, 7)
  for (k in seq_along(plots)) {
    gp <- ggplotGrob(plots[[k]])
    pushViewport(viewport(layout.pos.row = plot_rows[k], layout.pos.col = 1))
    grid.draw(gp)
    upViewport()
    gm <- align_marks(gp, ggplotGrob(build_strat_axis()))
    pushViewport(viewport(layout.pos.row = mark_rows[k], layout.pos.col = 1))
    grid.draw(gm)
    upViewport()
  }
  
  # Legend
  pushViewport(viewport(layout.pos.row = 8, layout.pos.col = 1))
  grid.draw(shared_legend)
  upViewport()
  upViewport()
  dev.off()
}

pdf_path <- file.path(out_dir_pdf,
                      "Proteasome_Profile.pdf")
png_path <- file.path(out_dir_png,
                      "Proteasome_Profile.png")

render_figure(pdf_path, function(p, ...) cairo_pdf(p, ..., onefile = TRUE))
render_figure(png_path, function(p, ...) png(p, ..., units = "in", res = 300, type = "cairo"))

message("Figure 1 (bar + planned contrasts) saved.")

# ==============================================================================
# 10. KINETIC PLOT (RFU vs Time) - rate visualization
# ==============================================================================

# Mean ± 95% CI at each timepoint
kinetic_stats <- merged_data %>%
  filter(!is.na(Substrate), Condition == "Total") %>%
  group_by(Group, Substrate, Time_Min) %>%
  summarise(
    Mean_RFU = mean(RFU, na.rm = TRUE),
    SD_RFU = sd(RFU, na.rm = TRUE),
    n = sum(!is.na(RFU)),
    .groups = "drop"
  ) %>%
  mutate(
    SEM = SD_RFU / sqrt(n),
    # 95% CI (wide band)
    t_val = qt(0.975, df = if_else(n > 1, n - 1, 1L)),
    CI95_Margin = replace_na(SEM * t_val, 0),
    # Narrower band: SEM only (recommended for visual clarity with 6 overlapping curves)
    SEM_Margin = replace_na(SEM, 0),
    Upper_CI = Mean_RFU + CI95_Margin,
    Lower_CI = Mean_RFU - CI95_Margin,
    Upper_SEM = Mean_RFU + SEM_Margin,
    Lower_SEM = Mean_RFU - SEM_Margin,
    Group = factor(Group, levels = group_order)
  )

# Add factorial cols for stratified faceting (LPS strata)
kinetic_stats <- kinetic_stats %>%
  mutate(
    LPS_strata = if_else(str_detect(as.character(Group), "LPS"),
                         "+ LPS", "No LPS"),
    LPS_strata = factor(LPS_strata, levels = c("No LPS", "+ LPS")),
    # Virus_label for cleaner facet headers if needed
    Virus_label = case_when(
      str_detect(as.character(Group), "A/PR/8/34") ~ "A/PR/8/34",
      str_detect(as.character(Group), "A/NWS/33") ~ "A/NWS/33",
      TRUE ~ "No virus"
    )
  )

# Save kinetic stats to CSV for reference
write.csv(kinetic_stats,
          file.path(data_out_dir, "Proteasome_KineticStats.csv"),
          row.names = FALSE)

create_kinetic_plot <- function(df, substrate_name, title_text = NULL,
                                band_type = c("SEM", "CI95"),
                                facet_by_lps = TRUE) {
  band_type <- match.arg(band_type)
  sub_df <- df %>% filter(Substrate == substrate_name)
  plot_title <- if (is.null(title_text)) substrate_name else title_text
  
  if (band_type == "SEM") {
    ymin_col <- "Lower_SEM"; ymax_col <- "Upper_SEM"
    band_lbl <- "Mean \u00B1 SEM"
  } else {
    ymin_col <- "Lower_CI"; ymax_col <- "Upper_CI"
    band_lbl <- "Mean \u00B1 95% CI"
  }
  
  p <- ggplot(sub_df, aes(x = Time_Min, y = Mean_RFU,
                          color = Group, fill = Group, group = Group)) +
    geom_ribbon(aes(ymin = .data[[ymin_col]], ymax = .data[[ymax_col]]),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = plot_title,
      x = "Time (min)",
      y = paste0(substrate_name, " — RFU (", band_lbl, ")")
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.margin = margin(10, 10, 5, 10),
      legend.position = "none",
      strip.background = element_rect(fill = "#F4F4F4", color = NA),
      strip.text = element_text(face = "bold", size = 10)
    )
  
  if (facet_by_lps) {
    p <- p + facet_wrap(~ LPS_strata, nrow = 1, scales = "fixed")
  }
  p
}

# Safe NULL-coalesce helper (used later in interaction plot)
`%||%` <- function(a, b) if (is.null(a)) b else a

# >>> KINETIC BAND TYPE: "SEM" (recommended, narrower) or "CI95" (wider) <<<
KINETIC_BAND <- "SEM"

p_kin_anw  <- create_kinetic_plot(kinetic_stats, "Ac-ANW-AMC",
                                  title_text = "Immunoproteasome activity",
                                  band_type = KINETIC_BAND, facet_by_lps = TRUE)
p_kin_llvy <- create_kinetic_plot(kinetic_stats, "Suc-LLVY-AMC",
                                  title_text = "Standard proteasome activity",
                                  band_type = KINETIC_BAND, facet_by_lps = TRUE)

# ==============================================================================
# 11. INTERACTION PLOT (Factorial: LPS x Virus) - visualizes sheets 4 & 5
# ==============================================================================

# Interaction EMMs for each metric from the factorial fit
get_interaction_emm <- function(fac_result, on_log_scale = FALSE) {
  int_df <- fac_result$marginals %>% filter(Factor == "LPS:Virus")
  if (nrow(int_df) == 0) return(NULL)
  int_df <- int_df %>%
    separate(Level, into = c("LPS", "Virus"), sep = "_", remove = FALSE) %>%
    mutate(
      LPS = factor(LPS, levels = c("No", "Yes"),
                   labels = c("No LPS", "+ LPS")),
      # Full virus names (not abbreviations)
      Virus = factor(Virus,
                     levels = c("None", "PR8", "NWS"),
                     labels = c("No virus", "A/PR/8/34", "A/NWS/33"))
    )
  # If model fit on log scale, back-transform for display
  if (on_log_scale) {
    int_df <- int_df %>%
      mutate(EMM_disp = exp(EMM),
             LCL_disp = exp(LCL),
             UCL_disp = exp(UCL),
             SE_disp = (UCL_disp - LCL_disp) / (2 * 1.96))  # approx SE on display scale
  } else {
    int_df <- int_df %>%
      mutate(EMM_disp = EMM,
             LCL_disp = LCL,
             UCL_disp = UCL,
             SE_disp = SE)
  }
  # For ribbon (band), need Virus as numeric for smooth fill
  int_df <- int_df %>%
    mutate(Virus_num = as.numeric(Virus))
  int_df
}

# Annotation text from the factorial table
make_factorial_annot <- function(fac_result) {
  f_tab <- fac_result$factorial
  if (nrow(f_tab) == 0) return("")
  lines <- vapply(seq_len(nrow(f_tab)), function(i) {
    term_lbl <- f_tab$Term[i]
    # Compact term labels
    term_lbl <- gsub("LPS:Virus", "LPS \u00D7 Virus", term_lbl)
    p_val <- f_tab$P_Value[i]
    p_fmt <- if (is.na(p_val)) "NA" else if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
    eta_val <- f_tab$Partial_Eta2[i]
    eta_fmt <- if (is.na(eta_val)) "" else sprintf(", \u03B7\u00B2p=%.2f", eta_val)
    sprintf("%s: F(%g,%g)=%.2f, p=%s%s%s",
            term_lbl, f_tab$df[i],
            fac_result$fit$df.residual,
            f_tab$F_stat[i], p_fmt, eta_fmt,
            ifelse(f_tab$Significance[i] == "ns" || f_tab$Significance[i] == "",
                   "", paste0(" ", f_tab$Significance[i])))
  }, character(1))
  paste(lines, collapse = "\n")
}

# Y-axis labels per metric (substrate name)
ylab_map <- list(
  "Ac-ANW-AMC" = "Ac-ANW-AMC (RFU/min)",
  "Suc-LLVY-AMC" = "Suc-LLVY-AMC (RFU/min)",
  "Ratio_ANW_LLVY (log)" = "Ac-ANW-AMC / Suc-LLVY-AMC (ratio)"
)

create_interaction_plot <- function(fac_result, on_log_scale, title_text,
                                    lps_colors = lps_pair_colors) {
  int_df <- get_interaction_emm(fac_result, on_log_scale = on_log_scale)
  if (is.null(int_df)) return(NULL)
  
  annot <- make_factorial_annot(fac_result)
  
  # Positioning: bring annotation up clearly above the data ribbon
  y_max <- max(int_df$UCL_disp, na.rm = TRUE)
  y_min <- min(int_df$LCL_disp, na.rm = TRUE)
  y_range <- y_max - y_min
  # Annotation block sits in dedicated headroom (top ~25-30% of plot)
  y_annot     <- y_max + y_range * 0.40   # top of text
  y_plot_top  <- y_annot + y_range * 0.08
  y_plot_bot  <- y_min - y_range * 0.05
  
  ylab <- ylab_map[[fac_result$factorial$Metric[1]]] %||% "Value"
  
  ggplot(int_df, aes(x = Virus_num, y = EMM_disp,
                     color = LPS, fill = LPS, group = LPS)) +
    # Soft ribbon (option 1: bandwidth)
    geom_ribbon(aes(ymin = LCL_disp, ymax = UCL_disp),
                alpha = 0.18, color = NA) +
    # Slim trend line
    geom_line(linewidth = 0.9) +
    # Slim capped errorbar — visual precision marker (option 2 refinement)
    geom_errorbar(aes(ymin = LCL_disp, ymax = UCL_disp),
                  width = 0.08, linewidth = 0.4, alpha = 0.7) +
    # Round filled point (cleaner than shape=21)
    geom_point(size = 3.0, shape = 16) +
    # x-axis: discrete labels at integer positions
    scale_x_continuous(
      breaks = unique(int_df$Virus_num),
      labels = levels(int_df$Virus),
      expand = expansion(add = 0.3)
    ) +
    scale_color_manual(values = lps_colors, name = NULL) +
    scale_fill_manual(values = lps_colors, name = NULL) +
    # Annotation block — clearly above curves. Size 2.4 fits A4-landscape
    # narrow panels (~2.5 inch wide). Mono font keeps columns aligned.
    annotate("text", x = min(int_df$Virus_num) - 0.25, y = y_annot,
             label = annot, hjust = 0, vjust = 1,
             size = 2.4, family = "mono",
             lineheight = 1.1) +
    coord_cartesian(ylim = c(y_plot_bot, y_plot_top), clip = "off") +
    labs(title = title_text, x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(
      axis.title.y = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 9, face = "plain"),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "top",
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(8, 8, 4, 8)
    )
}

p_int_anw   <- create_interaction_plot(fac_anw, on_log_scale = FALSE,
                                       "Immunoproteasome activity")
p_int_llvy  <- create_interaction_plot(fac_llvy, on_log_scale = FALSE,
                                       "Standard proteasome activity")
p_int_ratio <- create_interaction_plot(fac_ratio, on_log_scale = TRUE,
                                       "Immunoproteasome shift (ratio)")

# ==============================================================================
# 12. RENDER PAGE 2: Kinetic (top) + Interaction (bottom)
# ==============================================================================
# Helper: prepend a bold panel label (a, b, c, ...) to a plot in top-left
add_panel_label <- function(p, label) {
  p + labs(tag = label) +
    theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0.02, 0.985)
    )
}

# Apply panel labels (a-e). Top row: kinetics (a, b). Bottom row: interaction (c, d, e).
p_kin_anw_lbl   <- add_panel_label(p_kin_anw,   "a")
p_kin_llvy_lbl  <- add_panel_label(p_kin_llvy,  "b")
p_int_anw_lbl   <- add_panel_label(p_int_anw,   "c")
p_int_llvy_lbl  <- add_panel_label(p_int_llvy,  "d")
p_int_ratio_lbl <- add_panel_label(p_int_ratio, "e")

render_figure_p2 <- function(out_path, dev_fun) {
  # A4 landscape (11.69 x 8.27 inch) for publication standard.
  dev_fun(out_path, width = 11.69, height = 8.27)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 5, ncol = 17,
    heights = c(0.7,   # title box
                3.3,   # kinetic plots
                0.35,  # spacer between rows
                3.3,   # interaction plots
                0.7),  # shared legend
    widths = rep(1, 17)
  )))
  
  # Title box (full width)
  print(draw_box("Kinetics & Factorial (LPS \u00D7 Virus) Decomposition",
                 fill_c = "#FBEFE3"),
        vp = viewport(layout.pos.row = 1, layout.pos.col = 1:17))
  
  # ---- Top row: 2 kinetic plots (each = 8 cols, gap = 1 col) ----
  print(p_kin_anw_lbl,
        vp = viewport(layout.pos.row = 2, layout.pos.col = 1:8))
  print(p_kin_llvy_lbl,
        vp = viewport(layout.pos.row = 2, layout.pos.col = 10:17))
  
  # ---- Bottom row: 3 interaction plots (each = 5 cols, gaps = 1 col each) ----
  print(p_int_anw_lbl,
        vp = viewport(layout.pos.row = 4, layout.pos.col = 1:5))
  print(p_int_llvy_lbl,
        vp = viewport(layout.pos.row = 4, layout.pos.col = 7:11))
  print(p_int_ratio_lbl,
        vp = viewport(layout.pos.row = 4, layout.pos.col = 13:17))
  
  # Shared 6-group legend
  pushViewport(viewport(layout.pos.row = 5, layout.pos.col = 1:17))
  grid.draw(shared_legend)
  upViewport()
  upViewport()
  dev.off()
}

pdf_path_p2 <- file.path(out_dir_pdf,
                         "Proteasome_Kinetics_Factorial.pdf")
png_path_p2 <- file.path(out_dir_png,
                         "Proteasome_Kinetics_Factorial.png")

# IMPORTANT: use type="cairo" for PNG so Unicode chars (η²p, ×, ±) render correctly.
# The default bitmap font in non-cairo PNG can't display these glyphs.
render_figure_p2(pdf_path_p2, function(p, ...) cairo_pdf(p, ..., onefile = TRUE))
render_figure_p2(png_path_p2, function(p, ...) png(p, ..., units = "in", res = 300, type = "cairo"))

message("Figure 2 (kinetic + interaction) saved.")
message("All outputs:")
message("  Figure 1 PDF: ", pdf_path)
message("  Figure 1 PNG: ", png_path)
message("  Figure 2 PDF: ", pdf_path_p2)
message("  Figure 2 PNG: ", png_path_p2)
message("  Stats log:    ", stats_path)
message("Done.")