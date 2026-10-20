# ==============================================================================
# 02_Rate-Estimation-QC.R
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================
library(tidyverse)
library(readxl)
library(lubridate)
library(hms)
library(openxlsx)
library(ggpubr)
library(grid)

PROT_DIR <- "Proteasome"   # proteasome data folder (contains the plate-reader workbook)
data_file <- file.path(PROT_DIR, "260515_Virus-HMC3_Proteasome-activity.xlsx")
out_dir <- file.path(PROT_DIR, "results", "rate_compare")
data_out_dir <- file.path(PROT_DIR, "processed")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(data_out_dir)) dir.create(data_out_dir, recursive = TRUE)

group_order <- c(
  "Naive", "A/PR/8/34", "A/NWS/33",
  "LPS", "LPS+A/PR/8/34", "LPS+A/NWS/33"
)

# Tuning parameters for linear window detection
R2_THRESHOLD <- 0.98       # Median R^2 across Total wells must exceed this
MIN_WINDOW_MIN <- 30       # Window must be at least this many minutes long
TIME_FULL_START <- 0
TIME_FULL_END   <- 240

# ==============================================================================
# 2. READ DATA
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

merged_data <- inner_join(tidy_data, master_map, by = "Well") %>%
  filter(!is.na(Substrate))

# Get available timepoints (assume same grid across wells)
all_times <- sort(unique(merged_data$Time_Min))
cat("Total timepoints:", length(all_times),
    "  Range:", min(all_times), "to", max(all_times), "min\n")

# ==============================================================================
# 3. RATE CALCULATION FUNCTIONS
# ==============================================================================
calc_delta <- function(df, t_start, t_end) {
  df_win <- df %>% filter(Time_Min >= t_start & Time_Min <= t_end) %>% arrange(Time_Min)
  if (nrow(df_win) < 2) return(NA_real_)
  dr <- df_win$RFU[nrow(df_win)] - df_win$RFU[1]
  dt <- df_win$Time_Min[nrow(df_win)] - df_win$Time_Min[1]
  if (dt == 0) return(NA_real_)
  dr / dt
}

calc_lm <- function(df, t_start, t_end) {
  # Returns list(slope, r2, n)
  df_win <- df %>% filter(Time_Min >= t_start & Time_Min <= t_end) %>% arrange(Time_Min)
  if (nrow(df_win) < 3 || sd(df_win$RFU, na.rm = TRUE) == 0) {
    return(list(slope = NA_real_, r2 = NA_real_, n = nrow(df_win)))
  }
  fit <- tryCatch(lm(RFU ~ Time_Min, data = df_win), error = function(e) NULL)
  if (is.null(fit)) return(list(slope = NA_real_, r2 = NA_real_, n = nrow(df_win)))
  s <- summary(fit)
  list(slope = unname(coef(fit)["Time_Min"]),
       r2 = s$r.squared,
       n = nrow(df_win))
}

# ==============================================================================
# 4. LINEAR WINDOW PER SUBSTRATE
# ==============================================================================

detect_linear_window <- function(df_substrate,
                                 r2_thresh = R2_THRESHOLD,
                                 min_len = MIN_WINDOW_MIN,
                                 step_min = 10) {
  # df_substrate: long-format with columns Well, Time_Min, RFU (Total only)
  wells <- unique(df_substrate$Well)
  t_grid <- seq(TIME_FULL_START, TIME_FULL_END, by = step_min)
  
  # Build candidate windows (t_start < t_end, length >= min_len)
  candidates <- expand_grid(t_start = t_grid, t_end = t_grid) %>%
    filter(t_end - t_start >= min_len)
  
  # Eval each candidate
  res_list <- lapply(seq_len(nrow(candidates)), function(i) {
    ts <- candidates$t_start[i]; te <- candidates$t_end[i]
    r2s <- vapply(wells, function(w) {
      sub <- df_substrate %>% filter(Well == w)
      calc_lm(sub, ts, te)$r2
    }, numeric(1))
    data.frame(
      t_start = ts, t_end = te,
      length_min = te - ts,
      median_r2 = median(r2s, na.rm = TRUE),
      min_r2 = min(r2s, na.rm = TRUE),
      mean_r2 = mean(r2s, na.rm = TRUE),
      n_wells = sum(!is.na(r2s))
    )
  })
  eval_tab <- bind_rows(res_list)
  
  # Pick best
  qualifying <- eval_tab %>% filter(median_r2 > r2_thresh)
  if (nrow(qualifying) == 0) {
    # Fallback: relax threshold and report
    warning(sprintf("No window meets R2 > %.2f. Falling back to highest median R2.",
                    r2_thresh))
    best <- eval_tab %>% arrange(desc(median_r2), desc(length_min)) %>% slice(1)
  } else {
    best <- qualifying %>% arrange(desc(length_min), t_start) %>% slice(1)
  }
  list(best = best, all = eval_tab)
}

# Run detection per substrate (Total-condition wells only — these are the
# biologically meaningful wells we want to characterize)
substrates <- unique(merged_data$Substrate)
windows <- list()
window_search_logs <- list()
for (sub in substrates) {
  cat("\nDetecting linear window for substrate:", sub, "\n")
  df_sub <- merged_data %>% filter(Substrate == sub, Condition == "Total")
  det <- detect_linear_window(df_sub)
  windows[[sub]] <- det$best
  window_search_logs[[sub]] <- det$all %>% mutate(Substrate = sub)
  cat(sprintf("  -> Window: [%.0f, %.0f] min (length=%.0f, median R2=%.4f, min R2=%.4f)\n",
              det$best$t_start, det$best$t_end, det$best$length_min,
              det$best$median_r2, det$best$min_r2))
}

window_summary <- bind_rows(lapply(names(windows), function(s) {
  cbind(data.frame(Substrate = s), windows[[s]])
}))

# ==============================================================================
# 5. COMPUTE 3 RATE METHODS PER WELL
# ==============================================================================
compute_all_methods <- function(df, sub_type) {
  win <- windows[[sub_type]]
  t_lin_start <- win$t_start
  t_lin_end   <- win$t_end
  
  # M1. delta full
  m1_delta_full <- calc_delta(df, TIME_FULL_START, TIME_FULL_END)
  # M2. lm full
  m2 <- calc_lm(df, TIME_FULL_START, TIME_FULL_END)
  # M3. lm linear window
  m3 <- calc_lm(df, t_lin_start, t_lin_end)
  
  data.frame(
    rate_delta_full = m1_delta_full,
    rate_lm_full = m2$slope,  r2_lm_full = m2$r2,  n_lm_full = m2$n,
    rate_lm_linwin = m3$slope, r2_lm_linwin = m3$r2, n_lm_linwin = m3$n,
    linwin_start = t_lin_start, linwin_end = t_lin_end
  )
}

per_well_rates <- merged_data %>%
  group_by(Group, Substrate, Condition, Well) %>%
  nest() %>%
  mutate(rates = map2(data, Substrate, compute_all_methods)) %>%
  select(-data) %>%
  unnest(rates) %>%
  ungroup()

# ==============================================================================
# 6. APPLY BACKGROUND SUBTRACTION & COMPUTE SPECIFIC ACTIVITY FOR EACH METHOD
# ==============================================================================
methods <- c("delta_full", "lm_full", "lm_linwin")

apply_bg_sub <- function(df_method, rate_col) {
  bg <- df_method %>%
    filter(Condition == "Background") %>%
    group_by(Group, Substrate) %>%
    summarise(Mean_BG = mean(.data[[rate_col]], na.rm = TRUE), .groups = "drop")
  
  df_method %>%
    filter(Condition == "Total") %>%
    left_join(bg, by = c("Group", "Substrate")) %>%
    mutate(Mean_BG = replace_na(Mean_BG, 0),
           Specific_Activity = .data[[rate_col]] - Mean_BG) %>%
    select(Group, Substrate, Well, Rate = !!rate_col, Mean_BG, Specific_Activity) %>%
    mutate(Method = !!rate_col)
}

specific_activity_all <- bind_rows(
  apply_bg_sub(per_well_rates, "rate_delta_full"),
  apply_bg_sub(per_well_rates, "rate_lm_full"),
  apply_bg_sub(per_well_rates, "rate_lm_linwin")
) %>%
  mutate(
    Method = dplyr::case_match(Method,
                               "rate_delta_full" ~ "M1_delta_full",
                               "rate_lm_full"    ~ "M2_lm_full",
                               "rate_lm_linwin"  ~ "M3_lm_linwin",
                               .default = Method
    ),
    Group = factor(Group, levels = group_order)
  )

# ==============================================================================
# 7. METHOD COMPARISON: CORRELATION & EFFECT ON GROUP MEANS
# ==============================================================================
cor_summary <- specific_activity_all %>%
  select(Substrate, Group, Well, Method, Specific_Activity) %>%
  pivot_wider(names_from = Method, values_from = Specific_Activity) %>%
  group_by(Substrate) %>%
  summarize(
    cor_M1_M2 = cor(M1_delta_full, M2_lm_full, use = "pairwise.complete.obs"),
    cor_M1_M3 = cor(M1_delta_full, M3_lm_linwin, use = "pairwise.complete.obs"),
    cor_M2_M3 = cor(M2_lm_full, M3_lm_linwin, use = "pairwise.complete.obs"),
    .groups = "drop"
  )

group_summary <- specific_activity_all %>%
  group_by(Substrate, Method, Group) %>%
  summarize(
    N = sum(!is.na(Specific_Activity)),
    Mean = mean(Specific_Activity, na.rm = TRUE),
    SD = sd(Specific_Activity, na.rm = TRUE),
    SEM = SD / sqrt(N),
    CV = SD / abs(Mean),
    .groups = "drop"
  )

cv_summary <- group_summary %>%
  group_by(Substrate, Method) %>%
  summarize(
    Mean_CV = mean(CV, na.rm = TRUE),
    Median_CV = median(CV, na.rm = TRUE),
    .groups = "drop"
  )

signal_to_noise <- group_summary %>%
  group_by(Substrate, Method) %>%
  summarize(
    Range_means = max(Mean) - min(Mean),
    Pooled_SD = sqrt(mean(SD^2, na.rm = TRUE)),
    Signal_Noise = Range_means / Pooled_SD,
    .groups = "drop"
  )

# 7e. R^2 quality (only for lm methods)
r2_quality <- per_well_rates %>%
  filter(Condition == "Total") %>%
  group_by(Substrate) %>%
  summarize(
    Median_R2_lm_full = median(r2_lm_full, na.rm = TRUE),
    Min_R2_lm_full = min(r2_lm_full, na.rm = TRUE),
    Median_R2_lm_linwin = median(r2_lm_linwin, na.rm = TRUE),
    Min_R2_lm_linwin = min(r2_lm_linwin, na.rm = TRUE),
    .groups = "drop"
  )

# ==============================================================================
# 8. EXPORT EXCEL
# ==============================================================================
wb <- createWorkbook()

addWorksheet(wb, "0_Methods_Notes")
writeData(wb, "0_Methods_Notes", data.frame(
  Section = c("M1_delta_full", "M2_lm_full", "M3_lm_linwin",
              "Window detection", "Decision rule"),
  Description = c(
    "Simple (RFU_end - RFU_start) / total_time over 0-240 min. Original method.",
    "Slope from lm(RFU ~ Time_Min) over full 0-240 min window. Smooths noise but still affected by lag phase and saturation.",
    "Slope from lm() over auto-detected LINEAR window per substrate. Window picked s.t. median R^2 across all Total wells > 0.98, length >= 30 min, longest. Most rigorous: avoids lag and saturation artefacts.",
    sprintf("R^2 threshold = %.2f; min window length = %d min. Window picked PER SUBSTRATE using Total wells only.", R2_THRESHOLD, MIN_WINDOW_MIN),
    "Choose method with: (a) highest Signal_Noise in sheet 5, (b) lowest within-group CV in sheet 4, (c) high R^2 in sheet 6. M3 is expected to win for assays with lag or saturation."
  )
))

addWorksheet(wb, "1_Linear_Windows")
writeData(wb, "1_Linear_Windows", window_summary)

addWorksheet(wb, "2_PerWell_Rates")
writeData(wb, "2_PerWell_Rates", per_well_rates)

addWorksheet(wb, "3_Method_Correlations")
writeData(wb, "3_Method_Correlations", cor_summary)

addWorksheet(wb, "4_Group_Means_by_Method")
writeData(wb, "4_Group_Means_by_Method", group_summary)

addWorksheet(wb, "5_CV_and_SignalNoise")
writeData(wb, "5_CV_and_SignalNoise",
          left_join(cv_summary, signal_to_noise, by = c("Substrate", "Method")))

addWorksheet(wb, "6_R2_Quality")
writeData(wb, "6_R2_Quality", r2_quality)

addWorksheet(wb, "7_SpecificActivity_Long")
writeData(wb, "7_SpecificActivity_Long", specific_activity_all)

# Reorder sheets to put Methods_Notes first
worksheetOrder(wb) <- order(names(wb))

excel_path <- file.path(out_dir, "Proteasome_RateMethod_Compare.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)
message("Comparison Excel saved: ", excel_path)

# ==============================================================================
# 9. VISUAL COMPARISON FIGURE
# ==============================================================================
method_colors <- c("M1_delta_full" = "#9AA0A6",
                   "M2_lm_full" = "#4285F4",
                   "M3_lm_linwin" = "#C0392B")

group_colors <- c(
  "Naive"          = "#D5D8DC",
  "A/PR/8/34"      = "#A8C7F0",
  "A/NWS/33"       = "#F2B4B4",
  "LPS"            = "#6E7376",
  "LPS+A/PR/8/34"  = "#4285F4",
  "LPS+A/NWS/33"   = "#C0392B"
)

make_kinetic_panel <- function(sub_name) {
  df <- merged_data %>%
    filter(Substrate == sub_name, Condition == "Total") %>%
    mutate(Group = factor(Group, levels = group_order))
  win <- windows[[sub_name]]
  ggplot(df, aes(x = Time_Min, y = RFU,
                 color = Group, group = Well)) +
    annotate("rect", xmin = win$t_start, xmax = win$t_end,
             ymin = -Inf, ymax = Inf, fill = "#FFF3CD", alpha = 0.4) +
    annotate("segment",
             x = win$t_start, xend = win$t_start, y = -Inf, yend = Inf,
             color = "#856404", linetype = "dashed", linewidth = 0.5) +
    annotate("segment",
             x = win$t_end, xend = win$t_end, y = -Inf, yend = Inf,
             color = "#856404", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.5, alpha = 0.7) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(
      title = paste0(sub_name, " — Linear window: [",
                     win$t_start, ", ", win$t_end, "] min  ",
                     sprintf("(median R\u00B2=%.3f)", win$median_r2)),
      x = "Time (min)", y = "RFU"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold"))
}

make_method_bar <- function(sub_name) {
  df <- group_summary %>%
    filter(Substrate == sub_name) %>%
    mutate(Group = factor(Group, levels = group_order))
  ggplot(df, aes(x = Group, y = Mean, fill = Method)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(Mean - SEM, 0), ymax = Mean + SEM),
                  position = position_dodge(width = 0.8),
                  width = 0.2, linewidth = 0.4) +
    scale_fill_manual(values = method_colors,
                      labels = c("M1: delta full",
                                 "M2: lm full",
                                 "M3: lm linear-window")) +
    labs(title = paste0(sub_name, " — Group means by method"),
         x = NULL, y = "Specific Activity (RFU/min)", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
          legend.position = "top",
          plot.title = element_text(size = 10, face = "bold"))
}

sn_plot <- signal_to_noise %>%
  ggplot(aes(x = Method, y = Signal_Noise, fill = Method)) +
  geom_col(width = 0.6, alpha = 0.85) +
  facet_wrap(~ Substrate, scales = "free_y") +
  scale_fill_manual(values = method_colors) +
  labs(title = "Signal-to-Noise (range of group means / pooled SD)",
       subtitle = "Higher is better — more power for pairwise tests",
       x = NULL, y = "Signal / Noise") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        plot.title = element_text(size = 11, face = "bold"))

fig <- ggpubr::ggarrange(
  make_kinetic_panel(substrates[1]),
  make_kinetic_panel(substrates[2]),
  make_method_bar(substrates[1]),
  make_method_bar(substrates[2]),
  sn_plot, NULL,
  ncol = 2, nrow = 3,
  heights = c(1, 1, 1),
  labels = c("a", "b", "c", "d", "e", ""),
  font.label = list(size = 14, face = "bold")
)

png_path <- file.path(out_dir, "Proteasome_RateMethod_Compare.png")
pdf_path <- file.path(out_dir, "Proteasome_RateMethod_Compare.pdf")

ggsave(png_path, fig, width = 12, height = 12, dpi = 300, type = "cairo")
ggsave(pdf_path, fig, width = 12, height = 12, device = cairo_pdf)

message("Comparison figure saved:")
message("  PNG: ", png_path)
message("  PDF: ", pdf_path)

# ==============================================================================
# 10. PRINT VERDICT TO CONSOLE
# ==============================================================================
cat("\n", strrep("=", 78), "\n", sep = "")
cat("RATE METHOD COMPARISON — SUMMARY\n")
cat(strrep("=", 78), "\n\n", sep = "")

cat("Linear windows detected:\n")
for (s in substrates) {
  w <- windows[[s]]
  cat(sprintf("  %s: [%g, %g] min, length=%g, median R2=%.4f\n",
              s, w$t_start, w$t_end, w$length_min, w$median_r2))
}

cat("\nMethod correlations (per substrate):\n")
print(cor_summary, row.names = FALSE)

cat("\nSignal-to-Noise by method:\n")
print(signal_to_noise %>%
        arrange(Substrate, desc(Signal_Noise)),
      row.names = FALSE)

cat("\nR2 quality (lm methods only):\n")
print(r2_quality, row.names = FALSE)

cat("\n", strrep("-", 78), "\n", sep = "")
cat("HEURISTIC: Pick the method with HIGHEST Signal_Noise per substrate.\n")
cat("If methods agree (cor > 0.99), simpler is fine. If M3 wins notably,\n")
cat("update RATE_METHOD <- \"regression\" in the main analysis script and\n")
cat("plug in the windows reported above.\n")
cat(strrep("=", 78), "\n", sep = "")