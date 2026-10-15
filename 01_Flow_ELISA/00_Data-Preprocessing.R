
################################################################################
# 00_Data-Preprocessing.R — Data Ingestion & Schema Setup
# Load data from flow cytometry (as "% of total live cells") and ELISA, then transform into a tidy, long-format data frame.
# Report sample sizes (n) for each experimental group, identify missing values, determine the range of values ​​for each biomarker, and provide a list of all biomarkers.
# V1 = A/PR/8/34, V2 = A/NWS/33
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

set.seed(2026)

# --Path setup -----
base_dir <- "Flow-ELISA" 
out_dir  <- file.path(base_dir, "R_analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load % alive cells ------
facs_raw <- read_excel(file.path(data_dir, "FACS_ForR.xlsx"), sheet = 1)
cat("FACS raw dim:", nrow(facs_raw), "x", ncol(facs_raw), "\n")
cat("FACS columns:", paste(names(facs_raw), collapse = ", "), "\n\n")


flow_markers <- c("CD40", "CD68", "CD86", "CD206", "HLA_I", "HLA_II", "NP")
flow_wide <- facs_raw %>%
  select(sample_id, timepoint, group, rep, all_of(flow_markers))

cat("Flow wide dim (% total live):", nrow(flow_wide), "x", ncol(flow_wide), "\n")

# --- 2. Load ELISA ------
elisa_raw <- read_excel(file.path(data_dir, "ELISA_ForR.xlsx"), sheet = 1)
cat("ELISA raw dim:", nrow(elisa_raw), "x", ncol(elisa_raw), "\n")
cat("ELISA columns:", paste(names(elisa_raw), collapse = ", "), "\n\n")

elisa_markers <- c("CCL2", "IL_1b", "IL_6", "IL_10", "TNFa", "IL_12p70")
elisa_wide <- elisa_raw %>%
  select(sample_id, group, timepoint, rep, all_of(elisa_markers))

# --- 3. Factorial design columns ------
add_design <- function(df) {
  df %>%
    mutate(
      Timepoint = as.integer(str_extract(timepoint, "\\d+")),
      LPS = case_when(
        group %in% c("LPS", "LPS_PR8", "LPS_S33") ~ "Yes",
        TRUE ~ "No"
      ),
      Virus = case_when(
        group %in% c("C_PR8", "LPS_PR8") ~ "V1",
        group %in% c("C_S33", "LPS_S33") ~ "V2",
        group == "PolyIC" ~ "PolyIC",
        TRUE ~ "None"
      ),
      Bio_Rep = as.integer(rep),
      Group = group  # keep original group label
    )
}

flow_wide <- add_design(flow_wide)
elisa_wide <- add_design(elisa_wide)

# --- 4. Pivot to tidy long ------
# Flow
flow_long <- flow_wide %>%
  pivot_longer(
    cols = all_of(flow_markers),
    names_to = "Marker",
    values_to = "Value"
  ) %>%
  mutate(
    Platform = "flow",
    Unit = "pct_total_live",
    # Standardize marker names for display
    Marker = recode(Marker,
      "HLA_I" = "HLA-I",
      "HLA_II" = "HLA-II"
    )
  ) %>%
  select(Group, Timepoint, LPS, Virus, Bio_Rep, Platform, Marker, Value, Unit)

# ELISA
elisa_long <- elisa_wide %>%
  pivot_longer(
    cols = all_of(elisa_markers),
    names_to = "Marker",
    values_to = "Value"
  ) %>%
  mutate(
    Platform = "elisa",
    Unit = "pg_per_mL",
    # Standardize marker names for display
    Marker = recode(Marker,
      "IL_1b"   = "IL-1b",
      "IL_6"    = "IL-6",
      "IL_10"   = "IL-10",
      "TNFa"    = "TNF-a",
      "IL_12p70" = "IL-12p70"
    )
  ) %>%
  select(Group, Timepoint, LPS, Virus, Bio_Rep, Platform, Marker, Value, Unit)

# Combine
df_all <- bind_rows(flow_long, elisa_long)
cat("Combined tidy df dim:", nrow(df_all), "x", ncol(df_all), "\n\n")

# --- 5. Separate Poly(I:C) ----
df_factorial <- df_all %>% filter(Virus != "PolyIC")
df_polyic    <- df_all %>% filter(Virus == "PolyIC")

cat("Factorial rows:", nrow(df_factorial), "\n")
cat("Poly(I:C) rows:", nrow(df_polyic), "\n\n")

# --6. Checking data -----
# 6a. Cell counts (n per group × tp × marker)
cat("── n per cell (factorial, expect 5?) ──────────────\n")
cell_n <- df_factorial %>%
  group_by(Group, Timepoint, Platform, Marker) %>%
  summarise(n = sum(!is.na(Value)), n_miss = sum(is.na(Value)), .groups = "drop")

# Summary
cat("Min n:", min(cell_n$n), " Max n:", max(cell_n$n), "\n")
if (any(cell_n$n != 5)) {
  cat("⚠️  Cells with n ≠ 5:\n")
  print(cell_n %>% filter(n != 5))
} else {
  cat("✓  All cells have n = 5\n")
}
cat("\n")

# 6b. Missing values
cat("── Missing values ────────────────────────────────\n")
miss_total <- sum(is.na(df_all$Value))
cat("Total NAs in combined df:", miss_total, "\n")
if (miss_total > 0) {
  miss_detail <- df_all %>%
    filter(is.na(Value)) %>%
    group_by(Platform, Marker, Group, Timepoint) %>%
    summarise(n_miss = n(), .groups = "drop")
  cat("Missing detail:\n")
  print(miss_detail)
}
cat("\n")

# 6c. Poly(I:C) cell counts
cat("── Poly(I:C) n per cell ──────────────────────────\n")
poly_n <- df_polyic %>%
  group_by(Timepoint, Marker, Platform) %>%
  summarise(n = sum(!is.na(Value)), .groups = "drop")
cat("Min n:", min(poly_n$n), " Max n:", max(poly_n$n), "\n")
if (any(poly_n$n != 5)) {
  cat("⚠️  Poly(I:C) cells with n ≠ 5:\n")
  print(poly_n %>% filter(n != 5))
} else {
  cat("✓  All Poly(I:C) cells have n = 5\n")
}
cat("\n")

# 6d. Range per marker
cat("── Range per marker ──────────────────────────────\n")
marker_range <- df_all %>%
  group_by(Platform, Marker, Unit) %>%
  summarise(
    min_val  = min(Value, na.rm = TRUE),
    max_val  = max(Value, na.rm = TRUE),
    mean_val = mean(Value, na.rm = TRUE),
    sd_val   = sd(Value, na.rm = TRUE),
    n_obs    = sum(!is.na(Value)),
    .groups = "drop"
  ) %>%
  arrange(Platform, Marker)
print(marker_range, n = Inf)
cat("\n")

# 6e. Marker list
cat("── Marker list ───────────────────────────────────\n")
cat("Flow markers (% total live):", paste(sort(unique(flow_long$Marker)), collapse = ", "), "\n")
cat("ELISA markers (pg/mL):", paste(sort(unique(elisa_long$Marker)), collapse = ", "), "\n")
cat("\n")

# 6f. Design structure
cat("── Design structure (factorial) ───────────────────\n")
cat("Groups:", paste(sort(unique(df_factorial$Group)), collapse = ", "), "\n")
cat("Timepoints:", paste(sort(unique(df_factorial$Timepoint)), collapse = ", "), "\n")
cat("LPS levels:", paste(sort(unique(df_factorial$LPS)), collapse = ", "), "\n")
cat("Virus levels:", paste(sort(unique(df_factorial$Virus)), collapse = ", "), "\n")
cat("Bio_Rep range:", range(df_factorial$Bio_Rep), "\n\n")

# Cross-tab: Group × Timepoint (should be 5 per cell for each marker)
cat("── Cross-tab: unique Bio_Reps per Group×Tp ───────\n")
xtab <- df_factorial %>%
  filter(Marker == "CD86") %>%
  group_by(Group, Timepoint) %>%
  summarise(n = n_distinct(Bio_Rep), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = n)
print(xtab)
cat("\n")

# ── 7. Save tidy data ─────────────────────────────────────────────────────────
write.csv(df_all, file.path(out_dir, "tidy_all.csv"), row.names = FALSE)
write.csv(df_factorial, file.path(out_dir, "tidy_factorial.csv"), row.names = FALSE)
write.csv(df_polyic, file.path(out_dir, "tidy_polyic.csv"), row.names = FALSE)

cat("── Files saved ───────────────────────────────────\n")
cat("  tidy_all.csv       :", nrow(df_all), "rows\n")
cat("  tidy_factorial.csv :", nrow(df_factorial), "rows\n")
cat("  tidy_polyic.csv    :", nrow(df_polyic), "rows\n")
cat("\n═══ 00_Data-Preprocessing COMPLETE ═══\n")
