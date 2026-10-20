# =============================================================================
# 02_Preprocessing-Pipeline.R  —  RT-qPCR preprocessing (amplification → fold change)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(plotly)
  library(htmlwidgets)
})

# =============================================================================
# SS 1  GLOBAL CONFIGURATION
# =============================================================================

PROJ_DIR <- "./2606_PCR_VirusHMC3"
DATA_DIR <- file.path(PROJ_DIR, "DATA", "Excel")
OUT_DIR  <- file.path(PROJ_DIR, "analysis_master")

# Dual reference genes (geNORM validated)
REFERENCE_GENES <- c("GAPDH", "RPS11")

# Genes purged from the main 8-target pipeline.
# SNCA is kept here — it is processed separately in Stage 6.
PURGED_GENES <- c("SNCA", "SNCA-140", "RPL32", "RPS")

# Standard-curve efficiencies from 01_Efficiency-Calibration.R
# CSV columns: Gene, Efficiency, R2
CALIB_EFF_PATH <- "./2606_PCR_VirusHMC3/analysis/final-figures/final_efficiencies.csv"

eff_table <- tryCatch(
  read_csv(CALIB_EFF_PATH, show_col_types = FALSE) %>%
    mutate(
      Gene       = trimws(as.character(Gene)),
      Efficiency = 1 + eff_canon / 100,   # convert % -> Pfaffl E (e.g. 90.58% -> 1.906)
      R2         = r2_canon
    ) %>%
    select(Gene, Efficiency, R2),
  error = function(e) {
    stop("FATAL: Cannot read efficiency file at:\n  ", CALIB_EFF_PATH,
         "\n  Original error: ", conditionMessage(e))
  }
)
message("Loaded efficiencies for ", nrow(eff_table), " genes from: ", basename(CALIB_EFF_PATH))

# Instrument export aliases -> canonical gene names
GENE_ALIASES <- c("SNCA-140" = "SNCA", "OSA1" = "OAS1",
                   "SPL32"    = "RPL32", "SPS11" = "RPS11")

# Split-control calibrators
CTRL_4H <- "G1"
CTRL_8H <- "G8"

# ── Automated Melt QC thresholds ──────────────────────────────────────────────
TM_DEVIATION_TOL     <- 1.0    # degrees C
MIN_PEAK_HEIGHT      <- 0.02   # -dF/dT; below this = no real product
SECONDARY_PEAK_RATIO <- 0.50   # secondary peak >= 50 pct of primary = distinct 2nd band

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Output directory : ", OUT_DIR)

# =============================================================================
# SS 2  HELPER FUNCTIONS
# =============================================================================

safe_read_sheet <- function(path, sheet, skip = 0) {
  tryCatch(
    read_excel(path, sheet = sheet, skip = skip,
               col_names = TRUE, .name_repair = "unique"),
    error = function(e) { warning("Sheet '", sheet, "' not found in ", basename(path)); NULL }
  )
}

apply_aliases <- function(x) {
  x <- trimws(as.character(x))
  ifelse(x %in% names(GENE_ALIASES), GENE_ALIASES[x], x)
}

is_purged <- function(x) trimws(x) %in% PURGED_GENES

get_eff <- function(gene) {
  g   <- trimws(as.character(gene))
  hit <- eff_table$Efficiency[eff_table$Gene == g]
  if (length(hit) == 1L && !is.na(hit)) return(hit)
  message("  [get_eff] '", g, "' not in calibration CSV — falling back to E = 2.0")
  2.0
}

# Decode "G11" -> Group_ID="G1", Bio_Rep="Rep1", Timepoint="4h"
decode_sample_df <- function(snames) {
  purrr::map_dfr(snames, function(s) {
    s <- trimws(as.character(s))
    if (!grepl("^G\\d{2,3}$", s))
      return(tibble(Group_ID = s, Group_Num = NA_integer_,
                    Bio_Rep  = NA_character_, Timepoint = NA_character_))
    digits    <- sub("^G", "", s)
    rep_digit <- as.integer(substr(digits, nchar(digits), nchar(digits)))
    grp_num   <- as.integer(substr(digits, 1L, nchar(digits) - 1L))
    tibble(
      Group_ID  = paste0("G", grp_num),
      Group_Num = grp_num,
      Bio_Rep   = paste0("Rep", rep_digit),
      Timepoint = if (!is.na(grp_num) && grp_num <= 7L) "4h" else "8h"
    )
  })
}

# Second Derivative Maximum Ct — threshold-free amplification onset
sdm_ct <- function(cycles, rn, fine_step = 0.05) {
  if (length(cycles) < 8L || all(is.na(rn))) return(NA_real_)
  rn[is.na(rn)] <- 0
  sp <- tryCatch(smooth.spline(cycles, rn, spar = 0.5), error = function(e) NULL)
  if (is.null(sp)) return(NA_real_)
  fine <- seq(min(cycles), max(cycles), by = fine_step)
  ys   <- predict(sp, fine)$y
  d2   <- diff(diff(ys)) / fine_step ^ 2
  d2x  <- fine[seq(2L, length(fine) - 1L)]
  pk   <- which.max(d2)
  if (!length(pk)) return(NA_real_)
  if (pk > 1L && pk < length(d2)) {
    y0 <- d2[pk - 1L]; y1 <- d2[pk]; y2 <- d2[pk + 1L]
    denom <- 2 * y1 - y0 - y2
    if (abs(denom) > .Machine$double.eps)
      return(d2x[pk] + fine_step * 0.5 * (y0 - y2) / denom)
  }
  d2x[pk]
}

# Pivot ABI wide melt sheet to long format
pivot_melt_sheet <- function(df, value_name) {
  df %>%
    select(`Well Location`, Target, starts_with("Reading")) %>%
    rename(Well = `Well Location`, Gene = Target) %>%
    mutate(Well = as.character(trimws(Well)), Gene = apply_aliases(Gene)) %>%
    filter(!is.na(Well), Well != "", !is_purged(Gene)) %>%
    pivot_longer(starts_with("Reading"), names_to = "Reading_Index", values_to = value_name) %>%
    mutate(
      Reading_Index    = as.integer(str_extract(Reading_Index, "[0-9]+")),
      !!value_name := suppressWarnings(as.numeric(.data[[value_name]]))
    ) %>%
    filter(!is.na(.data[[value_name]]))
}

# ── Multi-peak detector (see Melt QC threshold documentation in SS 1) ─────────
# Returns TRUE when a secondary peak is "real" (tall enough to be a distinct band).
detect_multi_peak <- function(temps, derivs,
                               min_height   = MIN_PEAK_HEIGHT,
                               valley_ratio = SECONDARY_PEAK_RATIO) {
  if (length(derivs) < 5L || all(is.na(derivs))) return(FALSE)
  # Light smoothing before peak counting to suppress single-point spikes
  smoothed <- tryCatch(
    predict(smooth.spline(temps, derivs, spar = 0.4), temps)$y,
    error = function(e) derivs
  )
  n <- length(smoothed)
  is_local_max <- c(FALSE,
                    smoothed[seq(2L, n - 1L)] > smoothed[seq(1L, n - 2L)] &
                    smoothed[seq(2L, n - 1L)] > smoothed[seq(3L, n)],
                    FALSE)
  cand_idx <- which(is_local_max & smoothed >= min_height)
  if (length(cand_idx) <= 1L) return(FALSE)
  heights   <- smoothed[cand_idx]
  primary_h <- max(heights)
  any(sort(heights, decreasing = TRUE)[-1L] >= valley_ratio * primary_h)
}

# ── geNORM dual-reference Pfaffl — Bio_Rep-averaged input ────────────────────
pfaffl_per_rep <- function(ct_df, ctrl_group) {
  ref_ctrl_mean <- ct_df %>%
    filter(`Target Name` %in% REFERENCE_GENES, Group_ID == ctrl_group) %>%
    group_by(`Target Name`) %>%
    summarise(Ct_ref_ctrl = mean(Mean_Ct, na.rm = TRUE), .groups = "drop")

  ref_rq <- ct_df %>%
    filter(`Target Name` %in% REFERENCE_GENES) %>%
    left_join(ref_ctrl_mean, by = "Target Name") %>%
    mutate(E_ref  = map_dbl(`Target Name`, get_eff),
           RQ_ref = E_ref ^ (Ct_ref_ctrl - Mean_Ct))

  nf_df <- ref_rq %>%
    group_by(Bio_Rep, Group_ID) %>%
    summarise(
      NF = exp(mean(log(pmax(RQ_ref, .Machine$double.eps)), na.rm = TRUE)),
      .groups = "drop"
    )

  tgt_ctrl_mean <- ct_df %>%
    filter(!(`Target Name` %in% REFERENCE_GENES), Group_ID == ctrl_group) %>%
    group_by(`Target Name`) %>%
    summarise(Ct_tgt_ctrl = mean(Mean_Ct, na.rm = TRUE), .groups = "drop")

  ct_df %>%
    filter(!(`Target Name` %in% REFERENCE_GENES)) %>%
    left_join(tgt_ctrl_mean, by = "Target Name") %>%
    left_join(nf_df,         by = c("Bio_Rep", "Group_ID")) %>%
    mutate(E_target   = map_dbl(`Target Name`, get_eff),
           FoldChange = (E_target ^ (Ct_tgt_ctrl - Mean_Ct)) / NF) %>%
    select(Bio_Rep, Group_ID, `Target Name`, FoldChange)
}

# ── geNORM dual-reference Pfaffl — well-level input ──────────────────────────

pfaffl_well_level <- function(ct_tbl, ctrl_group) {
  # Step 1: average reference gene Cts within Bio_Rep (collapses tech reps)
  ref_bio_avg <- ct_tbl %>%
    filter(`Target Name` %in% REFERENCE_GENES, !is.na(Ct_SDM)) %>%
    group_by(`Target Name`, Group_ID, Bio_Rep) %>%
    summarise(Mean_Ct_ref = mean(Ct_SDM, na.rm = TRUE), .groups = "drop")

  # Step 2: calibrator group-mean ref Ct (one scalar per reference gene)
  ref_ctrl_mean <- ref_bio_avg %>%
    filter(Group_ID == ctrl_group) %>%
    group_by(`Target Name`) %>%
    summarise(Ct_ref_ctrl = mean(Mean_Ct_ref, na.rm = TRUE), .groups = "drop")

  # Step 3: geNORM NF per Bio_Rep — geometric mean across all reference genes
  nf_df <- ref_bio_avg %>%
    left_join(ref_ctrl_mean, by = "Target Name") %>%
    mutate(E_ref  = map_dbl(`Target Name`, get_eff),
           RQ_ref = E_ref ^ (Ct_ref_ctrl - Mean_Ct_ref)) %>%
    group_by(Group_ID, Bio_Rep) %>%
    summarise(NF = exp(mean(log(pmax(RQ_ref, .Machine$double.eps)), na.rm = TRUE)),
              .groups = "drop")

  # Step 4: calibrator group-mean Ct for each target gene (one scalar per gene)
  tgt_ctrl_mean <- ct_tbl %>%
    filter(!(`Target Name` %in% REFERENCE_GENES),
           Group_ID == ctrl_group, !is.na(Ct_SDM)) %>%
    group_by(`Target Name`) %>%
    summarise(Ct_tgt_ctrl = mean(Ct_SDM, na.rm = TRUE), .groups = "drop")

  # Step 5: FC for every individual well using its own Ct_SDM
  ct_tbl %>%
    filter(!(`Target Name` %in% REFERENCE_GENES), !is.na(Ct_SDM)) %>%
    left_join(tgt_ctrl_mean, by = "Target Name") %>%
    left_join(nf_df,         by = c("Group_ID", "Bio_Rep")) %>%
    mutate(E_target   = map_dbl(`Target Name`, get_eff),
           FoldChange = (E_target ^ (Ct_tgt_ctrl - Ct_SDM)) / NF) %>%
    select(Well, `Target Name`, Group_ID, Bio_Rep, Timepoint,
           Ct_SDM, Exclude.amplicon, FoldChange)
}

# =============================================================================
# SS 3  FILE INGESTION LOOP
# =============================================================================

data_files <- list.files(DATA_DIR, pattern = "\\.xlsx?$",
                          full.names = TRUE, ignore.case = TRUE)
if (!length(data_files)) stop("No Excel files found in: ", DATA_DIR)

message("\n=== FILE SCAN: ", length(data_files), " file(s) found ===")
walk(data_files, ~ message("  ", basename(.x)))

master_amp  <- tibble()
snca_amp    <- tibble()   # SNCA collected separately (bypasses main purge)
master_melt <- tibble()

for (fpath in data_files) {
  fname    <- basename(fpath)
  is_check <- grepl("Check-sang1", fname, ignore.case = TRUE)
  is_rep5  <- grepl("REP5-remain", fname, ignore.case = TRUE)

  message(
    "\n--- ", fname,
    if (is_check) "  [CHECK-SANG1: amp only]"   else "",
    if (is_rep5)  "  [REP5-remain: Well C4 excluded]" else ""
  )

  ss <- safe_read_sheet(fpath, "Results", skip = 7)
  if (is.null(ss)) { warning("Skipping ", fname); next }

  sample_map <- ss %>%
    select(any_of(c("Well", "Sample Name", "Target Name"))) %>%
    filter(!is.na(Well), !is.na(`Target Name`), !is.na(`Sample Name`)) %>%
    mutate(Well          = as.character(trimws(Well)),
           `Target Name` = apply_aliases(`Target Name`),
           `Sample Name` = trimws(as.character(`Sample Name`))) %>%
    filter(!is_purged(`Target Name`)) %>%
    distinct() %>%
    bind_cols(decode_sample_df(.$`Sample Name`))

  # Well -> group mapping for main pipeline (purged genes excluded)
  well_map <- sample_map %>%
    select(Well, `Sample Name`, Group_ID, Group_Num, Bio_Rep, Timepoint) %>%
    distinct()

  snca_ss <- ss %>%
    select(any_of(c("Well", "Sample Name", "Target Name"))) %>%
    filter(!is.na(Well), !is.na(`Target Name`), !is.na(`Sample Name`)) %>%
    mutate(Well          = as.character(trimws(Well)),
           `Target Name` = apply_aliases(`Target Name`),
           `Sample Name` = trimws(as.character(`Sample Name`))) %>%
    filter(`Target Name` == "SNCA") %>%
    distinct()

  snca_well_map <- if (nrow(snca_ss) > 0) {
    snca_ss %>%
      bind_cols(decode_sample_df(.$`Sample Name`)) %>%
      select(Well, `Sample Name`, Group_ID, Group_Num, Bio_Rep, Timepoint) %>%
      distinct()
  } else {
    tibble(Well = character(), `Sample Name` = character(),
           Group_ID = character(), Group_Num = integer(),
           Bio_Rep  = character(), Timepoint  = character())
  }

  raw_amp <- safe_read_sheet(fpath, "Amplification Data", skip = 7)
  if (!is.null(raw_amp)) {

    # ── Main pipeline amplification chunk ──────────────────────────────────
    amp_chunk <- raw_amp %>%
      select(any_of(c("Well", "Cycle", "Target Name", "Rn", "DeltaRn", "ΔRn"))) %>%
      filter(!is.na(Well), !is.na(Cycle), !is.na(Rn)) %>%
      mutate(Well          = as.character(trimws(Well)),
             Cycle         = suppressWarnings(as.integer(Cycle)),
             Rn            = suppressWarnings(as.numeric(Rn)),
             `Target Name` = apply_aliases(`Target Name`)) %>%
      rename_with(~ "Delta_Rn", any_of(c("DeltaRn", "ΔRn"))) %>%
      { if (!"Delta_Rn" %in% names(.)) mutate(., Delta_Rn = NA_real_) else . } %>%
      mutate(Delta_Rn = suppressWarnings(as.numeric(Delta_Rn))) %>%
      filter(!is.na(Cycle), !is.na(Rn), !is_purged(`Target Name`)) %>%
      left_join(well_map, by = "Well") %>%
      filter(!is.na(Group_ID), !is.na(Bio_Rep)) %>%
      mutate(File = fname)

    if (is_rep5) {
      n_b       <- nrow(amp_chunk)
      amp_chunk <- amp_chunk %>% filter(Well != "C4")
      message("  [REP5-remain] Dropped Well C4 from main amp: ",
              n_b - nrow(amp_chunk), " rows removed")
    }
    if (is_check) amp_chunk <- filter(amp_chunk, Group_ID == "G1", Timepoint == "4h")

    master_amp <- bind_rows(master_amp, amp_chunk)
    message("  amp rows: ", nrow(amp_chunk))

    # ── SNCA exploratory chunk ─────────────────────────
    snca_chunk <- raw_amp %>%
      select(any_of(c("Well", "Cycle", "Target Name", "Rn"))) %>%
      filter(!is.na(Well), !is.na(Cycle), !is.na(Rn)) %>%
      mutate(Well          = as.character(trimws(Well)),
             Cycle         = suppressWarnings(as.integer(Cycle)),
             Rn            = suppressWarnings(as.numeric(Rn)),
             `Target Name` = apply_aliases(trimws(as.character(`Target Name`)))) %>%
      filter(!is.na(Cycle), !is.na(Rn), `Target Name` == "SNCA") %>%
      # Use snca_well_map (not well_map): SNCA wells are purged from sample_map
      # so they only appear in the unpurged rescan above.
      left_join(snca_well_map, by = "Well") %>%
      filter(!is.na(Group_ID), !is.na(Bio_Rep)) %>%
      mutate(File = fname)

    if (is_rep5) snca_chunk <- snca_chunk %>% filter(Well != "C4")
    if (is_check) snca_chunk <- snca_chunk %>% filter(Group_ID == "G1", Timepoint == "4h")
    if (nrow(snca_chunk) > 0) {
      snca_amp <- bind_rows(snca_amp, snca_chunk)
      message("  snca amp rows: ", nrow(snca_chunk))
    }
  }

  # ── Melt Data ──────────────────────────────────────────────────────────────
  if (!is_check) {
    rt <- safe_read_sheet(fpath, "Melt Region Temperature Data",  skip = 7)
    rn <- safe_read_sheet(fpath, "Melt Region Normalized Data",   skip = 7)
    if (!is.null(rt) && !is.null(rn)) {
      melt_chunk <- inner_join(
        pivot_melt_sheet(rt, "Temperature"),
        pivot_melt_sheet(rn, "Raw_Fluorescence"),
        by = c("Well", "Gene", "Reading_Index")
      ) %>%
        rename(`Target Name` = Gene) %>%
        left_join(sample_map %>% select(Well, `Sample Name`, Group_ID, Bio_Rep, Timepoint) %>% distinct(),
                  by = "Well") %>%
        filter(!is.na(Group_ID)) %>%
        mutate(File = fname)
      master_melt <- bind_rows(master_melt, melt_chunk)
      message("  melt rows: ", nrow(melt_chunk))
    }
  }
}

message("\nTotal amp rows  : ", nrow(master_amp))
message("Total snca rows : ", nrow(snca_amp))
message("Total melt rows : ", nrow(master_melt))

# De-duplicate: Check-sang1 takes priority for G1/4h
dedup <- function(df) {
  df %>%
    mutate(Is_Check = grepl("Check-sang1", File, ignore.case = TRUE)) %>%
    group_by(Well, `Target Name`, Cycle, Group_ID, Bio_Rep, Timepoint) %>%
    arrange(desc(Is_Check), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-Is_Check)
}
master_amp <- dedup(master_amp)
if (nrow(snca_amp) > 0) snca_amp <- dedup(snca_amp)

# =============================================================================
# SS 4  STAGE 2 — MELT CURVE QC + AUTOMATED FLAGGING
# =============================================================================
message("\n=== STAGE 2: MELT CURVE QC & AUTOMATED FLAGGING ===")

compute_melt_deriv <- function(df) {
  df %>%
    filter(Temperature >= 65, Temperature <= 89) %>%
    group_by(Well, `Target Name`, `Sample Name`, Group_ID, Bio_Rep, Timepoint) %>%
    arrange(Temperature, .by_group = TRUE) %>%
    group_modify(function(d, keys) {
      if (nrow(d) < 5) return(tibble())
      sp <- tryCatch(smooth.spline(d$Temperature, d$Raw_Fluorescence, spar = 0.75),
                     error = function(e) NULL)
      if (is.null(sp)) return(tibble())
      d$Smoothed <- predict(sp, d$Temperature)$y
      d          <- arrange(d, Temperature)
      d$Deriv    <- c(-diff(d$Smoothed) / diff(d$Temperature), NA_real_)
      d
    }) %>%
    ungroup()
}

melt_deriv_bio <- tibble()
melt_flags     <- tibble()

if (nrow(master_melt) > 0) {

  melt_bio <- master_melt %>% filter(!is.na(Bio_Rep))

  if (nrow(melt_bio) > 0) {
    melt_deriv_bio <- compute_melt_deriv(melt_bio)

    # ── Interactive biological melt plot ──────────────────────────────────────
    bio_groups   <- sort(unique(melt_deriv_bio$Group_ID))
    grp_pal_melt <- setNames(
      colorRampPalette(c(
        "#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928",
        "#A6CEE3","#B2DF8A","#FB9A99","#FDBF6F","#CAB2D6","#FFFF99",
        "#8DD3C7","#FB8072"))(length(bio_groups)),
      bio_groups
    )
    p_melt_bio <- ggplotly(
      ggplot(melt_deriv_bio %>% filter(!is.na(Deriv)),
             aes(x = Temperature, y = Deriv, color = Group_ID,
                 group = interaction(Well, Bio_Rep),
                 text  = paste0("Group: ", Group_ID, "<br>Rep: ", Bio_Rep,
                                "<br>Temp: ", round(Temperature, 1)))) +
        geom_line(alpha = 0.7, linewidth = 0.45) +
        facet_wrap(~ `Target Name`, scales = "free_y", ncol = 3) +
        scale_color_manual(values = grp_pal_melt) +
        labs(title    = "Melt Peaks (-dF/dT) - Biological Samples",
             subtitle = "Colour = Group ID | Multiple lines = replicates",
             x = "Temperature (C)", y = "-dRn/dT", color = "Group") +
        theme_bw(base_size = 10) +
        theme(legend.position = "right", strip.text = element_text(face = "bold", size = 8)),
      tooltip = "text"
    )
    saveWidget(p_melt_bio, file.path(OUT_DIR, "01_interactive_melt_peaks.html"),
               selfcontained = TRUE)
    message("  Saved: 01_interactive_melt_peaks.html  (",
            n_distinct(melt_deriv_bio$Well), " wells)")

    # ── Automated Melt QC flagging ─────────────────────────────────────────────
    # Step 1: primary peak (Tm + height) per well x gene
    well_peaks <- melt_deriv_bio %>%
      filter(!is.na(Deriv), Deriv > 0) %>%
      group_by(Well, `Target Name`, Group_ID, Bio_Rep) %>%
      summarise(
        Peak_Tm     = Temperature[which.max(Deriv)],
        Peak_Height = max(Deriv, na.rm = TRUE),
        .groups     = "drop"
      )

    # Step 2: expected Tm = median primary Tm across all wells per gene
    expected_tm <- well_peaks %>%
      group_by(`Target Name`) %>%
      summarise(Expected_Tm = median(Peak_Tm, na.rm = TRUE), .groups = "drop")

    # Step 3: multiple-peak detection per well x gene
    multi_peak_flags <- melt_deriv_bio %>%
      filter(!is.na(Deriv)) %>%
      group_by(Well, `Target Name`, Group_ID, Bio_Rep) %>%
      summarise(
        Has_Multi_Peak = detect_multi_peak(Temperature, Deriv),
        .groups        = "drop"
      )

    # Step 4: combine criteria -> Exclude.amplicon
    melt_flags <- well_peaks %>%
      left_join(expected_tm,    by = "Target Name") %>%
      left_join(multi_peak_flags, by = c("Well", "Target Name", "Group_ID", "Bio_Rep")) %>%
      mutate(
        Tm_Deviation     = abs(Peak_Tm - Expected_Tm),
        # Exclusion rule: ANY of the three criteria fires -> flag = 1
        Exclude.amplicon = as.integer(
          Tm_Deviation  > TM_DEVIATION_TOL   |   # criterion 1: Tm shift
          Has_Multi_Peak                     |   # criterion 2: secondary band
          Peak_Height   < MIN_PEAK_HEIGHT        # criterion 3: flat/no product
        )
      ) %>%
      select(Well, `Target Name`, Group_ID, Bio_Rep,
             Peak_Tm, Expected_Tm, Tm_Deviation, Peak_Height,
             Has_Multi_Peak, Exclude.amplicon)

    n_flagged <- sum(melt_flags$Exclude.amplicon, na.rm = TRUE)
    message("  Melt QC flags: ", n_flagged, " / ", nrow(melt_flags), " wells excluded")
    if (n_flagged > 0) {
      message("  Flagged wells (Exclude.amplicon = 1):")
      melt_flags %>%
        filter(Exclude.amplicon == 1) %>%
        select(Well, `Target Name`, Group_ID, Bio_Rep,
               Peak_Tm, Expected_Tm, Tm_Deviation, Has_Multi_Peak) %>%
        print(n = 30)
    }

    # Save melt QC flags for audit
    write_csv(melt_flags,
              file.path(OUT_DIR, "00_melt_qc_flags.csv"))
    message("  Saved: 00_melt_qc_flags.csv")
  }

  # ── NTC contamination check ─────────────────────────────────────────────────
  melt_ntc <- master_melt %>% filter(Group_ID == "NTC")
  if (nrow(melt_ntc) > 0) {
    melt_deriv_ntc <- compute_melt_deriv(melt_ntc)
    ntc_targets    <- sort(unique(melt_deriv_ntc$`Target Name`))
    p_melt_ntc <- ggplotly(
      ggplot(melt_deriv_ntc %>% filter(!is.na(Deriv)),
             aes(x = Temperature, y = Deriv, color = `Target Name`,
                 group = interaction(Well, File),
                 text  = paste0("Target: ", `Target Name`,
                                "<br>Well: ", Well,
                                "<br>Temp: ", round(Temperature, 1)))) +
        geom_line(alpha = 0.75, linewidth = 0.5) +
        facet_wrap(~ `Target Name`, scales = "free_y", ncol = 3) +
        scale_color_manual(values = setNames(
          colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(ntc_targets)),
          ntc_targets)) +
        labs(title    = "Melt Peaks (NTC) - Negative Controls Only",
             subtitle = "Any clear peak = contamination",
             x = "Temperature (C)", y = "-dRn/dT") +
        theme_bw(base_size = 10) +
        theme(legend.position = "none", strip.text = element_text(face = "bold", size = 8)),
      tooltip = "text"
    )
    saveWidget(p_melt_ntc, file.path(OUT_DIR, "01b_ntc_melt_peaks.html"),
               selfcontained = TRUE)
    message("  Saved: 01b_ntc_melt_peaks.html")
  }

} else {
  message("  No melt data available — skipping QC plots and flagging.")
}

# =============================================================================
# SS 5  STAGE 3 — BASELINE CORRECTION + SDM Ct + AMPLIFICATION QC PLOT
# =============================================================================
message("\n=== STAGE 3: BASELINE CORRECTION & SDM Ct ===")

amp_corrected <- master_amp %>%
  group_by(Well, `Target Name`, Group_ID, Bio_Rep) %>%
  arrange(Cycle, .by_group = TRUE) %>%
  mutate(
    Baseline_Median       = median(Rn[Cycle >= 3 & Cycle <= 15], na.rm = TRUE),
    Baseline_Corrected_Rn = Rn - Baseline_Median
  ) %>%
  ungroup()

ct_sdm <- amp_corrected %>%
  group_by(Well, `Target Name`, Group_ID, Bio_Rep, Timepoint) %>%
  summarise(Ct_SDM = sdm_ct(Cycle, Baseline_Corrected_Rn), .groups = "drop")

# Attach Melt QC flag; wells with no melt data (Check-sang1) are treated as clean
ct_table <- ct_sdm %>%
  left_join(
    melt_flags %>% select(Well, `Target Name`, Group_ID, Bio_Rep, Exclude.amplicon),
    by = c("Well", "Target Name", "Group_ID", "Bio_Rep")
  ) %>%
  mutate(Exclude.amplicon = coalesce(Exclude.amplicon, 0L))

message("  Ct_SDM valid: ", sum(!is.na(ct_table$Ct_SDM)),
        "  range: [", round(min(ct_table$Ct_SDM, na.rm = TRUE), 1),
        ", ",         round(max(ct_table$Ct_SDM, na.rm = TRUE), 1), "]")
message("  Exclude.amplicon=1: ", sum(ct_table$Exclude.amplicon, na.rm = TRUE),
        " wells  |  clean: ", sum(ct_table$Exclude.amplicon == 0L, na.rm = TRUE))

# Inter-replicate consistency
ct_consistency <- ct_table %>%
  group_by(Group_ID, Timepoint, `Target Name`) %>%
  summarise(
    N_total  = n(),
    N_clean  = sum(Exclude.amplicon == 0L),
    Mean_SDM = round(mean(Ct_SDM, na.rm = TRUE), 2),
    SD_SDM   = round(sd(Ct_SDM,   na.rm = TRUE), 3),
    .groups  = "drop"
  ) %>%
  arrange(Timepoint, `Target Name`, Group_ID)
message("\n  Inter-replicate Ct_SDM consistency:")
print(ct_consistency, n = 30)

# Smooth display curves
amp_display <- amp_corrected %>%
  group_by(Well, `Target Name`, Group_ID, Bio_Rep, Timepoint) %>%
  arrange(Cycle, .by_group = TRUE) %>%
  group_modify(function(df, keys) {
    if (nrow(df) < 6) return(df %>% mutate(Rn_Smooth = Baseline_Corrected_Rn))
    sp <- tryCatch(smooth.spline(df$Cycle, df$Baseline_Corrected_Rn, spar = 0.55),
                   error = function(e) NULL)
    df$Rn_Smooth <- if (!is.null(sp)) predict(sp, df$Cycle)$y else df$Baseline_Corrected_Rn
    df
  }) %>%
  ungroup()

# SDM marker points for display
sdm_pts <- amp_corrected %>%
  inner_join(ct_sdm %>% select(Well, `Target Name`, Ct_SDM),
             by = c("Well", "Target Name"), relationship = "many-to-many") %>%
  group_by(Well, `Target Name`, Group_ID, Bio_Rep, Timepoint, Ct_SDM) %>%
  group_modify(function(df, keys) {
    sp <- tryCatch(smooth.spline(df$Cycle, df$Baseline_Corrected_Rn, spar = 0.5),
                   error = function(e) NULL)
    tibble(SDM_Rn = if (!is.null(sp)) predict(sp, keys$Ct_SDM)$y else NA_real_)
  }) %>%
  ungroup() %>%
  filter(!is.na(SDM_Rn))

# Interactive amplification QC plot
all_amp_groups <- sort(unique(amp_display$Group_ID))
grp_pal_amp <- setNames(
  colorRampPalette(c(
    "#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928",
    "#A6CEE3","#B2DF8A","#FB9A99","#FDBF6F","#CAB2D6","#FFFF99",
    "#8DD3C7","#FB8072"))(length(all_amp_groups)),
  all_amp_groups
)

p_amp <- suppressWarnings(ggplotly(
  ggplot(amp_display %>% filter(Cycle <= 40),
         aes(x = Cycle, y = Rn_Smooth, color = Group_ID,
             group = interaction(Well, Bio_Rep),
             text  = paste0("Group: ", Group_ID, "<br>Rep: ", Bio_Rep,
                            "<br>Cycle: ", Cycle,
                            "<br>Rn: ", round(Rn_Smooth, 4)))) +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    geom_point(
      data = sdm_pts,
      aes(x = Ct_SDM, y = SDM_Rn, color = Group_ID,
          text = paste0("Group: ", Group_ID, "<br>Rep: ", Bio_Rep,
                        "<br>Ct_SDM: ", round(Ct_SDM, 2))),
      shape = 25, size = 1.8, fill = "#FF8C00", inherit.aes = FALSE
    ) +
    facet_grid(Timepoint ~ `Target Name`, scales = "free_y") +
    scale_color_manual(values = grp_pal_amp) +
    labs(title    = "Amplification Curves - All Groups & Replicates (SDM)",
         subtitle = "Colour = Group | Orange triangle = SDM Ct",
         x = "Cycle", y = "Baseline-Corrected Rn (smoothed)", color = "Group") +
    theme_bw(base_size = 9) +
    theme(legend.position = "right",
          strip.text.x    = element_text(size = 7, face = "bold"),
          strip.text.y    = element_text(size = 8, face = "bold")),
  tooltip = "text"
)) %>% layout(legend = list(orientation = "v"))

saveWidget(p_amp, file.path(OUT_DIR, "02_interactive_amplification_curves.html"),
           selfcontained = TRUE)
message("  Saved: 02_interactive_amplification_curves.html")

# =============================================================================
# SS 6  STAGE 4 — DUAL PFAFFL RUNS (CLEAN-ONLY vs. ALL-INCLUSIVE)
# =============================================================================
message("\n=== STAGE 4: DUAL PFAFFL EXPRESSION QUANTIFICATION ===")

# Average technical replicates -> one Mean_Ct per Bio_Rep x Group x Gene
make_avg_ct <- function(tbl) {
  tbl %>%
    filter(!is.na(Ct_SDM)) %>%
    group_by(Group_ID, Timepoint, `Target Name`, Bio_Rep) %>%
    summarise(Mean_Ct = mean(Ct_SDM, na.rm = TRUE), .groups = "drop")
}

# Two input sets differing only by the Exclude.amplicon filter
avg_sdm_clean <- make_avg_ct(ct_table %>% filter(Exclude.amplicon == 0L))
avg_sdm_all   <- make_avg_ct(ct_table)

message("  avg_sdm_clean rows: ", nrow(avg_sdm_clean),
        "  |  avg_sdm_all rows: ", nrow(avg_sdm_all))

compute_fc <- function(avg_tbl, label) {
  bind_rows(
    avg_tbl %>% filter(Timepoint == "4h") %>% pfaffl_per_rep(CTRL_4H),
    avg_tbl %>% filter(Timepoint == "8h") %>% pfaffl_per_rep(CTRL_8H)
  ) %>%
    left_join(avg_tbl %>% select(Group_ID, Bio_Rep, Timepoint) %>% distinct(),
              by = c("Group_ID", "Bio_Rep")) %>%
    mutate(Dataset = label)
}

fc_clean <- compute_fc(avg_sdm_clean, "Clean (Exclude.amplicon = 0)")
fc_all   <- compute_fc(avg_sdm_all,   "All Inclusive")

ct_clean_wells <- ct_table %>% filter(Exclude.amplicon == 0L)

fc_well_level <- bind_rows(
  pfaffl_well_level(ct_clean_wells %>% filter(Timepoint == "4h"), CTRL_4H),
  pfaffl_well_level(ct_clean_wells %>% filter(Timepoint == "8h"), CTRL_8H)
)

# Pull Batch (source file) and Sample Name from master_amp for full traceability
well_meta <- master_amp %>%
  select(Well, `Target Name`, Group_ID, Bio_Rep, `Sample Name`, File) %>%
  rename(Batch = File) %>%
  distinct()

per_rep_export <- fc_well_level %>%
  left_join(well_meta, by = c("Well", "Target Name", "Group_ID", "Bio_Rep")) %>%
  rename(Gene = `Target Name`) %>%
  select(Batch, Well, Gene, `Sample Name`, Group_ID, Bio_Rep, Timepoint,
         Exclude.amplicon, Ct_SDM, FoldChange) %>%
  filter(!is.na(FoldChange)) %>%
  arrange(Timepoint, Group_ID, Gene, Bio_Rep)

write_csv(per_rep_export, file.path(OUT_DIR, "04_per_rep_fc_sdm.csv"))
message("  Saved: 04_per_rep_fc_sdm.csv  (",
        nrow(per_rep_export), " well-level rows | ",
        n_distinct(per_rep_export$Gene), " genes | ",
        n_distinct(paste(per_rep_export$Group_ID, per_rep_export$Bio_Rep)),
        " unique bio reps)")

make_bio_stats <- function(fc_df) {
  fc_df %>%
    group_by(Dataset, Timepoint, Group_ID, `Target Name`) %>%
    summarise(
      N_reps  = sum(!is.na(FoldChange)),
      Mean_FC = mean(FoldChange, na.rm = TRUE),
      SD_FC   = sd(FoldChange,   na.rm = TRUE),
      SEM_FC  = SD_FC / sqrt(max(N_reps, 1L)),
      .groups = "drop"
    ) %>%
    arrange(Timepoint, `Target Name`, Group_ID)
}

bio_stats_clean <- make_bio_stats(fc_clean)
bio_stats_all   <- make_bio_stats(fc_all)

# Ct summary for CSV
ct_bio <- ct_table %>%
  group_by(Group_ID, Timepoint, `Target Name`) %>%
  summarise(
    N_total     = n(),
    N_clean     = sum(Exclude.amplicon == 0L),
    Mean_Ct_SDM = round(mean(Ct_SDM, na.rm = TRUE), 3),
    SD_Ct_SDM   = round(sd(Ct_SDM,   na.rm = TRUE), 3),
    .groups     = "drop"
  )

# Combined expression summary CSV
bind_rows(bio_stats_clean, bio_stats_all) %>%
  left_join(ct_bio, by = c("Group_ID", "Timepoint", "Target Name")) %>%
  arrange(Dataset, Timepoint, `Target Name`, Group_ID) %>%
  write_csv(file.path(OUT_DIR, "03_master_biological_expression_summary.csv"))
message("  Saved: 03_master_biological_expression_summary.csv")

# =============================================================================
# SS 7  STAGE 5 — TWO EXPRESSION PDFs (clean-only + all-inclusive)
# =============================================================================
message("\n=== STAGE 5: PUBLICATION FIGURES ===")

grp_levels    <- paste0("G", 1:14)
target_levels <- sort(unique(bio_stats_clean$`Target Name`[
  !(bio_stats_clean$`Target Name` %in% REFERENCE_GENES)]))

# Color palette: 4h (G1-G7) = blues, 8h (G8-G14) = reds
pal_4h  <- colorRampPalette(c("#08306B","#2171B5","#4292C6",
                               "#6BAED6","#9ECAE1","#C6DBEF","#DEEBF7"))(7)
pal_8h  <- colorRampPalette(c("#67000D","#A50026","#D73027",
                               "#F46D43","#FDAE61","#FEE090","#FFFFBF"))(7)
grp_pal <- setNames(c(pal_4h, pal_8h), grp_levels)

bar_theme <- theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "#EBEBEB"),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    panel.spacing.x  = unit(0.3, "lines"),
    panel.spacing.y  = unit(0.6, "lines"),
    legend.position  = "none"
  )

n_targets  <- length(target_levels)
n_cols_pdf <- 3L
pdf_w      <- n_cols_pdf * 5.5
pdf_h      <- ceiling(n_targets / n_cols_pdf) * 3.8 + 2.5

# Generic bar + jitter page builder
make_bar_page <- function(stats_df, reps_df, page_title, page_subtitle) {

  plot_data <- stats_df %>%
    filter(!(`Target Name` %in% REFERENCE_GENES),
           !is.na(Mean_FC), is.finite(Mean_FC)) %>%
    mutate(Group_ID      = factor(Group_ID,      levels = grp_levels),
           `Target Name` = factor(`Target Name`, levels = target_levels))

  rep_data <- reps_df %>%
    filter(!(`Target Name` %in% REFERENCE_GENES),
           !is.na(FoldChange), is.finite(FoldChange)) %>%
    mutate(Group_ID      = factor(Group_ID,      levels = grp_levels),
           `Target Name` = factor(`Target Name`, levels = target_levels))

  ggplot(plot_data, aes(x = Group_ID, y = Mean_FC, fill = Group_ID)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_col(width = 0.72, alpha = 0.88) +
    geom_errorbar(
      aes(ymin = pmax(Mean_FC - SEM_FC, 0), ymax = Mean_FC + SEM_FC),
      width = 0.30, linewidth = 0.5, color = "grey20"
    ) +
    geom_point(
      data        = rep_data,
      aes(x = Group_ID, y = FoldChange),
      position    = position_jitter(width = 0.15, seed = 42),
      size        = 1.5, alpha = 0.7, color = "grey15",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ `Target Name`, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = grp_pal) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(title    = page_title,
         subtitle = page_subtitle,
         x = "Treatment Group",
         y = "Fold Change (Pfaffl geNORM)") +
    bar_theme
}

sub_common <- paste0(
  "G1-G7 = 4h  |  G8-G14 = 8h  |  Calibrators: ", CTRL_4H, " (4h) & ", CTRL_8H, " (8h)",
  "  |  Ref: ", paste(REFERENCE_GENES, collapse = " + "),
  " (geNORM geometric mean)  |  Bars = mean +/- SEM  |  Dots = individual reps"
)

p_clean <- make_bar_page(
  bio_stats_clean, fc_clean,
  "[A] Pfaffl SDM Ct - CLEAN WELLS ONLY (Exclude.amplicon = 0)",
  paste0(sub_common, "  |  Flagged wells removed")
)
p_all <- make_bar_page(
  bio_stats_all, fc_all,
  "[B] Pfaffl SDM Ct - ALL INCLUSIVE (unfiltered)",
  paste0(sub_common, "  |  All wells included regardless of Melt QC flag")
)

pdf(file.path(OUT_DIR, "04a_expression_clean_only.pdf"),
    width = pdf_w, height = pdf_h, onefile = FALSE)
print(p_clean)
dev.off()
message("  Saved: 04a_expression_clean_only.pdf")

pdf(file.path(OUT_DIR, "04b_expression_all_inclusive.pdf"),
    width = pdf_w, height = pdf_h, onefile = FALSE)
print(p_all)
dev.off()
message("  Saved: 04b_expression_all_inclusive.pdf")

# =============================================================================
# SS 8  STAGE 6 — SNCA EXPLORATORY EXPRESSION
# =============================================================================
message("\n=== STAGE 6: SNCA EXPLORATORY ANALYSIS ===")

if (nrow(snca_amp) == 0) {
  message("  No SNCA amplification data found — skipping exploratory plot.")
} else {
  message("  SNCA amp rows: ", nrow(snca_amp),
          "  |  E_SNCA = ", get_eff("SNCA"), " (pilot; update after standard curve)")

  # Baseline correction
  snca_corrected <- snca_amp %>%
    group_by(Well, `Target Name`, Group_ID, Bio_Rep) %>%
    arrange(Cycle, .by_group = TRUE) %>%
    mutate(
      Baseline_Median       = median(Rn[Cycle >= 3 & Cycle <= 15], na.rm = TRUE),
      Baseline_Corrected_Rn = Rn - Baseline_Median
    ) %>%
    ungroup()

  # SDM Ct for SNCA
  snca_ct <- snca_corrected %>%
    group_by(Well, `Target Name`, Group_ID, Bio_Rep, Timepoint) %>%
    summarise(Ct_SDM = sdm_ct(Cycle, Baseline_Corrected_Rn), .groups = "drop")

  # Average technical reps
  snca_avg <- snca_ct %>%
    filter(!is.na(Ct_SDM)) %>%
    group_by(Group_ID, Timepoint, `Target Name`, Bio_Rep) %>%
    summarise(Mean_Ct = mean(Ct_SDM, na.rm = TRUE), .groups = "drop")

  message("  SNCA Ct_SDM valid: ", nrow(snca_avg),
          "  range: [", round(min(snca_avg$Mean_Ct, na.rm = TRUE), 1),
          ", ",         round(max(snca_avg$Mean_Ct, na.rm = TRUE), 1), "]")

  # Combine SNCA Cts with reference gene Cts from the all-inclusive avg_sdm_all.
  # Reference gene rows supply GAPDH and RPS11 Cts needed for geNORM normalization.
  ref_cts <- avg_sdm_all %>%
    filter(`Target Name` %in% REFERENCE_GENES)

  snca_combined <- bind_rows(snca_avg, ref_cts)

  # Pfaffl FC (SNCA treated as target; references supply NF via pfaffl_per_rep)
  snca_fc <- bind_rows(
    snca_combined %>% filter(Timepoint == "4h") %>% pfaffl_per_rep(CTRL_4H),
    snca_combined %>% filter(Timepoint == "8h") %>% pfaffl_per_rep(CTRL_8H)
  ) %>%
    left_join(snca_combined %>% select(Group_ID, Bio_Rep, Timepoint) %>% distinct(),
              by = c("Group_ID", "Bio_Rep")) %>%
    filter(`Target Name` == "SNCA")

  snca_stats <- snca_fc %>%
    group_by(Timepoint, Group_ID, `Target Name`) %>%
    summarise(
      N_reps  = sum(!is.na(FoldChange)),
      Mean_FC = mean(FoldChange, na.rm = TRUE),
      SD_FC   = sd(FoldChange,   na.rm = TRUE),
      SEM_FC  = SD_FC / sqrt(max(N_reps, 1L)),
      .groups = "drop"
    )

  snca_plot_df <- snca_stats %>%
    filter(!is.na(Mean_FC), is.finite(Mean_FC)) %>%
    mutate(Group_ID = factor(Group_ID, levels = grp_levels))

  snca_rep_df <- snca_fc %>%
    filter(!is.na(FoldChange), is.finite(FoldChange)) %>%
    mutate(Group_ID = factor(Group_ID, levels = grp_levels))

  p_snca <- ggplot(snca_plot_df, aes(x = Group_ID, y = Mean_FC, fill = Group_ID)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_col(width = 0.72, alpha = 0.88) +
    geom_errorbar(
      aes(ymin = pmax(Mean_FC - SEM_FC, 0), ymax = Mean_FC + SEM_FC),
      width = 0.30, linewidth = 0.5, color = "grey20"
    ) +
    geom_point(
      data        = snca_rep_df,
      aes(x = Group_ID, y = FoldChange),
      position    = position_jitter(width = 0.15, seed = 42),
      size        = 2.0, alpha = 0.75, color = "grey15",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = grp_pal) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(
      title    = "SNCA - Exploratory Relative Expression (SDM Ct, E = 2.0)",
      subtitle = paste0(
        "EXPLORATORY - no valid SNCA standard curve (E = 2.0 assumed)  |  G1-G7 = 4h  |  G8-G14 = 8h",
        "  |  Calibrators: ", CTRL_4H, " (4h) & ", CTRL_8H, " (8h)",
        "  |  Ref: ", paste(REFERENCE_GENES, collapse = " + "),
        " (geNORM NF)  |  Bars = mean +/- SEM  |  Dots = individual reps"
      ),
      x = "Treatment Group",
      y = "Fold Change (Pfaffl geNORM, SNCA E = 2.0)"
    ) +
    bar_theme +
    theme(legend.position = "none")

  pdf(file.path(OUT_DIR, "05_exploratory_snca_trend.pdf"),
      width = 10, height = 5.5, onefile = FALSE)
  print(p_snca)
  dev.off()
  message("  Saved: 05_exploratory_snca_trend.pdf")

  # Save SNCA Ct and FC summary
  snca_summary <- snca_stats %>%
    left_join(
      snca_ct %>%
        group_by(Group_ID, Timepoint, `Target Name`) %>%
        summarise(N_wells    = n(),
                  Mean_Ct    = round(mean(Ct_SDM, na.rm = TRUE), 3),
                  SD_Ct      = round(sd(Ct_SDM,   na.rm = TRUE), 3),
                  .groups    = "drop"),
      by = c("Group_ID", "Timepoint", "Target Name")
    ) %>%
    mutate(across(c(Mean_FC, SEM_FC), ~ round(.x, 3))) %>%
    arrange(Timepoint, Group_ID)

  write_csv(snca_summary,
            file.path(OUT_DIR, "05_snca_exploratory_summary.csv"))
  message("  Saved: 05_snca_exploratory_summary.csv")
}

# =============================================================================
# SS 9  CONSOLE SUMMARY
# =============================================================================
message("\n=== PIPELINE COMPLETE ===")
message("  Output directory: ", OUT_DIR)
message("  Files generated:")
walk(list.files(OUT_DIR, full.names = FALSE), ~ message("    ", .x))

message("\n  Clean expression summary (first 20 rows):")
print(
  bio_stats_clean %>%
    filter(!(`Target Name` %in% REFERENCE_GENES)) %>%
    select(Timepoint, Group_ID, `Target Name`, N_reps, Mean_FC, SEM_FC) %>%
    mutate(across(c(Mean_FC, SEM_FC), ~ round(.x, 3))) %>%
    head(20),
  n = 20
)
