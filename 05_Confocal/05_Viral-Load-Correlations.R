# ==============================================================================
# 05_Viral-Load-Correlations.R
# ==============================================================================

library(tidyverse)

# ─── CONFIGS ───
cells_csv <- file.path("Confocal", "results_segmentation", "master_cells.csv")
out_csv   <- file.path("Confocal", "results_confocal_stats", "data", "Figure5_PartialR.csv")

message("Loading cells data...")
if (!file.exists(cells_csv)) {
  stop("Cells CSV file not found: ", cells_csv)
}
cells_raw <- read_csv(cells_csv, show_col_types = FALSE)

cells <- cells_raw %>%
  mutate(
    antibody_clean = case_when(
      grepl("Syn1",              antibody, ignore.case = TRUE) ~ "Syn1",
      grepl("EPR|Abcam|C-?term", antibody, ignore.case = TRUE) ~ "EPR20535",
      TRUE ~ as.character(antibody)
    )
  )

configs <- list(
  list(name = "NC_ratio",       ab = "Syn1",     col = "NC_ratio"),
  list(name = "asyn_nuc_Syn1",  ab = "Syn1",     col = "asyn_nuclear_mean"),
  list(name = "asyn_cyto_Syn1", ab = "Syn1",     col = "asyn_cytoplasmic_mean"),
  list(name = "asyn_cyto_EPR", ab = "EPR20535", col = "asyn_cytoplasmic_mean"),
  list(name = "asyn_nuc_EPR",  ab = "EPR20535", col = "asyn_nuclear_mean")
)

results <- list()

for (cfg in configs) {
  for (vir in c("V1", "V2")) {
    sub_df <- cells %>%
      filter(
        condition == vir,
        antibody_clean == cfg$ab,
        as.logical(infected) == TRUE,
        !is.na(.data[[cfg$col]]),
        !is.na(iav_mean),
        iav_mean > 0
      )
    
    n_cells <- nrow(sub_df)
    n_fovs <- if (n_cells > 0) n_distinct(sub_df$image_id) else 0
    
    part_r  <- NA_real_
    ci_low  <- NA_real_
    ci_high <- NA_real_
    p_val   <- NA_real_
    note_msg <- ""
    
    if (n_cells < 15) {
      note_msg <- "n<15 insufficient"
    } else {
      x <- log10(sub_df$iav_mean)
      y <- sub_df[[cfg$col]]
      
      if (n_fovs <= 1) {
        test_res <- cor.test(x, y, method = "pearson")
        part_r  <- as.numeric(test_res$estimate)
        ci_low  <- as.numeric(test_res$conf.int[1])
        ci_high <- as.numeric(test_res$conf.int[2])
        p_val   <- as.numeric(test_res$p.value)
        note_msg <- "n_FOVs=1, raw Pearson (not partial)"
      } else if (n_cells - n_fovs - 1 <= 0) {
        note_msg <- "insufficient degrees of freedom for partial correlation"
      } else {
        res_y <- residuals(lm(y ~ image_id, data = sub_df))
        res_x <- residuals(lm(x ~ image_id, data = sub_df))
        
        part_r <- cor(res_x, res_y)
        
        r_clamped <- max(min(part_r, 0.9999), -0.9999)
        z  <- atanh(r_clamped)
        se <- 1 / sqrt(n_cells - n_fovs - 1)
        
        ci_low  <- tanh(z - 1.96 * se)
        ci_high <- tanh(z + 1.96 * se)
        p_val   <- 2 * (1 - pnorm(abs(z) / se))
        note_msg <- "partial Pearson controlled for image_id"
      }
    }
    
    results[[length(results) + 1]] <- tibble(
      Metric = cfg$name,
      Antibody_filter = cfg$ab,
      Virus = vir,
      n_cells = n_cells,
      n_FOVs = n_fovs,
      partial_r = part_r,
      r_CI_low = ci_low,
      r_CI_high = ci_high,
      p_value = p_val,
      Note = note_msg
    )
  }
}

res_df <- bind_rows(results)

out_dir <- dirname(out_csv)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

write_csv(res_df, out_csv)
message("\n✓ Correlation analysis completed successfully.")
message("  Saved results to: ", out_csv)

cat("\n======================================================================\n")
cat("                   FIGURE 5 PARTIAL CORRELATION SUMMARY               \n")
cat("======================================================================\n")
print(as.data.frame(res_df))
cat("======================================================================\n")
