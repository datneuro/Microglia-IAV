

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsignif)
  library(patchwork)
  library(ggridges)
  library(openxlsx)
  library(ggpubr)
  library(png)
  library(grid)
  library(cowplot)
})

# Inputs are the per-cell table and the statistics workbook written by
# 05_Confocal/. Paths may be given on the command line (see below) or through
# the CONFOCAL_OUT environment variable.
CONFOCAL_DIR <- Sys.getenv("CONFOCAL_DIR", unset = "Confocal")
CONFOCAL_OUT <- Sys.getenv("CONFOCAL_OUT", unset = file.path(CONFOCAL_DIR, "results_confocal_stats"))
SEG_DIR      <- file.path(CONFOCAL_DIR, "results_segmentation")
args <- commandArgs(trailingOnly = TRUE)
cells_csv   <- file.path(SEG_DIR, "master_cells.csv")
coloc_csv   <- file.path(SEG_DIR, "master_coloc.csv")
stats_dir   <- CONFOCAL_OUT
exclude_csv <- NULL

for (i in seq_along(args)) {
  if (args[i] == "--master-cells" && i < length(args)) cells_csv <- args[i + 1]
  if (args[i] == "--stats-dir"    && i < length(args)) stats_dir <- args[i + 1]
  if (args[i] == "--master-coloc" && i < length(args)) coloc_csv <- args[i + 1]
  if (args[i] == "--exclude"      && i < length(args)) exclude_csv <- args[i + 1]
}

# Locate stats Excel file

out_dir_final <- Sys.getenv("FIGURE_OUT", unset = file.path("Figures", "output"))
dir.create(out_dir_final, showWarnings = FALSE, recursive = TRUE)
show_scatter_stats <- FALSE
show_cd_stats <- FALSE
stats_files <- list.files(file.path(stats_dir, "data"), pattern = "^Confocal_Statistics_LMER.*\\.xlsx$", full.names = TRUE)
if (length(stats_files) == 0) stats_files <- list.files(file.path(stats_dir, "data"), pattern = "\\.xlsx$", full.names = TRUE)
stats_xlsx <- stats_files[which.max(file.info(stats_files)$mtime)]
# ─── COLOR PALETTES ───
strain_colors <- c(
    "A/PR/8/34" = "#4285F4",     # Bright blue
    "A/NWS/33"  = "#C0392B"      # Deep red
)

cond_state_colors <- c(
    "Naive"               = "#D5D8DC",   # Light gray
    "A/PR/8/34.Bystander" = "#A8C7F0",   # Pale blue (V1 Bys)
    "A/PR/8/34.Infected"  = "#4285F4",   # Bright blue (V1 Inf)
    "A/NWS/33.Bystander"  = "#F2B4B4",   # Pale red (V2 Bys)
    "A/NWS/33.Infected"   = "#C0392B"    # Deep red (V2 Inf)
)

condition_colors <- c(
    "A/PR/8/34" = "#4285F4",     # V1 Blue
    "A/NWS/33"  = "#C0392B"      # V2 Red
)

# ─── HELPER FUNCTIONS ───
p_to_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.0001) return("****")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("")
}

rename_condition <- function(x) {
    case_when(
        x == "V1" ~ "A/PR/8/34",
        x == "V2" ~ "A/NWS/33",
        TRUE ~ x
    )
}

rename_antibody <- function(x) {
    case_when(
        x == "Syn1" ~ "α-Synuclein (Syn1)",
        x == "EPR20535" ~ "α-Synuclein (EPR20535)",
        TRUE ~ x
    )
}

draw_key_split <- function(data, params, size) {
    fill_color <- data$fill
    if (is.na(fill_color)) {
        c1 <- "grey90"
        c2 <- "grey60"
    } else if (fill_color == "#4285F4" || fill_color == "#A8C7F0") {
        c1 <- "#A8C7F0"
        c2 <- "#4285F4"
    } else if (fill_color == "#C0392B" || fill_color == "#F2B4B4") {
        c1 <- "#F2B4B4"
        c2 <- "#C0392B"
    } else if (fill_color == "#D5D8DC" || fill_color == "Control" || fill_color == "Naive") {
        c1 <- "#E5E7E9"
        c2 <- "#979A9A"
    } else {
        c1 <- fill_color
        c2 <- fill_color
    }
    
    grid::grobTree(
        grid::polygonGrob(x = c(0, 0, 1), y = c(0, 1, 1), gp = grid::gpar(fill = c1, col = NA)),
        grid::polygonGrob(x = c(0, 1, 1), y = c(0, 0, 1), gp = grid::gpar(fill = c2, col = NA)),
        grid::rectGrob(gp = grid::gpar(fill = NA, col = "black", lwd = 0.5))
    )
}

# Grouped x-axis lines and labels
add_grouped_x_axis <- function(plot_obj, y_min, y_max, is_five_groups = FALSE) {
    line_grob_1 <- grid::segmentsGrob(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = unit(0, "npc") - unit(18, "points"),
        y1 = unit(0, "npc") - unit(18, "points"),
        gp = grid::gpar(col = "black", lwd = 1.2)
    )
    text_grob_1 <- grid::textGrob(
        label = "A/PR/8/34",
        x = unit(0.5, "npc"),
        y = unit(0, "npc") - unit(30, "points"),
        gp = grid::gpar(fontface = "bold", fontsize = 7.2, col = "black")
    )
    
    line_grob_2 <- grid::segmentsGrob(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = unit(0, "npc") - unit(18, "points"),
        y1 = unit(0, "npc") - unit(18, "points"),
        gp = grid::gpar(col = "black", lwd = 1.2)
    )
    text_grob_2 <- grid::textGrob(
        label = "A/NWS/33",
        x = unit(0.5, "npc"),
        y = unit(0, "npc") - unit(30, "points"),
        gp = grid::gpar(fontface = "bold", fontsize = 7.2, col = "black")
    )
    
    if (is_five_groups) {
        plot_obj <- plot_obj +
            annotation_custom(line_grob_1, xmin = 1.8, xmax = 3.2, ymin = -Inf, ymax = Inf) +
            annotation_custom(text_grob_1, xmin = 1.8, xmax = 3.2, ymin = -Inf, ymax = Inf) +
            annotation_custom(line_grob_2, xmin = 3.8, xmax = 5.2, ymin = -Inf, ymax = Inf) +
            annotation_custom(text_grob_2, xmin = 3.8, xmax = 5.2, ymin = -Inf, ymax = Inf)
    } else {
        plot_obj <- plot_obj +
            annotation_custom(line_grob_1, xmin = 0.8, xmax = 2.2, ymin = -Inf, ymax = Inf) +
            annotation_custom(text_grob_1, xmin = 0.8, xmax = 2.2, ymin = -Inf, ymax = Inf) +
            annotation_custom(line_grob_2, xmin = 2.8, xmax = 4.2, ymin = -Inf, ymax = Inf) +
            annotation_custom(text_grob_2, xmin = 2.8, xmax = 4.2, ymin = -Inf, ymax = Inf)
    }
    
    plot_obj + 
        coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
        theme(plot.margin = margin(t = 15, r = 5, b = 40, l = 5, unit = "pt"))
}

# ─── LOAD & PREPROCESS DATA ───
message("Loading cells data...")
cells_raw <- read_csv(cells_csv, show_col_types = FALSE)

# Auto-merge coloc file if columns are missing
if (!all(c("nuc_asyn_iav_costes_p", "nuc_asyn_iav_pearson_r", "nuc_asyn_dapi_pearson_r") %in% colnames(cells_raw))) {
  default_coloc <- coloc_csv
  if (file.exists(default_coloc)) {
    message("Coloc columns missing. Auto-merging from: ", default_coloc)
    coloc <- read_csv(default_coloc, show_col_types = FALSE)
    coloc_cols <- setdiff(colnames(coloc), c(colnames(cells_raw), "nucleus_id", "image_id"))
    cells_raw <- cells_raw %>%
      left_join(coloc %>% select(nucleus_id, image_id, all_of(coloc_cols)),
                by = c("nucleus_id", "image_id"))
  } else {
    warning("Coloc columns are missing and default master_coloc.csv could not be found.")
  }
}

if (!is.null(exclude_csv) && file.exists(exclude_csv)) {
  excl <- read_csv(exclude_csv, show_col_types = FALSE)
  n_before <- nrow(cells_raw)
  cells_raw <- cells_raw %>% filter(!image_id %in% excl$image_id)
  message("Excluded ", n_before - nrow(cells_raw), " cells from flagged FOVs.")
}

cells <- cells_raw %>%
    mutate(
        antibody = case_when(
            grepl("Syn1", antibody, ignore.case = TRUE) ~ "Syn1",
            grepl("EPR|Abcam|C-?term", antibody, ignore.case = TRUE) ~ "EPR20535",
            TRUE ~ as.character(antibody)
        ),
        condition = rename_condition(condition),
        condition = factor(condition, levels = c("Control", "A/PR/8/34", "A/NWS/33")),
        Antibody = rename_antibody(antibody),
        Antibody = factor(Antibody, levels = c("α-Synuclein (Syn1)", "α-Synuclein (EPR20535)")),
        IAV_state = factor(if_else(as.logical(infected), "Infected", "Bystander"),
            levels = c("Bystander", "Infected")
        ),
        cond_state = case_when(
            condition == "Control" ~ "Naive",
            TRUE ~ paste(condition, IAV_state, sep = ".")
        ),
        cond_state = factor(cond_state, levels = c("Naive", "A/PR/8/34.Bystander", "A/PR/8/34.Infected", "A/NWS/33.Bystander", "A/NWS/33.Infected"))
    ) %>%
    filter(!is.na(condition), !is.na(Antibody), !grepl("PFF", condition))

# ─── LOAD STATS RESULTS ───
m2v1_pw   <- read.xlsx(stats_xlsx, sheet = "M2_V1_Pairwise")
m2v2_pw   <- read.xlsx(stats_xlsx, sheet = "M2_V2_Pairwise")
m3b_anova <- read.xlsx(stats_xlsx, sheet = "M3B_ANOVA")
m4a_v1    <- tryCatch(read.xlsx(stats_xlsx, sheet = "M4A_V1_Slopes"), error = function(e) NULL)
m4a_v2    <- tryCatch(read.xlsx(stats_xlsx, sheet = "M4A_V2_Slopes"), error = function(e) NULL)
m4b_diff  <- tryCatch(read.xlsx(stats_xlsx, sheet = "M4B_SlopeDiff"), error = function(e) NULL)

# Helper function to get simple effect stars
get_m2_sig_label <- function(df_pw, metric_name, ab_name) {
    if (grepl("\\[", metric_name)) {
        row <- df_pw %>% filter(Metric == metric_name)
    } else {
        row <- df_pw %>% filter(Metric == metric_name, Antibody == ab_name)
    }
    if (nrow(row) == 0) return("")
    sig <- row$Significance[1]
    if (is.na(sig) || sig == "ns") return("")
    return(sig)
}

# Helper function to get simple effect p_adj and p_boot
get_m2_stats <- function(df_pw, metric_name, ab_name) {
    if (grepl("\\[", metric_name)) {
        row <- df_pw %>% filter(Metric == metric_name)
    } else {
        row <- df_pw %>% filter(Metric == metric_name, Antibody == ab_name)
    }
    if (nrow(row) == 0) return(list(p_adj = NA_real_, p_boot = NA_real_))
    return(list(p_adj = row$P_adjusted[1], p_boot = row$P_boot[1]))
}

# Formatting function for p-values
fmt_p_val <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("<0.001")
    return(sprintf("%.3f", p))
}

# Helper function to get Interaction stats
get_m3b_interaction <- function(df_anova, metric_val) {
    row <- df_anova %>% filter(Metric == metric_val, Term == "condition:IAV_state")
    if (nrow(row) == 0) return("")
    
    f_val <- as.numeric(row$F.value[1])
    p_val <- as.numeric(row$`Pr(>F)`[1])
    
    f_fmt <- if (is.na(f_val)) "NA" else sprintf("%.2f", f_val)
    
    return(sprintf("Interaction F=%s (p=%s)", f_fmt, fmt_p_val(p_val)))
}

# Helper to construct subtitle for bar plots incorporating bootstrap p-values & interaction F-stat
get_bar_subtitle <- function(df_pw_v1, df_pw_v2, df_anova, metric_filter, ab_short, show_interaction = TRUE) {
    stats_v1 <- get_m2_stats(df_pw_v1, metric_filter, ab_short)
    stats_v2 <- get_m2_stats(df_pw_v2, metric_filter, ab_short)
    
    annot_v1 <- get_m2_sig_label(df_pw_v1, metric_filter, ab_short)
    annot_v2 <- get_m2_sig_label(df_pw_v2, metric_filter, ab_short)
    
    # Interaction stats
    interaction_str <- if (show_interaction) get_m3b_interaction(df_anova, metric_filter) else ""
    
    # Logic for V1
    p_adj_v1 <- stats_v1$p_adj
    p_boot_v1 <- stats_v1$p_boot
    if (!is.na(p_adj_v1) && p_adj_v1 <= 0.05) {
        sig_v1 <- if (annot_v1 != "") paste0("", annot_v1) else ""
        str_v1 <- sprintf("p_adj = %s%s [boot p = %s]", fmt_p_val(p_adj_v1), sig_v1, fmt_p_val(p_boot_v1))
    } else {
        str_v1 <- sprintf("p_adj = %s", fmt_p_val(p_adj_v1))
    }
    
    # Logic for V2
    p_adj_v2 <- stats_v2$p_adj
    p_boot_v2 <- stats_v2$p_boot
    if (!is.na(p_adj_v2) && p_adj_v2 <= 0.05) {
        sig_v2 <- if (annot_v2 != "") paste0("", annot_v2) else ""
        str_v2 <- sprintf("p_adj = %s%s [boot p = %s]", fmt_p_val(p_adj_v2), sig_v2, fmt_p_val(p_boot_v2))
    } else {
        str_v2 <- sprintf("p_adj = %s", fmt_p_val(p_adj_v2))
    }
    
    if (show_interaction && interaction_str != "") {
        txt <- sprintf("A/PR/8/34 Bys vs Inf: %s\nA/NWS/33 Bys vs Inf: %s\n%s",
                       str_v1, str_v2, interaction_str)
    } else {
        txt <- sprintf("A/PR/8/34 Bys vs Inf: %s\nA/NWS/33 Bys vs Inf: %s",
                       str_v1, str_v2)
    }
    return(txt)
}

# Helper function to plot standard bar panels for A1-A3 & B1
plot_bar_panel <- function(df_cells, df_pw_v1, df_pw_v2, df_anova, metric_filter, raw_col, ab_short, title_text, ylab, show_interaction = TRUE, show_stats = TRUE, show_group_labels = TRUE) {
    # Filter cells (Naive is excluded)
    df_sub <- df_cells %>%
        filter(antibody == ab_short, cond_state != "Naive") %>%
        mutate(cond_state = factor(cond_state, levels = c("A/PR/8/34.Bystander", "A/PR/8/34.Infected", "A/NWS/33.Bystander", "A/NWS/33.Infected")))
    
    # Calculate summary stats (Mean, SEM)
    df_sum <- df_sub %>%
        group_by(cond_state) %>%
        summarize(
            M = mean(.data[[raw_col]], na.rm = TRUE),
            SE = sd(.data[[raw_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[raw_col]]))),
            .groups = "drop"
        )
    
    # Subtitle with boot p and interaction F-stat integration
    subtitle_txt <- if (show_stats) {
      get_bar_subtitle(df_pw_v1, df_pw_v2, df_anova, metric_filter, ab_short, show_interaction)
    } else {
      NULL
    }
    
    # Simple effects labels
    annot_v1 <- get_m2_sig_label(df_pw_v1, metric_filter, ab_short)
    annot_v2 <- get_m2_sig_label(df_pw_v2, metric_filter, ab_short)
    
    comparisons_list <- list(
        c("A/PR/8/34.Bystander", "A/PR/8/34.Infected"),
        c("A/NWS/33.Bystander", "A/NWS/33.Infected")
    )
    annotations_list <- c(annot_v1, annot_v2)
    
    # Outlier management
    max_errorbar <- max(df_sum$M + df_sum$SE, na.rm = TRUE)
    raw_max <- as.numeric(max(df_sub[[raw_col]], na.rm = TRUE))
    
    data_top <- max(raw_max, max_errorbar)
    bracket_y <- data_top * 1.15
    y_max <- data_top * 1.40
    if (raw_col == "n_nuclear_puncta") {
        y_max <- 400
        bracket_y <- 320
    }
    
    p <- ggplot(df_sub, aes(x = cond_state, y = .data[[raw_col]])) +
        geom_col(data = df_sum, aes(x = cond_state, y = M, fill = cond_state),
                 alpha = 0.6, width = 0.7, color = NA, key_glyph = draw_key_split, show.legend = TRUE) +
        geom_errorbar(data = df_sum, aes(x = cond_state, y = M,
                                         ymin = pmax(0, M - SE), ymax = M + SE,
                                         color = cond_state),
                      width = 0.22, linewidth = 0.6, show.legend = FALSE) +
        geom_jitter(data = df_sub %>% filter(.data[[raw_col]] <= y_max),
                    aes(color = cond_state),
                    shape = 21, fill = "white", size = 1.5, stroke = 0.5,
                    alpha = 0.4, width = 0.15, height = 0, show.legend = FALSE) +
        scale_fill_manual(values = cond_state_colors,
                          breaks = c("A/PR/8/34.Infected", "A/NWS/33.Infected"),
                          labels = c("A/PR/8/34", "A/NWS/33"),
                          name = NULL, drop = FALSE) +
        scale_color_manual(values = cond_state_colors, guide = "none") +
        scale_x_discrete(labels = c(
            "A/PR/8/34.Bystander" = "Bystander",
            "A/PR/8/34.Infected"  = "Infected",
            "A/NWS/33.Bystander"  = "Bystander",
            "A/NWS/33.Infected"   = "Infected"
        ), drop = FALSE) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, y_max), oob = scales::oob_keep) +
        labs(x = NULL, y = ylab, title = title_text, subtitle = subtitle_txt) +
        theme_classic(base_size = 12) +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(size = 9, face = "bold"),
            axis.title.y = element_text(face = "bold", size = 11),
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)),
            plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0, lineheight = 1.2, margin = margin(b = 22))
        )
    
    # Add brackets if significant
    keep_idx <- which(annotations_list != "ns" & annotations_list != "" & !is.na(annotations_list))
    if (length(keep_idx) > 0) {
        levels_vector <- c("A/PR/8/34.Bystander", "A/PR/8/34.Infected", "A/NWS/33.Bystander", "A/NWS/33.Infected")
        xmin_coords <- sapply(comparisons_list[keep_idx], function(x) match(x[1], levels_vector))
        xmax_coords <- sapply(comparisons_list[keep_idx], function(x) match(x[2], levels_vector))
        
        raw_range <- diff(range(df_sub[[raw_col]], na.rm = TRUE))
        visible_range <- y_max
        rel_tip <- if (!is.na(raw_range) && raw_range > 0) (0.02 * visible_range) / raw_range else 0.02
        
        p <- p + ggsignif::geom_signif(
            xmin = xmin_coords,
            xmax = xmax_coords,
            annotations = annotations_list[keep_idx],
            y_position = rep(bracket_y, length(keep_idx)),
            tip_length = rel_tip, textsize = 3.5, color = "black"
        )
    }
    
    if (show_group_labels) {
      p <- add_grouped_x_axis(p, 0, y_max, is_five_groups = FALSE)
    }
    return(p)
}


message("Generating Panel A...")
A1 <- plot_bar_panel(cells, m2v1_pw, m2v2_pw, m3b_anova, "asyn_nuc [Syn1]", "asyn_nuclear_mean", "Syn1",
                     "Nuclear α-Syn", "Intensity (A.U.)", show_stats = show_cd_stats, show_group_labels = FALSE)

A2 <- plot_bar_panel(cells, m2v1_pw, m2v2_pw, m3b_anova, "asyn_cyto [Syn1]", "asyn_cytoplasmic_mean", "Syn1",
                     "Perinuclear α-Syn", "Intensity (A.U.)", show_stats = show_cd_stats, show_group_labels = FALSE)

A3 <- plot_bar_panel(cells, m2v1_pw, m2v2_pw, m3b_anova, "n_puncta", "n_nuclear_puncta", "Syn1",
                     "Nuclear α-Syn puncta", "Count per cell", show_stats = show_cd_stats, show_group_labels = FALSE)


message("Generating Panel B...")
# B1: NC_ratio Syn1 bar plot
B1 <- plot_bar_panel(cells, m2v1_pw, m2v2_pw, m3b_anova, "NC_ratio", "NC_ratio", "Syn1",
                     "Nuclear α-Syn enrichment", "N/P ratio", show_interaction = FALSE, show_stats = show_cd_stats, show_group_labels = FALSE)

# B2: NC_ratio Syn1 Ridge Plot (with Naive Group added as baseline reference)
B2_df <- cells %>%
    filter(antibody == "Syn1") %>%
    mutate(
        cond_state = factor(cond_state, levels = c("Naive", "A/PR/8/34.Bystander", "A/PR/8/34.Infected", "A/NWS/33.Bystander", "A/NWS/33.Infected"))
    ) %>%
    filter(!is.na(cond_state))

B2 <- ggplot(B2_df, aes(x = NC_ratio, y = fct_rev(cond_state), fill = cond_state)) +
    geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01,
                        color = "black", linewidth = 0.4, show.legend = FALSE) +
    scale_fill_manual(values = cond_state_colors, drop = FALSE) +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_discrete(labels = c(
        "Naive"               = "Naive",
        "A/PR/8/34.Bystander" = "A/PR/8/34 Bystander",
        "A/PR/8/34.Infected"  = "A/PR/8/34 Infected",
        "A/NWS/33.Bystander"  = "A/NWS/33 Bystander",
        "A/NWS/33.Infected"   = "A/NWS/33 Infected"
    )) +
    labs(x = "N/P ratio (log₁₀)", y = NULL,
         title = "Nuclear α-Syn enrichment",
         subtitle = NULL) +
    theme_classic(base_size = 12) +
    theme(
        axis.text.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)),
        plot.subtitle = element_text(size = 8.5, color = "grey30", hjust = 0.5, margin = margin(b = 15))
    )

# B3: violin chromatin colocalization (coloc_dapi_r, Syn1)
plot_b3_violin <- function(df_cells, df_pw_v1, df_pw_v2, df_anova, show_stats = TRUE, show_group_labels = TRUE) {
    metric_filter <- "coloc_dapi_r"
    ab_short <- "Syn1"
    
    df_sub <- df_cells %>%
        filter(antibody == ab_short, cond_state != "Naive", !is.na(nuc_asyn_dapi_pearson_r)) %>%
        mutate(cond_state = factor(cond_state, levels = c("A/PR/8/34.Bystander", "A/PR/8/34.Infected", "A/NWS/33.Bystander", "A/NWS/33.Infected")))
    
    y_min_data <- as.numeric(quantile(df_sub$nuc_asyn_dapi_pearson_r, 0.02, na.rm = TRUE))
    y_max_data <- as.numeric(quantile(df_sub$nuc_asyn_dapi_pearson_r, 0.98, na.rm = TRUE))
    
    df_sub_filtered <- df_sub %>% 
        filter(nuc_asyn_dapi_pearson_r >= y_min_data, nuc_asyn_dapi_pearson_r <= y_max_data)
    
    y_span <- y_max_data - y_min_data
    y_min <- y_min_data - y_span * 0.30
    bracket_y <- y_max_data + y_span * 0.45
    y_max <- bracket_y + y_span * 0.50
    
    stats_v1 <- get_m2_stats(df_pw_v1, metric_filter, ab_short)
    stats_v2 <- get_m2_stats(df_pw_v2, metric_filter, ab_short)
    
    annot_v1 <- get_m2_sig_label(df_pw_v1, metric_filter, ab_short)
    annot_v2 <- get_m2_sig_label(df_pw_v2, metric_filter, ab_short)
    
    # Logic for V1
    p_adj_v1 <- stats_v1$p_adj
    p_boot_v1 <- stats_v1$p_boot
    if (!is.na(p_adj_v1) && p_adj_v1 <= 0.05) {
        sig_v1 <- if (annot_v1 != "") paste0("", annot_v1) else ""
        str_v1 <- sprintf("p_adj = %s%s [boot p = %s]", fmt_p_val(p_adj_v1), sig_v1, fmt_p_val(p_boot_v1))
    } else {
        str_v1 <- sprintf("p_adj = %s", fmt_p_val(p_adj_v1))
    }
    
    # Logic for V2
    p_adj_v2 <- stats_v2$p_adj
    p_boot_v2 <- stats_v2$p_boot
    if (!is.na(p_adj_v2) && p_adj_v2 <= 0.05) {
        sig_v2 <- if (annot_v2 != "") paste0("", annot_v2) else ""
        str_v2 <- sprintf("p_adj = %s%s [boot p = %s]", fmt_p_val(p_adj_v2), sig_v2, fmt_p_val(p_boot_v2))
    } else {
        str_v2 <- sprintf("p_adj = %s", fmt_p_val(p_adj_v2))
    }
    
    subtitle_txt <- if (show_stats) {
      sprintf("A/PR/8/34 Bys vs Inf: %s\nA/NWS/33 Bys vs Inf: %s", str_v1, str_v2)
    } else {
      NULL
    }
    
    p <- ggplot(df_sub_filtered, aes(x = cond_state, y = nuc_asyn_dapi_pearson_r)) +
        geom_violin(aes(fill = cond_state, color = cond_state), trim = FALSE, alpha = 0.3,
                    scale = "width", width = 0.75, linewidth = 0.4, show.legend = FALSE) +
        geom_boxplot(aes(color = cond_state), width = 0.12, fill = NA,
                     outlier.shape = NA, linewidth = 0.25, show.legend = FALSE) +
        geom_jitter(aes(color = cond_state), shape = 21, fill = "white", stroke = 0.5,
                    alpha = 0.5, width = 0.15, height = 0, size = 1.5, show.legend = FALSE) +
        scale_fill_manual(values = cond_state_colors, guide = "none", drop = FALSE) +
        scale_color_manual(values = cond_state_colors, guide = "none", drop = FALSE) +
        scale_x_discrete(labels = c(
            "A/PR/8/34.Bystander" = "Bystander",
            "A/PR/8/34.Infected"  = "Infected",
            "A/NWS/33.Bystander"  = "Bystander",
            "A/NWS/33.Infected"   = "Infected"
        ), drop = FALSE) +
        scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max), oob = scales::oob_keep) +
        labs(x = NULL, y = "Pearson r",
             title = "α-Syn & DAPI colocalization",
             subtitle = subtitle_txt) +
        theme_classic(base_size = 12) +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(size = 9, face = "bold"),
            axis.title.y = element_text(face = "bold", size = 11),
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)),
            plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0, lineheight = 1.2, margin = margin(b = 22))
        )
    
    comps_5c2 <- list(
        c("A/PR/8/34.Bystander", "A/PR/8/34.Infected"),
        c("A/NWS/33.Bystander", "A/NWS/33.Infected")
    )
    annots_5c2 <- c(annot_v1, annot_v2)
    keep_5c2 <- which(annots_5c2 != "ns" & annots_5c2 != "" & !is.na(annots_5c2))
    
    if (length(keep_5c2) > 0) {
        levels_vector <- c("A/PR/8/34.Bystander", "A/PR/8/34.Infected", "A/NWS/33.Bystander", "A/NWS/33.Infected")
        xmin_coords <- sapply(comps_5c2[keep_5c2], function(x) match(x[1], levels_vector))
        xmax_coords <- sapply(comps_5c2[keep_5c2], function(x) match(x[2], levels_vector))
        
        raw_range <- diff(range(df_sub_filtered$nuc_asyn_dapi_pearson_r, na.rm = TRUE))
        visible_range <- y_max - y_min
        rel_tip <- if (!is.na(raw_range) && raw_range > 0) (0.02 * visible_range) / raw_range else 0.02
        
        p <- p + ggsignif::geom_signif(
            xmin = xmin_coords,
            xmax = xmax_coords,
            annotations = annots_5c2[keep_5c2],
            y_position = rep(bracket_y, length(keep_5c2)),
            tip_length = rel_tip, textsize = 2.6, color = "black"
        )
    }
    
    if (show_group_labels) {
      p <- add_grouped_x_axis(p, y_min, y_max, is_five_groups = FALSE)
    }
    return(p)
}

B3 <- plot_b3_violin(cells, m2v1_pw, m2v2_pw, m3b_anova, show_stats = show_cd_stats, show_group_labels = FALSE)

# ==============================================================================
# PANEL C: Dose-response (Infected cells, condition V1+V2, iav_mean > 0)
# ==============================================================================
message("Generating Panel C...")

# Sub-helper to compute partial Pearson r dynamically inside script
compute_partial_r_stats <- function(sub_df, raw_col) {
  n_cells <- nrow(sub_df)
  n_fovs <- if (n_cells > 0) n_distinct(sub_df$image_id) else 0
  
  if (n_cells < 15) {
    return(list(r = NA_real_, p = NA_real_, note = ""))
  }
  
  x <- log10(sub_df$iav_mean)
  y <- sub_df[[raw_col]]
  
  if (n_fovs <= 1) {
    test_res <- cor.test(x, y, method = "pearson")
    return(list(r = as.numeric(test_res$estimate), p = as.numeric(test_res$p.value), note = ""))
  } else if (n_cells - n_fovs - 1 <= 0) {
    return(list(r = NA_real_, p = NA_real_, note = ""))
  } else {
    res_y <- residuals(lm(y ~ image_id, data = sub_df))
    res_x <- residuals(lm(x ~ image_id, data = sub_df))
    part_r <- cor(res_x, res_y)
    
    r_clamped <- max(min(part_r, 0.9999), -0.9999)
    z <- atanh(r_clamped)
    se <- 1 / sqrt(n_cells - n_fovs - 1)
    p_val <- 2 * (1 - pnorm(abs(z) / se))
    return(list(r = part_r, p = p_val, note = ""))
  }
}

# Helper to format r stats string
get_r_string <- function(r_stats) {
  if (is.na(r_stats$r)) return("NA")
  
  p <- r_stats$p
  stars <- if (is.na(p)) ""
           else if (p < 0.001) "***"
           else if (p < 0.01)  "**"
           else if (p < 0.05)  "*"
           else                ""
  
  return(sprintf("r=%.2f%s%s", r_stats$r, stars, r_stats$note))
}

get_r_parts <- function(r_stats) {
  if (is.na(r_stats$r)) {
    return(list(r_text = "r=NA", stars = ""))
  }

  p <- r_stats$p
  stars <- if (is.na(p)) ""
           else if (p < 0.001) "***"
           else if (p < 0.01)  "**"
           else if (p < 0.05)  "*"
           else                ""

  list(r_text = sprintf("r=%.2f", r_stats$r), stars = stars)
}

# Helper to retrieve slopes and slope diffs from Excel sheets & combine with partial r
get_m4_subtitle_compact <- function(cells_df, df_m4a_v1, df_m4a_v2, df_m4b, metric_val, raw_col, ab_filter) {
  # V1 cells
  df_v1_cells <- cells_df %>%
    filter(condition == "A/PR/8/34", IAV_state == "Infected", !is.na(iav_mean), iav_mean > 0)
  if (!is.null(ab_filter)) df_v1_cells <- df_v1_cells %>% filter(antibody == ab_filter)
  r_stats_v1 <- compute_partial_r_stats(df_v1_cells, raw_col)
  r_str_v1 <- get_r_string(r_stats_v1)
  
  # V2 cells
  df_v2_cells <- cells_df %>%
    filter(condition == "A/NWS/33", IAV_state == "Infected", !is.na(iav_mean), iav_mean > 0)
  if (!is.null(ab_filter)) df_v2_cells <- df_v2_cells %>% filter(antibody == ab_filter)
  r_stats_v2 <- compute_partial_r_stats(df_v2_cells, raw_col)
  r_str_v2 <- get_r_string(r_stats_v2)
  
  # Fetch slope stats
  row_v1 <- df_m4a_v1 %>% filter(Metric == metric_val)
  slope_v1   <- if(nrow(row_v1) > 0) row_v1$Estimate[1] else NA_real_
  p_v1       <- if(nrow(row_v1) > 0) row_v1$`Pr(>|t|)`[1] else NA_real_
  
  row_v2 <- df_m4a_v2 %>% filter(Metric == metric_val)
  slope_v2   <- if(nrow(row_v2) > 0) row_v2$Estimate[1] else NA_real_
  p_v2       <- if(nrow(row_v2) > 0) row_v2$`Pr(>|t|)`[1] else NA_real_
  
  row_diff <- df_m4b %>% filter(Metric == metric_val)
  p_diff      <- if(nrow(row_diff) > 0) row_diff$P_adjusted[1] else NA_real_
  sig_diff    <- if(nrow(row_diff) > 0) row_diff$Significance[1] else ""
  
  fmt_val <- function(v, digits = 3) {
    if (is.na(v)) return("NA")
    return(sprintf(paste0("%.", digits, "f"), v))
  }
  
  sig_diff_clean <- if (is.na(sig_diff) || sig_diff == "ns" || sig_diff == "") "" else paste0(" [", sig_diff, "]")
  
  txt <- sprintf("PR8 b=%s (p=%s), %s | NWS b=%s (p=%s), %s | Delta p=%s%s",
                 fmt_val(slope_v1), fmt_p_val(p_v1), r_str_v1,
                 fmt_val(slope_v2), fmt_p_val(p_v2), r_str_v2,
                 fmt_p_val(p_diff), sig_diff_clean)
  return(txt)
}

plot_c_scatter <- function(df_cells, df_m4a_v1, df_m4a_v2, df_m4b, metric_val, raw_col, ab_filter, title_text, ylab, show_stats = FALSE) {
    df_sub <- df_cells %>%
        filter(condition %in% c("A/PR/8/34", "A/NWS/33"), IAV_state == "Infected", !is.na(iav_mean), iav_mean > 0)
    
    if (!is.null(ab_filter)) {
        df_sub <- df_sub %>% filter(antibody == ab_filter)
    }
    
    df_sub <- df_sub %>%
        mutate(
            log_iav = log10(iav_mean),
            condition = factor(condition, levels = c("A/PR/8/34", "A/NWS/33"))
        ) %>%
        filter(!is.na(log_iav), !is.na(.data[[raw_col]]))
    
    y_max_c <- as.numeric(quantile(df_sub[[raw_col]], 0.95, na.rm = TRUE))
    
    subtitle_txt <- if (show_stats) {
      get_m4_subtitle_compact(df_cells, df_m4a_v1, df_m4a_v2, df_m4b, metric_val, raw_col, ab_filter)
    } else {
      NULL
    }

    x_rng <- range(df_sub$log_iav, na.rm = TRUE)
    x_span <- max(0.1, diff(x_rng))
    cond_levels <- c("A/PR/8/34", "A/NWS/33")
    ann_rows <- lapply(seq_along(cond_levels), function(i) {
      cond_i <- cond_levels[i]
      sub_i <- df_sub %>% filter(condition == cond_i)
      parts <- get_r_parts(compute_partial_r_stats(sub_i, raw_col))
      tibble(
        condition = cond_i,
        x = x_rng[1] + 0.04 * x_span,
        y = Inf,
        vjust_row = 1.05 + (i - 1) * 1.35,
        r_label = paste0(parts$r_text, ifelse(parts$stars == "", "", paste0(" ", parts$stars)))
      )
    })
    ann_df <- bind_rows(ann_rows)
    
    p <- ggplot(df_sub, aes(x = log_iav, y = .data[[raw_col]], color = condition)) +
        geom_point(data = df_sub %>% filter(.data[[raw_col]] <= y_max_c), alpha = 0.2, size = 1.5, aes(fill = condition), shape = 21, stroke = 0.3) +
        geom_smooth(data = df_sub %>% filter(.data[[raw_col]] <= y_max_c), method = "lm", formula = y ~ x, se = TRUE, aes(fill = condition), alpha = 0.15, linewidth = 0.8) +
        geom_label(data = ann_df, aes(x = x, y = y, label = r_label, color = condition, vjust = vjust_row), hjust = 0, size = 2.0, fontface = "bold", fill = scales::alpha("white", 0.78), linewidth = 0, label.padding = unit(0.04, "lines"), show.legend = FALSE) +
        scale_color_manual(values = condition_colors, name = NULL) +
        scale_fill_manual(values = condition_colors, name = NULL) +
        scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
        labs(
            x = expression(bold(log[10](IAV))), y = ylab,
            title = title_text,
            subtitle = subtitle_txt
        ) +
        theme_classic(base_size = 12) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(t = 5, b = 5)),
            plot.subtitle = element_text(size = 5.4, color = "grey30", hjust = 0, lineheight = 1.0, margin = margin(b = 6)),
            axis.title = element_text(face = "bold", size = 11),
            plot.margin = margin(t = 10, r = 10, b = 25, l = 10, unit = "pt")
        ) +
        coord_cartesian(ylim = c(0, y_max_c), clip = "off")
    
    return(p)
}

plot_c_scatter_pair <- function(df_cells, title_text, ab_filter, metric_nuc, col_nuc, metric_peri, col_peri, show_stats = FALSE) {
  build_one <- function(metric_val, raw_col, compartment_label) {
    df_sub <- df_cells %>%
      filter(condition %in% c("A/PR/8/34", "A/NWS/33"), IAV_state == "Infected", antibody == ab_filter, !is.na(iav_mean), iav_mean > 0) %>%
      mutate(
        log_iav = log10(iav_mean),
        condition = factor(condition, levels = c("A/PR/8/34", "A/NWS/33")),
        y = .data[[raw_col]],
        compartment = compartment_label
      ) %>%
      filter(!is.na(log_iav), !is.na(y))

    y_cap <- as.numeric(quantile(df_sub$y, 0.95, na.rm = TRUE))
    df_sub %>% filter(y <= y_cap) %>% select(log_iav, condition, y, compartment)
  }

  df_plot <- bind_rows(
    build_one(metric_nuc, col_nuc, "Nuclear"),
    build_one(metric_peri, col_peri, "Perinuclear")
  )

  subtitle_txt <- if (show_stats) {
    txt_nuc <- get_m4_subtitle_compact(df_cells, m4a_v1, m4a_v2, m4b_diff, metric_nuc, col_nuc, ab_filter)
    txt_peri <- get_m4_subtitle_compact(df_cells, m4a_v1, m4a_v2, m4b_diff, metric_peri, col_peri, ab_filter)
    paste0("Nuclear: ", txt_nuc, "\nPerinuclear: ", txt_peri)
  } else {
    NULL
  }

  build_ann <- function(metric_val, raw_col, compartment_label) {
    df_sub <- df_cells %>%
      filter(condition %in% c("A/PR/8/34", "A/NWS/33"), IAV_state == "Infected", antibody == ab_filter, !is.na(iav_mean), iav_mean > 0) %>%
      mutate(log_iav = log10(iav_mean), y = .data[[raw_col]], condition = factor(condition, levels = c("A/PR/8/34", "A/NWS/33"))) %>%
      filter(!is.na(log_iav), !is.na(y))

    # r uses all infected cells; the 95th-percentile cap is for display only
    y_cap <- as.numeric(quantile(df_sub$y, 0.95, na.rm = TRUE))
    x_rng <- range(df_sub$log_iav[df_sub$y <= y_cap], na.rm = TRUE)
    x_span <- max(0.1, diff(x_rng))

    cond_levels <- c("A/PR/8/34", "A/NWS/33")
    ann_rows <- lapply(seq_along(cond_levels), function(i) {
      cond_i <- cond_levels[i]
      sub_i <- df_sub %>% filter(condition == cond_i)
      parts <- get_r_parts(compute_partial_r_stats(sub_i, "y"))
      tibble(
        compartment = compartment_label,
        condition = cond_i,
        x = x_rng[1] + 0.04 * x_span,
        y = Inf,
        vjust_row = 1.05 + (i - 1) * 1.35,
        r_label = paste0(parts$r_text, ifelse(parts$stars == "", "", paste0(" ", parts$stars)))
      )
    })
    bind_rows(ann_rows)
  }

  ann_df <- bind_rows(
    build_ann(metric_nuc, col_nuc, "Nuclear"),
    build_ann(metric_peri, col_peri, "Perinuclear")
  )

  ggplot(df_plot, aes(x = log_iav, y = y, color = condition)) +
    geom_point(alpha = 0.2, size = 1.5, aes(fill = condition), shape = 21, stroke = 0.25) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(fill = condition), alpha = 0.15, linewidth = 0.7) +
    geom_label(data = ann_df, aes(x = x, y = y, label = r_label, color = condition, vjust = vjust_row), hjust = 0, size = 1.8, fontface = "bold", fill = scales::alpha("white", 0.78), linewidth = 0, label.padding = unit(0.04, "lines"), inherit.aes = FALSE, show.legend = FALSE) +
    facet_wrap(~compartment, nrow = 1) +
    scale_color_manual(values = condition_colors, name = "Virus strain") +
    scale_fill_manual(values = condition_colors, name = "Virus strain") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    labs(
      x = expression(bold(log[10](IAV))),
      y = "Intensity (A.U.)",
      title = title_text,
      subtitle = subtitle_txt
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(t = 5, b = 5)),
      plot.subtitle = element_text(size = 5.2, color = "grey30", hjust = 0, lineheight = 1.0, margin = margin(b = 5)),
      axis.title = element_text(face = "bold", size = 10),
      plot.margin = margin(t = 10, r = 10, b = 25, l = 10, unit = "pt")
    )
}

C1 <- plot_c_scatter(cells, m4a_v1, m4a_v2, m4b_diff, "NC_ratio", "NC_ratio", "Syn1",
                     "N/P & IAV correlation", "N/P ratio", show_stats = show_scatter_stats)

E1 <- plot_c_scatter_pair(cells,
                          title_text = "α-Syn & IAV correlation",
                          ab_filter = "Syn1",
                          metric_nuc = "asyn_nuc [Syn1]",
                          col_nuc = "asyn_nuclear_mean",
                          metric_peri = "asyn_cyto [Syn1]",
                          col_peri = "asyn_cytoplasmic_mean",
                          show_stats = show_scatter_stats)

E2 <- plot_c_scatter_pair(cells,
                          title_text = "<i>a-Syn</i><sub style='font-size:6pt;'>EPR20535</sub> & IAV correlation",
                          ab_filter = "EPR20535",
                          metric_nuc = "asyn_nuc [EPR20535]",
                          col_nuc = "asyn_nuclear_mean",
                          metric_peri = "asyn_cyto [EPR20535]",
                          col_peri = "asyn_cytoplasmic_mean",
                          show_stats = show_scatter_stats)

# ==============================================================================

# ==============================================================================
# SHRINK THEME & LOAD IMAGES
# ==============================================================================
shrink_theme <- theme(
  axis.text.y = element_text(size = 4),
  axis.text.x = element_text(size = 6.5, face = "bold", angle = 30, hjust = 1),
  axis.title.x = element_text(size = 6.2, face = "bold", margin = margin(t = 2)),
  axis.title.y = element_text(size = 6.2, face = "bold", margin = margin(r = 2)),
  plot.title = ggtext::element_textbox_simple(
    size = 9.0,
    face = "bold",
    halign = 0.5,
    width = unit(1, "npc"),
    fill = "grey92",
    box.color = "grey75",
    r = unit(3, "pt"),
    padding = margin(1.5, 2.5, 1.5, 2.5),
    margin = margin(t = 1, b = 3)
  ),
  plot.subtitle = element_text(size = 3.6, color = "grey30", margin = margin(b = 1)),
  plot.title.position = "plot",
  legend.position = "none", 
  plot.margin = margin(5, 5, 24, 5)
)
A1 <- A1 + shrink_theme
A2 <- A2 + shrink_theme
A3 <- A3 + shrink_theme
B1 <- B1 + shrink_theme
B2 <- B2 + shrink_theme
B3 <- B3 + shrink_theme
C1 <- C1 + shrink_theme
E1 <- E1 + shrink_theme
E2 <- E2 + shrink_theme
C1 <- C1 + theme(plot.margin = margin(5, 5, 8, 5))
E1 <- E1 + theme(plot.margin = margin(5, 5, 8, 5))
E2 <- E2 + theme(plot.margin = margin(5, 5, 8, 5))
E2 <- E2 + theme(
  plot.title = ggtext::element_textbox_simple(
    size = 8.2,
    face = "bold",
    halign = 0.5,
    width = unit(1, "npc"),
    fill = "grey92",
    box.color = "grey75",
    r = unit(3, "pt"),
    padding = margin(1.5, 2.5, 1.5, 2.5),
    margin = margin(t = 1, b = 3)
  )
)
# Panels a-b: written by Figures/Figure-05_Representative-Images.py
png_dir <- file.path(out_dir_final, "Figure5_Representative")
imgA <- readPNG(file.path(png_dir, "PanelA_combined_orthogonal.png"))
imgB <- readPNG(file.path(png_dir, "PanelB_combined_zoom.png"))
pA_img <- ggdraw() + draw_image(imgA, scale = 1.00) + theme(plot.margin = margin(4, 8, 3, 8))
pB_img <- ggdraw() + draw_image(imgB, scale = 1.03) + theme(plot.margin = margin(4, 8, 3, 15.7))
assign_tag <- function(p, tag_text, y_pos = 1, x_pos = -0.015) {
  p + labs(tag = tag_text) +
    theme(
      plot.tag = element_text(size = 18, face = "bold", hjust = 0, vjust = 1),
      plot.tag.position = c(x_pos, y_pos)
    )
}
pA_img <- assign_tag(pA_img, "a", y_pos = 1.00, x_pos = -0.06)
pB_img <- assign_tag(pB_img, "b", y_pos = 1.00, x_pos = -0.06)
A1 <- assign_tag(A1, "c", y_pos = 1.03, x_pos = -0.206)
B1 <- assign_tag(B1, "d", y_pos = 1.03, x_pos = -0.206)
E1 <- assign_tag(E1, "e1", y_pos = 1.03, x_pos = -0.206)
E2 <- assign_tag(E2, "e2", y_pos = 1.03, x_pos = -0.206)
C1 <- assign_tag(C1, "e3", y_pos = 1.03, x_pos = -0.206)
A1 <- A1 + theme(plot.margin = margin(5, 5, 10, 5))
A2 <- A2 + theme(plot.margin = margin(5, 5, 10, 5))
A3 <- A3 + theme(plot.margin = margin(5, 5, 10, 5))
B1 <- B1 + theme(plot.margin = margin(5, 5, 10, 5), plot.title = ggtext::element_textbox_simple(size=9.0, face="bold", halign=0.5, width=unit(1,"npc"), fill="grey92", box.color="grey75", r=unit(3,"pt"), padding=margin(1.5, 2.5, 1.5, 2.5), margin=margin(t=1, b=8)))
B2 <- B2 + theme(plot.margin = margin(5, 5, 10, 5), axis.text.y = element_text(size = 6.0), plot.title = ggtext::element_textbox_simple(size=9.0, face="bold", halign=0.5, width=unit(1,"npc"), fill="grey92", box.color="grey75", r=unit(3,"pt"), padding=margin(1.5, 2.5, 1.5, 2.5), margin=margin(t=1, b=8)))
B3 <- B3 + theme(plot.margin = margin(5, 5, 10, 5), plot.title = ggtext::element_textbox_simple(size=9.0, face="bold", halign=0.5, width=unit(1,"npc"), fill="grey92", box.color="grey75", r=unit(3,"pt"), padding=margin(1.5, 2.5, 1.5, 2.5), margin=margin(t=1, b=8)))
row3_core <- plot_grid(A1, NULL, A2, NULL, A3, ncol = 5, rel_widths = c(1, 0.22, 1, 0.22, 1), align = "h", axis = "b")
row4_core <- plot_grid(B1, NULL, B2, NULL, B3, ncol = 5, rel_widths = c(1, 0.22, 1, 0.22, 1), align = "h", axis = "b")
row5_core <- plot_grid(E1, NULL, E2, NULL, C1, ncol = 5, rel_widths = c(1, 0.22, 1, 0.22, 1), align = "h", axis = "b")

row3 <- plot_grid(NULL, row3_core, NULL, ncol = 3, rel_widths = c(0.040, 0.86, 0.100))
row4 <- plot_grid(NULL, row4_core, NULL, ncol = 3, rel_widths = c(0.040, 0.86, 0.100))
row5 <- plot_grid(NULL, row5_core, NULL, ncol = 3, rel_widths = c(0.040, 0.86, 0.100))

legend_items <- tibble(
  virus = c("A/PR/8/34", "A/NWS/33"),
  x0 = c(2.2, 4.9),
  light = c("#A8C7F0", "#F2B4B4"),
  dark = c("#4285F4", "#C0392B")
)

tri_upper <- legend_items %>%
  transmute(
    grp = paste0(virus, "_u"), fill = light,
    x1 = x0, y1 = 0.82,
    x2 = x0, y2 = 1.18,
    x3 = x0 + 0.34, y3 = 1.18
  ) %>%
  pivot_longer(cols = c(x1, y1, x2, y2, x3, y3), names_to = c("axis", "pt"), names_pattern = "([xy])(.)") %>%
  pivot_wider(names_from = axis, values_from = value) %>%
  select(grp, fill, x, y)

tri_lower <- legend_items %>%
  transmute(
    grp = paste0(virus, "_l"), fill = dark,
    x1 = x0, y1 = 0.82,
    x2 = x0 + 0.34, y2 = 0.82,
    x3 = x0 + 0.34, y3 = 1.18
  ) %>%
  pivot_longer(cols = c(x1, y1, x2, y2, x3, y3), names_to = c("axis", "pt"), names_pattern = "([xy])(.)") %>%
  pivot_wider(names_from = axis, values_from = value) %>%
  select(grp, fill, x, y)

tri_df <- bind_rows(tri_upper, tri_lower)

band_rect <- legend_items %>%
  transmute(xmin = x0 + 0.50, xmax = x0 + 0.90, ymin = 0.90, ymax = 1.10, fill = scales::alpha(dark, 0.20))

line_seg <- legend_items %>%
  transmute(x = x0 + 0.56, xend = x0 + 0.88, y = 1.00, yend = 1.00, col = dark)

label_df <- legend_items %>%
  transmute(x = x0 + 1.06, y = 1.00, label = virus)

row6_core <- ggplot() +
  geom_polygon(data = tri_df, aes(x = x, y = y, group = grp, fill = fill), color = NA) +
  geom_rect(data = band_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), color = NA) +
  geom_segment(data = line_seg, aes(x = x, xend = xend, y = y, yend = yend, color = col), linewidth = 0.8) +
  geom_text(data = label_df, aes(x = x, y = y, label = label), hjust = 0, vjust = 0.5, size = 2.2) +
  scale_fill_identity() +
  scale_color_identity() +
  coord_cartesian(xlim = c(1.6, 7.4), ylim = c(0.76, 1.24), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

row6 <- plot_grid(NULL, row6_core, NULL, ncol = 3, rel_widths = c(0.07, 0.86, 0.07))

pA_row <- plot_grid(NULL, pA_img, NULL, ncol = 3, rel_widths = c(0.040, 0.86, 0.100))
pB_row <- plot_grid(NULL, pB_img, NULL, ncol = 3, rel_widths = c(0.040, 0.86, 0.100))
final_fig <- plot_grid(
  pA_row,
  pB_row,
  row3,
  row4,
  row5,
  row6,
  ncol = 1,
  rel_heights = c(11.2, 8.0, 4.35, 4.35, 4.5, 0.42)
)
pdf_path <- file.path(out_dir_final, "Figure5.pdf")
png_path <- file.path(out_dir_final, "Figure5.png")
ggsave(pdf_path, final_fig, width = 8.27, height = 11.69, device = cairo_pdf)
ggsave(png_path, final_fig, width = 8.27, height = 11.69, dpi = 1000, bg = "white")
message("Done! Saved to: ", pdf_path)
