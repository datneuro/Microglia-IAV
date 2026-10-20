

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggsignif)
  library(openxlsx)
  library(grid)
})

# ==============================================================================
# §1  PATHS CONFIGURATION
# ==============================================================================
# Output directory for the assembled figure.
OUTDIR <- Sys.getenv("FIGURE_OUT", unset = file.path("Figures", "output"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

QPCR_BASE   <- Sys.getenv("RTQPCR_OUT", unset = file.path("2606_PCR_VirusHMC3", "analysis_master"))
QPCR_CSV    <- file.path(QPCR_BASE, "04_per_rep_fc_sdm.csv")
QPCR_XLSX   <- file.path(QPCR_BASE, "04_master_qpcr_stats_report.xlsx")

PROT_BASE   <- Sys.getenv("PROTEASOME_OUT", unset = file.path("Proteasome", "processed"))
PROT_SA_CSV <- file.path(PROT_BASE, "Proteasome_SpecificActivity.csv")
PROT_RT_CSV <- file.path(PROT_BASE, "Proteasome_Ratio.csv")
PROT_XLSX   <- file.path(PROT_BASE, "Proteasome_Statistics_Log.xlsx")

STRAIN_V1 <- "A/PR/8/34"
STRAIN_V2 <- "A/NWS/33"

GROUP_MAP_QPCR <- tribble(
  ~Group_ID, ~Timepoint, ~LPS,      ~Virus,  ~ShortLabel,                    ~Label,                         ~in_factorial,
  "G1",  "4h",  "No LPS",  "None", "Ctrl",                            "4h-Ctrl",                            TRUE,
  "G2",  "4h",  "No LPS",  "V1",   STRAIN_V1,                         paste0("4h-", STRAIN_V1),             TRUE,
  "G3",  "4h",  "No LPS",  "V2",   STRAIN_V2,                         paste0("4h-", STRAIN_V2),             TRUE,
  "G4",  "4h",  "LPS",     "None", "LPS",                             "4h-LPS",                             TRUE,
  "G5",  "4h",  "LPS",     "V1",   paste0("LPS+", STRAIN_V1),         paste0("4h-LPS+", STRAIN_V1),         TRUE,
  "G6",  "4h",  "LPS",     "V2",   paste0("LPS+", STRAIN_V2),         paste0("4h-LPS+", STRAIN_V2),         TRUE,
  "G7",  "4h",  NA,        NA,     "PolyI:C",                         "4h-PolyI:C",                         FALSE,
  "G8",  "8h",  "No LPS",  "None", "Ctrl",                            "8h-Ctrl",                            TRUE,
  "G9",  "8h",  "No LPS",  "V1",   STRAIN_V1,                         paste0("8h-", STRAIN_V1),             TRUE,
  "G10", "8h",  "No LPS",  "V2",   STRAIN_V2,                         paste0("8h-", STRAIN_V2),             TRUE,
  "G11", "8h",  "LPS",     "None", "LPS",                             "8h-LPS",                             TRUE,
  "G12", "8h",  "LPS",     "V1",   paste0("LPS+", STRAIN_V1),         paste0("8h-LPS+", STRAIN_V1),         TRUE,
  "G13", "8h",  "LPS",     "V2",   paste0("LPS+", STRAIN_V2),         paste0("8h-LPS+", STRAIN_V2),         TRUE,
  "G14", "8h",  NA,        NA,     "PolyI:C",                         "8h-PolyI:C",                         FALSE
)

GRP_LEVELS_QPCR <- paste0("G", 1:14)

GROUP_COLORS_QPCR <- c(
  "G1"  = "#D5D8DC",  "G8"  = "#D5D8DC",   # light gray  — Ctrl
  "G2"  = "#A8C7F0",  "G9"  = "#A8C7F0",   # light blue  — V1
  "G3"  = "#F2B4B4",  "G10" = "#F2B4B4",   # light red   — V2
  "G4"  = "#6E7376",  "G11" = "#6E7376",   # dark gray   — LPS
  "G5"  = "#4285F4",  "G12" = "#4285F4",   # dark blue   — LPS+V1
  "G6"  = "#C0392B",  "G13" = "#C0392B",   # dark red    — LPS+V2
  "G7"  = "#2ECC71",  "G14" = "#2ECC71"    # stable green — PolyI:C
)

LEGEND_BREAKS_QPCR <- c("G1",  "G2",  "G3",  "G7")
LEGEND_LABELS_QPCR <- c("Naive", STRAIN_V1, STRAIN_V2, "Poly(I:C)")

draw_key_split_grp <- function(data, params, size) {
  fc <- data$fill
  pair <- if (is.na(fc))                               c("grey90",  "grey60") else
          if (fc %in% c("#D5D8DC", "#6E7376"))         c("#D5D8DC", "#6E7376") else
          if (fc %in% c("#A8C7F0", "#4285F4"))         c("#A8C7F0", "#4285F4") else
          if (fc %in% c("#F2B4B4", "#C0392B"))         c("#F2B4B4", "#C0392B") else
                                                       c(fc, fc)
  grid::grobTree(
    grid::polygonGrob(x = c(0, 0, 1), y = c(0, 1, 1),
                      gp = grid::gpar(fill = pair[1], col = NA)),
    grid::polygonGrob(x = c(0, 1, 1), y = c(0, 0, 1),
                      gp = grid::gpar(fill = pair[2], col = NA)),
    grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey50", lwd = 0.5))
  )
}

# Load qPCR Data
raw_qpcr <- read_csv(QPCR_CSV, show_col_types = FALSE)
if ("Gene" %in% colnames(raw_qpcr) && !"Target_Name" %in% colnames(raw_qpcr))
  raw_qpcr <- rename(raw_qpcr, Target_Name = Gene)
if ("Method" %in% colnames(raw_qpcr))
  raw_qpcr <- filter(raw_qpcr, grepl("SDM", Method))
if (!"Timepoint" %in% colnames(raw_qpcr))
  raw_qpcr <- mutate(raw_qpcr, Timepoint = if_else(Group_ID %in% paste0("G", 1:7), "4h", "8h"))

df_clean_qpcr <- raw_qpcr %>%
  filter(!is.na(FoldChange), is.finite(FoldChange), FoldChange > 0,
         !Target_Name %in% c("GAPDH", "RPS11")) %>%
  left_join(
    GROUP_MAP_QPCR %>% select(Group_ID, LPS, Virus, ShortLabel, Label, in_factorial),
    by = "Group_ID"
  ) %>%
  mutate(
    Group_ID  = factor(Group_ID,  levels = GRP_LEVELS_QPCR),
    Timepoint = factor(Timepoint, levels = c("4h", "8h")),
    Bio_Rep   = factor(Bio_Rep)
  )

contrasts_df <- openxlsx::read.xlsx(QPCR_XLSX, sheet = "2_All_Contrasts")
contrasts_df$P_adjusted_Holm <- as.numeric(contrasts_df$P_adjusted_Holm)
grp_order_qpcr <- setNames(seq_along(GRP_LEVELS_QPCR), GRP_LEVELS_QPCR)

signif_annot_base <- contrasts_df %>%
  filter(!is.na(P_adjusted_Holm), P_adjusted_Holm < 0.05, !Significance %in% c("ns", "")) %>%
  select(Target_Name = Gene, Timepoint, x = Group_ref, xend = Group_test, label = Significance) %>%
  mutate(span = grp_order_qpcr[xend] - grp_order_qpcr[x]) %>%
  group_by(Target_Name, Timepoint) %>%
  arrange(abs(span), .by_group = TRUE) %>%
  mutate(bracket_rank = row_number()) %>%
  ungroup() %>%
  select(-span)

# Per-gene original floor and ceiling from 04_plot_only_qpcr.R
Y_GENE_FLOOR_QPCR <- c(
  "IL1B"  = -6.0,
  "IL6"   = -3.0,
  "IRF7"  = -4.0,
  "MX1"   = -4.0,
  "OAS1"  = -4.0,
  "PSMB5" = -4.0,
  "PSMB8" = -4.0,
  "RSAD2" = -4.0
)
Y_GENE_CEIL_QPCR <- c(
  "IL1B"  =  6.0,
  "IL6"   =  5.5,
  "IRF7"  = 12.0,
  "MX1"   = 12.0,
  "OAS1"  = 12.0,
  "PSMB5" =  5.0,
  "PSMB8" =  9.0,
  "RSAD2" = 25.0
)

build_gene_plot <- function(gene) {
  jitter_g <- df_clean_qpcr %>%
    filter(Target_Name == gene, !is.na(FoldChange), FoldChange > 0) %>%
    mutate(.y_val = log2(FoldChange))

  if (nrow(jitter_g) == 0) return(ggplot() + labs(title = gene))

  y_vals    <- jitter_g$.y_val
  y_raw_min <- min(y_vals, na.rm = TRUE)
  y_raw_max <- max(y_vals, na.rm = TRUE)
  y_data_span <- max(y_raw_max - y_raw_min, 0.5)

  floor_val <- { tmp <- Y_GENE_FLOOR_QPCR[gene]; if (!is.na(tmp)) tmp else -0.5 }
  ceil_val  <- { tmp <- Y_GENE_CEIL_QPCR[gene];  if (!is.na(tmp)) tmp else  0.5 }
  y_bot0    <- min(y_raw_min - y_data_span * 0.28, floor_val)
  y_top0    <- max(y_raw_max + y_data_span * 0.28, ceil_val)

  # Enforce minimum total range of 5 log2 units like original
  if (y_top0 - y_bot0 < 5.0) {
    y_mid  <- (y_top0 + y_bot0) / 2
    y_bot0 <- y_mid - 2.5
    y_top0 <- y_mid + 2.5
  }

  sig_base <- filter(signif_annot_base, Target_Name == gene)
  if (nrow(sig_base) > 0) {
    # brk_base starts safely above highest data point and violin top tail
    brk_base  <- y_raw_max + max(y_data_span * 0.22, 0.80)
    n_brk_max <- max(sig_base$bracket_rank)
    # Distribute brackets across available space up to y_top0; minimum step 0.65
    brk_step  <- max((y_top0 - brk_base) / (n_brk_max + 1), 0.65)
    annot_df  <- sig_base %>%
      group_by(Timepoint) %>%
      mutate(y_pos = brk_base + brk_step * bracket_rank) %>%
      ungroup()
    y_top <- max(y_top0, max(annot_df$y_pos) + brk_step * 1.25)
  } else {
    annot_df <- sig_base
    y_top    <- y_top0
  }
  y_bot <- y_bot0

  y_plot_span <- y_top - y_bot
  bar_y <- y_bot - 0.055 * y_plot_span
  txt_y <- y_bot - 0.120 * y_plot_span

  suppressWarnings({
    ggplot(jitter_g, aes(x = Group_ID, y = .y_val, fill = Group_ID, color = Group_ID)) +
      annotate("rect", xmin = 0.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
               fill = "#F2F4F8", alpha = 0.40) +
      annotate("rect", xmin = 7.5, xmax = 14.5, ymin = -Inf, ymax = Inf,
               fill = "#FDF2E9", alpha = 0.40) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.35) +
      geom_vline(xintercept = 7.5, linetype = "dotted", color = "grey60", linewidth = 0.4) +
      # Layer 1: Smooth violin with adjust = 2.5, trim = FALSE (visible rounded body even for close dots!)
      geom_violin(aes(fill = Group_ID), alpha = 0.40, color = NA,
                  scale = "width", width = 0.95, trim = FALSE, adjust = 2.5) +
      # Layer 2: Compact, narrow boxplot inside violin (width = 0.16)
      geom_boxplot(aes(fill = Group_ID, color = Group_ID),
                   alpha = 0.50, width = 0.16, outlier.shape = NA, linewidth = 0.35) +
      # Layer 3: Jitter replicates inside boxplot
      geom_jitter(shape = 21, size = 1.15, stroke = 0.35, fill = "white",
                  alpha = 0.85, width = 0.05, height = 0) +
      annotate("text", x = 4.0,  y = y_top, label = "\u2190 4h \u2192",
               size = 2.4, color = "grey20", vjust = 1.05) +
      annotate("text", x = 11.0, y = y_top, label = "\u2190 8h \u2192",
               size = 2.4, color = "grey20", vjust = 1.05) +
      {
        if (nrow(annot_df) > 0)
          ggsignif::geom_signif(
            data = annot_df,
            aes(xmin = x, xmax = xend, annotations = label, y_position = y_pos),
            manual      = TRUE,
            tip_length  = 0.010,
            vjust       = 0.40,
            size        = 0.28,
            textsize    = 2.8,
            color       = "grey30",
            inherit.aes = FALSE
          )
        else list()
      } +
      # 4h +LPS bar & text (BELOW x-axis line at y_bot)
      annotate("segment", x = 3.6,  xend = 6.4,  y = bar_y, yend = bar_y,
               color = "#1A1A1A", linewidth = 0.8) +
      annotate("text",   x = 4.75, y = txt_y, label = "+ LPS",
               size = 1.65, fontface = "bold", color = "#1A1A1A", vjust = 1.0) +
      # 8h +LPS bar & text (BELOW x-axis line at y_bot)
      annotate("segment", x = 10.6, xend = 13.4, y = bar_y, yend = bar_y,
               color = "#1A1A1A", linewidth = 0.8) +
      annotate("text",   x = 11.75, y = txt_y, label = "+ LPS",
               size = 1.65, fontface = "bold", color = "#1A1A1A", vjust = 1.0) +
      # 4h Poly(I:C) bar & single-line text (BELOW x-axis line at y_bot)
      annotate("segment", x = 6.6,  xend = 7.4,  y = bar_y, yend = bar_y,
               color = "#2ECC71", linewidth = 0.8) +
      annotate("text",   x = 7.20, y = txt_y, label = "Poly(I:C)",
               size = 1.55, color = "#2ECC71", fontface = "bold", vjust = 1.0) +
      # 8h Poly(I:C) bar & single-line text (BELOW x-axis line at y_bot)
      annotate("segment", x = 13.6, xend = 14.4, y = bar_y, yend = bar_y,
               color = "#2ECC71", linewidth = 0.8) +
      annotate("text",   x = 14.20, y = txt_y, label = "Poly(I:C)",
               size = 1.55, color = "#2ECC71", fontface = "bold", vjust = 1.0) +
      scale_fill_manual(values = GROUP_COLORS_QPCR, guide = "none") +
      scale_color_manual(values = GROUP_COLORS_QPCR, guide = "none") +
      scale_x_discrete(labels = NULL, expand = expansion(add = 0.85)) +
      scale_y_continuous(expand = expansion(mult = c(0.0, 0.0))) +
      # Axis is strictly anchored at y_bot, with clip="off" to render annotations below axis
      coord_cartesian(ylim = c(y_bot, y_top), clip = "off") +
      labs(
        x = NULL,
        y = bquote(bolditalic(.(gene)) ~ bold(mRNA ~ (log[2]*FC)))
      ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        axis.text.y     = element_text(size = 6.5),
        axis.title.y    = element_text(size = 7.5, face = "bold"),
        # Generous top and bottom margins for outside annotations and inter-plot gaps
        plot.margin     = margin(t = 10, r = 6, b = 26, l = 6, unit = "pt"),
        legend.position = "none"
      )
  })
}

# Horizontal qPCR split-key legend (no title)
legend_df_qpcr <- tibble(
  x        = factor(LEGEND_BREAKS_QPCR, levels = GRP_LEVELS_QPCR),
  y        = 1,
  Group_ID = factor(LEGEND_BREAKS_QPCR, levels = GRP_LEVELS_QPCR)
)

p_legend_src_qpcr <- ggplot(legend_df_qpcr, aes(x = x, y = y, fill = Group_ID)) +
  geom_tile(width = 0.7, height = 0.5, key_glyph = draw_key_split_grp) +
  scale_fill_manual(
    values = GROUP_COLORS_QPCR,
    breaks = LEGEND_BREAKS_QPCR,
    labels = LEGEND_LABELS_QPCR,
    name   = NULL
  ) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8.0),
    legend.key.size   = unit(0.45, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width  = unit(0.50, "cm"),
    legend.spacing.x  = unit(0.35, "cm"),
    legend.key        = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin       = margin(t = 2, r = 0, b = 4, l = 0, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title = NULL))

p_qpcr_legend <- cowplot::get_legend(p_legend_src_qpcr)

target_genes <- c("IL1B", "IL6", "IRF7", "RSAD2", "MX1", "OAS1", "PSMB5", "PSMB8")
plot_list_qpcr <- lapply(target_genes, build_gene_plot)
names(plot_list_qpcr) <- target_genes

# Separate Row 1 & Row 2 with generous gap between rows (Point 1)
row_qpcr_1 <- cowplot::plot_grid(plotlist = plot_list_qpcr[1:4], ncol = 4, align = "hv")
row_qpcr_2 <- cowplot::plot_grid(plotlist = plot_list_qpcr[5:8], ncol = 4, align = "hv")

# Balanced layout: maximum plot height, crisp clean gaps, no excess empty space
panel_a_content <- cowplot::plot_grid(
  NULL,          # Top headroom for panel label "a"
  row_qpcr_1,
  NULL,          # Crisp gap between qPCR row 1 and row 2 (0.07)
  row_qpcr_2,
  NULL,          # Small gap before legend (0.03)
  p_qpcr_legend,
  ncol = 1,
  rel_heights = c(0.05, 0.95, 0.07, 0.95, 0.03, 0.10)
)

# Label "a" placed in dedicated headroom
panel_a_labeled <- cowplot::plot_grid(
  panel_a_content,
  labels = "a",
  label_size = 16,
  label_fontface = "bold",
  label_x = 0.01,
  label_y = 1.0
)

group_order_prot <- c(
  "Naive", "A/PR/8/34", "A/NWS/33",
  "LPS", "LPS+A/PR/8/34", "LPS+A/NWS/33"
)

palette_grbr <- c(
  "Naive"          = "#D5D8DC", # Light gray
  "A/PR/8/34"      = "#A8C7F0", # Pale blue
  "A/NWS/33"       = "#F2B4B4", # Pale red
  "LPS"            = "#6E7376", # Dark gray
  "LPS+A/PR/8/34"  = "#4285F4", # Bright blue
  "LPS+A/NWS/33"   = "#C0392B"  # Deep red
)

raw_sa <- read.csv(PROT_SA_CSV) %>%
  mutate(Group = factor(Group, levels = group_order_prot))
raw_rt <- read.csv(PROT_RT_CSV) %>%
  mutate(Group = factor(Group, levels = group_order_prot))

pc_prot <- openxlsx::read.xlsx(PROT_XLSX, sheet = "3_Planned_Contrasts")

get_sig_prot <- function(metric_name) {
  sig <- pc_prot %>%
    filter(Metric == metric_name, P_adjusted < 0.05)
  if (nrow(sig) == 0) return(list(pairs = list(), labels = character()))
  list(
    pairs = lapply(seq_len(nrow(sig)), function(i) c(sig$Group1[i], sig$Group2[i])),
    labels = sig$Significance
  )
}

sig_anw   <- get_sig_prot("Ac-ANW-AMC")
sig_llvy  <- get_sig_prot("Suc-LLVY-AMC")
sig_ratio <- get_sig_prot("Ratio_ANW_LLVY (log)")

# Exact original titles (title_text = NULL), exact original y-axis,
# and stratification indicators strictly OUTSIDE / BELOW the x-axis line at y = 0
create_bar_original <- function(df, ycol, ylab, sig_obj, title_text = NULL) {
  sdf <- df %>%
    group_by(Group) %>%
    summarise(
      m = mean(.data[[ycol]], na.rm = TRUE),
      s = sd(.data[[ycol]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[ycol]]))),
      .groups = "drop"
    )
  
  mx <- max(df[[ycol]], na.rm = TRUE)
  ns <- length(sig_obj$pairs)
  y_top <- if (ns > 0) mx * (1.18 + 0.12 * ns) else mx * 1.15
  
  # Negative offsets: strictly BELOW y = 0 (outside the plot)
  bar_y <- -0.075 * y_top
  txt_y <- -0.150 * y_top
  
  p <- ggplot(df, aes(x = Group, y = .data[[ycol]])) +
    geom_col(data = sdf, aes(x = Group, y = m, fill = Group),
             alpha = 0.6, width = 0.68, color = NA) +
    geom_errorbar(data = sdf,
                  aes(x = Group, y = m,
                      ymin = pmax(m - s, 0), ymax = m + s,
                      color = Group),
                  width = 0.22, linewidth = 0.5) +
    geom_jitter(aes(color = Group),
                shape = 21, fill = "white", size = 2.0, stroke = 0.6,
                alpha = 0.85, width = 0.14) +
    # Stratification indicators outside x-axis (below 0)
    annotate("segment", x = 0.8, xend = 3.2, y = bar_y, yend = bar_y,
             color = "#1A1A1A", linewidth = 0.8) +
    annotate("text", x = 2.0, y = txt_y, label = "No LPS",
             size = 2.6, fontface = "bold", color = "#1A1A1A", vjust = 1.0) +
    annotate("segment", x = 3.8, xend = 6.2, y = bar_y, yend = bar_y,
             color = "#1A1A1A", linewidth = 0.8) +
    annotate("text", x = 5.0, y = txt_y, label = "+ LPS",
             size = 2.6, fontface = "bold", color = "#1A1A1A", vjust = 1.0) +
    scale_fill_manual(values = palette_grbr, drop = FALSE) +
    scale_color_manual(values = palette_grbr, drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.5), drop = FALSE) +
    # Axis anchored at y = 0, annotations drawn in bottom margin via clip="off"
    coord_cartesian(ylim = c(0, y_top), clip = "off") +
    labs(y = ylab, x = NULL, title = title_text) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.text.y     = element_text(size = 7.0),
      axis.title.y    = element_text(face = "bold", size = 8.2),
      plot.title      = element_text(face = "bold", size = 9.0, hjust = 0.5),
      # Ample bottom margin for outside stratification text and gap between panels
      plot.margin     = margin(t = 6, r = 10, b = 28, l = 10, unit = "pt"),
      legend.position = "none"
    )
  
  if (ns > 0) {
    p <- p + ggsignif::geom_signif(
      comparisons   = sig_obj$pairs,
      annotations   = sig_obj$labels,
      step_increase = 0.10,
      tip_length    = 0.015,
      vjust         = 0.35,
      textsize      = 2.8,
      color         = "grey20"
    )
  }
  p
}

df_anw  <- raw_sa %>% filter(Substrate == "Ac-ANW-AMC")
df_llvy <- raw_sa %>% filter(Substrate == "Suc-LLVY-AMC")

p_bar_anw   <- create_bar_original(df_anw, "Specific_Activity",
                                   ylab = "Ac-ANW-AMC (RFU/min)", sig_obj = sig_anw, title_text = NULL)
p_bar_llvy  <- create_bar_original(df_llvy, "Specific_Activity",
                                   ylab = "Suc-LLVY-AMC (RFU/min)", sig_obj = sig_llvy, title_text = NULL)
p_ratio_df  <- raw_rt %>% mutate(Ratio_ANW_LLVY = Ratio_ANW_LLVY)
p_bar_ratio <- create_bar_original(p_ratio_df, "Ratio_ANW_LLVY",
                                   ylab = "Ac-ANW-AMC / Suc-LLVY-AMC (ratio)", sig_obj = sig_ratio, title_text = NULL)

# Point 3: Diagonal split 2-color square key legend for Panel b (matching Panel a style)
legend_df_prot <- tibble(
  x     = factor(c("Naive", "A/PR/8/34", "A/NWS/33"), levels = c("Naive", "A/PR/8/34", "A/NWS/33")),
  y     = 1,
  Group = factor(c("Naive", "A/PR/8/34", "A/NWS/33"), levels = c("Naive", "A/PR/8/34", "A/NWS/33"))
)

p_legend_src_prot <- ggplot(legend_df_prot, aes(x = x, y = y, fill = Group)) +
  geom_tile(width = 0.7, height = 0.5, key_glyph = draw_key_split_grp) +
  scale_fill_manual(
    values = c("Naive" = "#D5D8DC", "A/PR/8/34" = "#A8C7F0", "A/NWS/33" = "#F2B4B4"),
    breaks = c("Naive", "A/PR/8/34", "A/NWS/33"),
    labels = c("Naive", STRAIN_V1, STRAIN_V2),
    name   = NULL
  ) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8.0),
    legend.key.size   = unit(0.45, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width  = unit(0.50, "cm"),
    legend.spacing.x  = unit(0.35, "cm"),
    legend.key        = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin       = margin(t = 2, r = 0, b = 4, l = 0, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title = NULL))

shared_prot_legend <- cowplot::get_legend(p_legend_src_prot)

bars_row <- cowplot::plot_grid(
  p_bar_anw, p_bar_llvy, p_bar_ratio,
  ncol = 3, align = "h"
)

panel_b_content <- cowplot::plot_grid(
  NULL,          # Top headroom for panel label "b"
  bars_row,
  NULL,          # Vertical gap before legend
  shared_prot_legend,
  ncol = 1,
  rel_heights = c(0.06, 1.0, 0.04, 0.12)
)

panel_b_labeled <- cowplot::plot_grid(
  panel_b_content,
  labels = "b",
  label_size = 16,
  label_fontface = "bold",
  label_x = 0.01,
  label_y = 1.0
)

# ==============================================================================
# §4  PANEL C: PROTEASOME INTERACTION (ORIGINAL TITLES & Y-AXIS)
# ==============================================================================
fac_prot <- openxlsx::read.xlsx(PROT_XLSX, sheet = "4_Factorial_ANOVA")
mm_prot  <- openxlsx::read.xlsx(PROT_XLSX, sheet = "5_Marginal_Means")

make_factorial_annot_original <- function(metric_name) {
  f_tab <- fac_prot %>% filter(Metric == metric_name)
  if (nrow(f_tab) == 0) return("")
  lines <- vapply(seq_len(nrow(f_tab)), function(i) {
    term_lbl <- f_tab$Term[i]
    term_lbl <- gsub("LPS:Virus", "LPS \u00D7 Virus", term_lbl)
    p_val <- f_tab$P_Value[i]
    p_fmt <- if (is.na(p_val)) "NA" else if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
    eta_val <- f_tab$Partial_Eta2[i]
    eta_fmt <- if (is.na(eta_val)) "" else sprintf(", \u03B7\u00B2p=%.2f", eta_val)
    sig_star <- f_tab$Significance[i]
    star_fmt <- if (is.na(sig_star) || sig_star %in% c("ns", "")) "" else paste0(" ", sig_star)
    sprintf("%s: F(%g,12)=%.2f, p=%s%s%s",
            term_lbl, f_tab$df[i], f_tab$F_stat[i], p_fmt, eta_fmt, star_fmt)
  }, character(1))
  paste(lines, collapse = "\n")
}

get_int_df <- function(metric_name, on_log_scale = FALSE) {
  int_df <- mm_prot %>%
    filter(Metric == metric_name, Factor == "LPS:Virus") %>%
    separate(Level, into = c("LPS", "Virus"), sep = "_", remove = FALSE) %>%
    mutate(
      LPS = factor(LPS, levels = c("No", "Yes"), labels = c("No LPS", "+ LPS")),
      # Point 4: Change "No virus" to "Naive"
      Virus = factor(Virus, levels = c("None", "PR8", "NWS"),
                     labels = c("Naive", "A/PR/8/34", "A/NWS/33"))
    )
  if (on_log_scale) {
    int_df <- int_df %>%
      mutate(EMM_disp = exp(EMM),
             LCL_disp = exp(LCL),
             UCL_disp = exp(UCL))
  } else {
    int_df <- int_df %>%
      mutate(EMM_disp = EMM,
             LCL_disp = LCL,
             UCL_disp = UCL)
  }
  int_df %>% mutate(Virus_num = as.numeric(Virus))
}

lps_pair_colors <- c("No LPS" = "#6E7376", "+ LPS" = "#C0392B")

create_interaction_plot_original <- function(metric_name, on_log_scale, title_text, ylab) {
  int_df <- get_int_df(metric_name, on_log_scale)
  annot  <- make_factorial_annot_original(metric_name)
  
  y_max <- max(int_df$UCL_disp, na.rm = TRUE)
  y_min <- min(int_df$LCL_disp, na.rm = TRUE)
  y_range <- y_max - y_min
  y_annot    <- y_max + y_range * 0.44
  y_plot_top <- y_annot + y_range * 0.08
  y_plot_bot <- y_min - y_range * 0.08
  
  ggplot(int_df, aes(x = Virus_num, y = EMM_disp, color = LPS, fill = LPS, group = LPS)) +
    geom_ribbon(aes(ymin = LCL_disp, ymax = UCL_disp), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(aes(ymin = LCL_disp, ymax = UCL_disp), width = 0.08, linewidth = 0.4, alpha = 0.7) +
    geom_point(size = 2.4, shape = 16) +
    scale_x_continuous(
      breaks = unique(int_df$Virus_num),
      labels = levels(int_df$Virus),
      expand = expansion(add = 0.3)
    ) +
    scale_color_manual(values = lps_pair_colors, name = NULL) +
    scale_fill_manual(values = lps_pair_colors, name = NULL) +
    annotate("text", x = min(int_df$Virus_num) - 0.25, y = y_annot,
             label = annot, hjust = 0, vjust = 1,
             size = 1.80, family = "mono", lineheight = 1.06) +
    coord_cartesian(ylim = c(y_plot_bot, y_plot_top), clip = "off") +
    labs(title = title_text, x = NULL, y = ylab) +
    theme_classic(base_size = 9) +
    theme(
      axis.title.y    = element_text(face = "bold", size = 8.2),
      axis.text.x     = element_text(size = 7.0, face = "plain"),
      axis.text.y     = element_text(size = 7.0),
      plot.title      = element_text(face = "bold", size = 9.0, hjust = 0.5),
      legend.position = "none",
      plot.margin     = margin(t = 6, r = 10, b = 6, l = 10, unit = "pt")
    )
}

p_int_anw   <- create_interaction_plot_original(
  "Ac-ANW-AMC", FALSE,
  title_text = "Immunoproteasome activity",
  ylab       = "Ac-ANW-AMC (RFU/min)"
)
p_int_llvy  <- create_interaction_plot_original(
  "Suc-LLVY-AMC", FALSE,
  title_text = "Standard proteasome activity",
  ylab       = "Suc-LLVY-AMC (RFU/min)"
)
p_int_ratio <- create_interaction_plot_original(
  "Ratio_ANW_LLVY (log)", TRUE,
  title_text = "Immunoproteasome shift (ratio)",
  ylab       = "Ac-ANW-AMC / Suc-LLVY-AMC (ratio)"
)

shared_int_legend <- cowplot::get_legend(
  ggplot(data.frame(LPS = factor(c("No LPS", "+ LPS"), levels = c("No LPS", "+ LPS")), val = 1),
         aes(x = LPS, y = val, color = LPS)) +
    geom_point(size = 3.0) +
    geom_line(aes(group = 1), linewidth = 0.8) +
    scale_color_manual(values = lps_pair_colors, name = NULL) +
    theme_classic(base_size = 8.5) +
    theme(
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 8.0, face = "plain"),
      legend.key.size   = unit(0.40, "cm"),
      legend.spacing.x  = unit(0.20, "cm"),
      plot.margin       = margin(t = 2, r = 0, b = 4, l = 0, unit = "pt")
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL))
)

int_row <- cowplot::plot_grid(
  p_int_anw, p_int_llvy, p_int_ratio,
  ncol = 3, align = "h"
)

panel_c_content <- cowplot::plot_grid(
  NULL,          # Top headroom for panel label "c"
  int_row,
  NULL,          # Vertical gap before legend
  shared_int_legend,
  ncol = 1,
  rel_heights = c(0.06, 1.0, 0.04, 0.12)
)

panel_c_labeled <- cowplot::plot_grid(
  panel_c_content,
  labels = "c",
  label_size = 16,
  label_fontface = "bold",
  label_x = 0.01,
  label_y = 1.0
)

# ==============================================================================
# §5  ASSEMBLE FULL FIGURE 2 (A4 PORTRAIT: 8.27 × 11.69 in)
# ==============================================================================
# Generous vertical space for Panel a so violins and brackets have full room
final_fig2 <- cowplot::plot_grid(
  panel_a_labeled,
  NULL,              # Crisp gap between Panel a and b (0.05)
  panel_b_labeled,
  NULL,              # Crisp gap between Panel b and c (0.05)
  panel_c_labeled,
  ncol = 1,
  rel_heights = c(1.58, 0.05, 0.85, 0.05, 0.90)
)

# Add outer margins all around the entire figure to prevent any edge clipping
final_fig2_padded <- cowplot::plot_grid(
  NULL, final_fig2, NULL,
  ncol = 3,
  rel_widths = c(0.015, 0.97, 0.015)
)

# Target file paths
out_pdf      <- file.path(OUTDIR, "Figure2.pdf")
out_png_600  <- file.path(OUTDIR, "Figure2.png")
out_png_1200 <- file.path(OUTDIR, "Figure2_1200dpi.png")

# Export via a temporary file, then move into place (avoids file-lock issues on
# cloud-synced directories).
message("\n=== SAVING FIGURE 2 ARTIFACTS ===")

message("Saving vector PDF (cairo_pdf)...")
tmp_pdf <- tempfile(fileext = ".pdf")
ggsave(tmp_pdf, final_fig2_padded, width = 8.27, height = 11.69, device = cairo_pdf)
file.copy(tmp_pdf, out_pdf, overwrite = TRUE)
unlink(tmp_pdf)

message("Saving high-resolution PNG (600 dpi)...")
tmp_png_600 <- tempfile(fileext = ".png")
ggsave(tmp_png_600, final_fig2_padded, width = 8.27, height = 11.69, dpi = 600, bg = "white")
file.copy(tmp_png_600, out_png_600, overwrite = TRUE)
unlink(tmp_png_600)

message("Saving ultra high-resolution PNG (1200 dpi)...")
tmp_png_1200 <- tempfile(fileext = ".png")
ggsave(tmp_png_1200, final_fig2_padded, width = 8.27, height = 11.69, dpi = 1200, bg = "white")
file.copy(tmp_png_1200, out_png_1200, overwrite = TRUE)
unlink(tmp_png_1200)

message("\nSUCCESS: Figure2 created successfully at:\n  ",
        out_pdf, "\n  ", out_png_600, "\n  ", out_png_1200)
