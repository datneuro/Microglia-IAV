
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggrepel)
  library(magick)
  library(ggtext)
  library(Hmisc)
  library(shadowtext)
  library(forcats)
  library(grid)
})

# ---- Directories ------------------------------------------------------------
OUT1   <- Sys.getenv("FLOW_ELISA_OUT",      unset = file.path("Flow-ELISA", "R_analysis", "output"))
OUT2   <- Sys.getenv("BYSTANDER_OUT",       unset = file.path("Flow-ELISA", "R-analysis-figure_VS-bystander", "output"))
OUTDIR <- Sys.getenv("FIGURE_OUT",          unset = file.path("Figures", "output"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── 2. Shared Style & Palettes ──────────────────────────────────────────────
GROUP_COLORS <- c(C="#D5D8DC", C_PR8="#A8C7F0", C_S33="#F2B4B4",
                  LPS="#6E7376", LPS_PR8="#4285F4", LPS_S33="#C0392B",
                  PolyIC="#2ECC71")
VIRUS_COLORS <- c(None="#6E7376", V1="#4285F4", V2="#C0392B")
VIRUS_DISPLAY <- c(None="Control", V1="A/PR/8/34", V2="A/NWS/33")
STRAIN_V1 <- "A/PR/8/34"
STRAIN_V2 <- "A/NWS/33"
CONTRAST_LABS <- c("Ctrl-V1" = paste0(STRAIN_V1, " vs. Control"),
                   "Ctrl-V2" = paste0(STRAIN_V2, " vs. Control"),
                   "V1-V2"   = paste0(STRAIN_V1, " vs. ", STRAIN_V2))
CONTRAST_COLORS <- setNames(c(unname(VIRUS_COLORS["V1"]), unname(VIRUS_COLORS["V2"]), "#8E44AD"),
                            CONTRAST_LABS)
tp_fac <- function(x) factor(paste0(gsub("h$","",x), "h"), levels = c("4h", "8h", "24h"))
pfmt <- function(p) if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
HM_LOW <- "#114B9A"; HM_MID <- "#FFFFFF"; HM_HIGH <- "#C0392B"

theme_paper <- function(base = 9.5) {
  theme_classic(base_size = base) +
    theme(
      plot.title    = element_text(face = "bold", size = base, hjust = 0.5, margin = margin(t = 2, b = 2)),
      plot.subtitle = element_markdown(size = base - 3, color = "grey40", hjust = 0.5, lineheight = 1.1, margin = margin(b = 4)),
      axis.title    = element_text(face = "bold", size = base - 0.5),
      axis.text     = element_text(size = base - 2),
      axis.text.x   = element_markdown(size = base - 2, face = "bold"),
      strip.text    = element_textbox_simple(fill = "grey92", box.colour = NA, color = "black", halign = 0.5, valign = 0.5, r = unit(3, "pt"), padding = margin(2, 2, 2, 2), margin = margin(1, 1, 1, 1), face = "bold", size = base - 1.5),
      strip.background = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = base - 1.5),
      legend.text      = element_text(size = base - 2.5),
      legend.margin    = margin(t = 1, b = 0)
    )
}

# ── 3. Panel A
SCHEME_IMAGE_PATHS <- c(
  file.path(OUTDIR, "Fig1a_scheme.png"),
  file.path(OUTDIR, "20260929 Fig2 new.png"),   # file name of the published schematic
  file.path(OUTDIR, "Fig1a_scheme.pdf"),
  file.path(OUTDIR, "scheme.png"),
  file.path(OUTDIR, "schematic.png")
)
found_scheme <- SCHEME_IMAGE_PATHS[file.exists(SCHEME_IMAGE_PATHS)]

if (length(found_scheme) > 0) {
  message("Loading custom scheme image: ", found_scheme[1])
  scheme_img <- magick::image_read(found_scheme[1])
  # Auto-trim outer white borders so diagram fills space cleanly while strictly preserving aspect ratio
  scheme_trimmed <- magick::image_trim(scheme_img)
  plot_scheme <- cowplot::ggdraw() +
    cowplot::draw_image(scheme_trimmed, scale = 0.94, halign = 0, valign = 0.5)
} else {
  steps_df <- tibble(
    step = 1:4,
    x = c(1.3, 3.8, 6.2, 8.8),
    title = c(
      "1 · Microglia Model",
      "2 · Infection & Priming",
      "3 · Timecourse Dishes",
      "4 · Multi-omic Phenotyping"
    ),
    detail = c(
      "<b>HMC3 human microglia</b><br>Quiescent baseline<br>Matched replicates (n = 5)",
      "<b>IAV strains:</b> PR8 vs NWS<br>± LPS priming (100 ng/mL)<br>Differential viral replication",
      "<b>4 h, 8 h, 24 h post-infection</b><br>Parallel independent cultures<br>Dynamic activation tracking",
      "<b>Flow Cytometry (Surface + NP)</b><br><b>ELISA (Secreted Cytokines)</b><br>Cell-autonomous vs paracrine"
    ),
    fill_col = c("#EBF5FB", "#FDEDEC", "#FEF9E7", "#EAFAF1"),
    border_col = c("#5DADE2", "#EC7063", "#F5B041", "#58D68D")
  )
  arrows_df <- tibble(
    x = c(2.45, 4.95, 7.35),
    xend = c(2.65, 5.15, 7.55),
    y = rep(1.05, 3)
  )

  plot_scheme <- ggplot() +
    geom_rect(data = steps_df, aes(xmin = x - 1.15, xmax = x + 1.15, ymin = 0.22, ymax = 1.95),
              fill = steps_df$fill_col, color = steps_df$border_col, linewidth = 0.55) +
    geom_richtext(data = steps_df, aes(x = x, y = 1.74, label = title),
                  size = 3.2, fontface = "bold", fill = NA, label.color = NA, hjust = 0.5) +
    geom_richtext(data = steps_df, aes(x = x, y = 0.96, label = detail),
                  size = 2.45, lineheight = 1.18, fill = NA, label.color = NA, hjust = 0.5, color = "grey20") +
    geom_segment(data = arrows_df, aes(x = x, y = y, xend = xend, yend = y),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 color = "grey45", linewidth = 0.8) +
    coord_cartesian(xlim = c(0, 10.1), ylim = c(0.1, 2.05), expand = FALSE) +
    theme_void() +
    theme(
      plot.margin = margin(2, 6, 2, 6),
      plot.background = element_rect(fill = "white", color = "grey85", linewidth = 0.3)
    )
}

# ── 4. Panel B: Forest Plot (formerly Panel a) ──────────────────────────────
idx_levels <- c("Surface_M1", "Secreted_M1", "M2", "Decoupling_Index")
forest <- read_csv(file.path(OUT1, "L2_index_contrasts.csv"), show_col_types = FALSE) %>%
  filter(Contrast %in% c("Ctrl-V1", "Ctrl-V2", "V1-V2"), Index %in% idx_levels) %>%
  mutate(
    flip = if_else(Contrast %in% c("Ctrl-V1", "Ctrl-V2"), -1, 1),
    estimate = estimate * flip,
    new_lo = if_else(flip < 0, CI_high * flip, CI_low),
    new_hi = if_else(flip < 0, CI_low * flip, CI_high),
    CI_low = new_lo, CI_high = new_hi,
    Contrast = factor(CONTRAST_LABS[Contrast], levels = rev(CONTRAST_LABS)),
    Timepoint = tp_fac(Timepoint),
    Index = factor(Index, levels = rev(idx_levels)),
    star = if_else(p_holm < 0.05, "*", "")
  ) %>%
  mutate(y_label = case_when(
    grepl("Surface_M1", Index) ~ "<span style='font-size:9pt;'><b>Surface M1</b></span><br><span style='font-size:6.5pt; color:#555555;'><i>Virus Effect:</i> F = 60.94,<br>&eta;<sup>2</sup><sub>p</sub> = 0.70, <i>p</i> < 0.001</span>",
    grepl("Secreted_M1", Index) ~ "<span style='font-size:9pt;'><b>Secreted M1</b></span><br><span style='font-size:6.5pt; color:#555555;'><i>Time × Virus:</i> F = 18.71,<br>&eta;<sup>2</sup><sub>p</sub> = 0.55, <i>p</i> < 0.001</span>",
    grepl("M2", Index) ~ "<span style='font-size:9pt;'><b>M2</b></span><br><span style='font-size:6.5pt; color:#555555;'><i>Time Effect:</i> F = 1.13,<br>&eta;<sup>2</sup><sub>p</sub> = 0.05, <i>p</i> = 0.328</span>",
    grepl("Decoupling", Index) ~ "<span style='font-size:9pt;'><b>Decoupling Index</b></span><br><span style='font-size:6.5pt; color:#555555;'><i>Virus Effect:</i> F = 21.69,<br>&eta;<sup>2</sup><sub>p</sub> = 0.49, <i>p</i> < 0.001</span>"
  ))

dodge <- position_dodge(width = 0.6)
plot_forest <- ggplot(forest, aes(x = estimate, y = Index, color = Contrast)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = "grey88", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.45) +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high), position = dodge, size = 0.1, linewidth = 0.3) +
  geom_text(aes(x = CI_high, label = star), position = dodge, hjust = -0.4, vjust = 0.80, size = 3.2, show.legend = FALSE) +
  facet_wrap(~ Timepoint, nrow = 1) + scale_color_manual(values = CONTRAST_COLORS) +
  scale_y_discrete(labels = setNames(forest$y_label[match(rev(idx_levels), forest$Index)], rev(idx_levels))) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = "Standardized effect (estimate ± 95% CI)", y = NULL) + theme_paper() +
  theme(
    axis.text.y = ggtext::element_markdown(lineheight = 1.15, hjust = 1, margin = margin(r = 4)),
    axis.text.x = ggtext::element_markdown(size = 8.5, face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    panel.spacing.x = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.justification = "left",
    legend.key.size = unit(0.35, "lines"),
    legend.text = element_text(size = 6, face = "plain"),
    legend.margin = margin(t = 2, b = 0)
  ) +
  guides(color = guide_legend(title = NULL, nrow = 2, reverse = TRUE, override.aes = list(size = 0.5, linewidth = 0.5)))

# ── 5. Panel C: Trajectories (formerly Panel b) ─────────────────────────────
anova_idx <- read_csv(file.path(OUT1, "L2_index_anova.csv"), show_col_types = FALSE)
contr_idx <- read_csv(file.path(OUT1, "L2_index_contrasts.csv"), show_col_types = FALSE)
di_v <- anova_idx %>% filter(Index == "Decoupling_Index", Term == "Virus")
pi_lps <- anova_idx %>% filter(Index == "Polarization_Index", Term == "LPS")
di24 <- contr_idx %>% filter(Index == "Decoupling_Index", Timepoint == 24)
v1_delta <- (di24 %>% filter(Contrast == "Ctrl-V1"))$estimate[1]; v1_p <- (di24 %>% filter(Contrast == "Ctrl-V1"))$p_holm[1]
v2_delta <- (di24 %>% filter(Contrast == "Ctrl-V2"))$estimate[1]; v2_p <- (di24 %>% filter(Contrast == "Ctrl-V2"))$p_holm[1]

metric_html <- function(F, eta, p) sprintf("<i>F</i> = %.2f, <i>&eta;</i><sup>2</sup><sub><i>p</i></sub> = %.2f, <i>p</i> %s", F, eta, pfmt(p))
lab_di <- paste0("<b>Virus Effect</b> · ", metric_html(di_v$F_hc3, di_v$partial_eta2, di_v$p_hc3))
lab_pi <- paste0("<b>LPS Priming dominates</b> · <i>&eta;</i><sup>2</sup><sub><i>p</i></sub> = ", sprintf("%.2f", pi_lps$partial_eta2), ", <i>p</i> ", pfmt(pi_lps$p_hc3), "<br>Virus Effect: <i>p</i> > 0.05")

comp <- read_csv(file.path(OUT1, "L2_composites.csv"), show_col_types = FALSE) %>%
  mutate(Virus = factor(Virus, levels = c("None", "V1", "V2")), LPS = factor(LPS, levels = c("No", "Yes")),
         Timepoint_num = as.numeric(Timepoint), Timepoint_fac = tp_fac(Timepoint))
idx_lab <- c(Decoupling_Index = "Decoupling<br>Index (DI)", Polarization_Index = "Polarization<br>Index (PI)")
di_pi <- comp %>% select(Group, Timepoint_num, Virus, LPS, Bio_Rep, Decoupling_Index, Polarization_Index) %>%
  pivot_longer(c(Decoupling_Index, Polarization_Index), names_to = "Index", values_to = "Value") %>%
  mutate(Index = factor(Index, levels = c("Decoupling_Index", "Polarization_Index")))
ytop <- di_pi %>% group_by(Index) %>% summarise(ymax = max(Value), .groups = "drop")
anova_anno <- tibble(Index = factor(c("Decoupling_Index", "Polarization_Index"), levels = levels(di_pi$Index)),
                     label = c(lab_di, lab_pi)) %>% left_join(ytop, by = "Index")
end24 <- comp %>% filter(Timepoint_num == 24, Group %in% c("C_PR8", "C_S33")) %>% group_by(Group, Virus) %>% summarise(y = mean(Decoupling_Index), .groups = "drop")
tail_anno <- tibble(Index = factor("Decoupling_Index", levels = levels(di_pi$Index)), Virus = factor(c("V1", "V2"), levels = c("None", "V1", "V2")),
                    x = 24, y = c(end24$y[end24$Virus == "V1"], end24$y[end24$Virus == "V2"]), ny = c(0.65, -0.65),
                    label = c(sprintf("A/PR/8/34: Δ = %.2f, <i>p</i><sub>Holm</sub> %s", v1_delta, pfmt(v1_p)), sprintf("A/NWS/33: Δ = %.2f, <i>p</i><sub>Holm</sub> %s", v2_delta, pfmt(v2_p))))

set.seed(42)
plot_traj <- ggplot(di_pi, aes(Timepoint_num, Value, color = Virus, fill = Virus, linetype = LPS, group = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.45, height = 0), alpha = 0.16, size = 0.75, shape = 16, show.legend = FALSE) +
  stat_summary(data = filter(di_pi, LPS == "No"), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.13, color = NA, show.legend = FALSE) +
  stat_summary(aes(linewidth = LPS), fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point", shape = 21, fill = "white", size = 1.4, stroke = 0.55, show.legend = FALSE) +
  ggtext::geom_richtext(data = anova_anno, aes(x = 4, y = ymax, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.4, fill = alpha("white", 0.15), label.color = NA) +
  ggtext::geom_richtext(data = tail_anno, aes(x = x, y = y + ny, label = label, color = Virus), inherit.aes = FALSE, hjust = 1, size = 2.4, fill = alpha("white", 0.25), label.color = NA, show.legend = FALSE) +
  facet_grid(Index ~ ., scales = "free_y", labeller = labeller(Index = idx_lab), switch = "y") +
  scale_x_continuous(breaks = c(4, 8, 24), expand = expansion(mult = c(0.07, 0.05))) +
  scale_color_manual(values = VIRUS_COLORS, labels = VIRUS_DISPLAY, name = "Strain") +
  scale_fill_manual(values = VIRUS_COLORS, labels = VIRUS_DISPLAY, name = "Strain") +
  scale_linetype_manual(values = c(No = "solid", Yes = "22"), labels = c(No="No LPS", Yes="LPS Priming"), name = "Priming") +
  scale_linewidth_manual(values = c(No = 0.55, Yes = 0.2), guide = "none") +
  labs(x = "Time post-infection (h)", y = NULL) + theme_paper() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_textbox_simple(size = 9, face = "bold", orientation = "left-rotated", fill = "grey92", box.colour = NA, color = "black", r = unit(3, "pt"), padding = margin(2, 2, 2, 2), margin = margin(t=1, b=1, l=1, r=1), halign = 0.5, valign = 0.5),
    strip.background.y = element_blank(),
    panel.spacing.y = unit(0.35, "lines"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "pt"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(l = 8),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.text.x = ggtext::element_markdown(size = 8.5, face = "bold"),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.key.size = unit(0.35, "lines"),
    legend.text = element_text(size = 6, face = "plain"),
    legend.margin = margin(t = 2, b = 0)
  )

# ── 6. Panel D: Biplot (formerly Panel c) ────────────────────────────────────
CORE_GRP <- c("C", "C_PR8", "C_S33")
GROUP_DISPLAY <- c(C = "Control", C_PR8 = STRAIN_V1, C_S33 = STRAIN_V2)
TIME_ALPHA <- c("4h" = 0.10, "8h" = 0.22, "24h" = 0.36)
scores <- read_csv(file.path(OUT1, "L3_pca_scores.csv"), show_col_types = FALSE) %>%
  filter(LPS == "No") %>% mutate(Group = factor(Group, levels = CORE_GRP), Timepoint = tp_fac(Timepoint))
loadings <- read_csv(file.path(OUT1, "L3_pca_loadings.csv"), show_col_types = FALSE) %>%
  mutate(
    Platform = case_when(
      Marker %in% c("CD40","CD86","HLA-I","HLA-II","CD206") ~ "Surface (Flow)",
      Marker %in% c("TNF-a","IL-1b","IL-6","IL-10","IL-12p70","CCL2") ~ "Secreted (ELISA)",
      TRUE ~ "Other (Flow)"
    ),
    arrow_col = c("Surface (Flow)"="#4682B4", "Secreted (ELISA)"="#DC143C", "Other (Flow)"="#7F8C8D")[Platform],
    lab = recode(Marker, "IL-1b"="IL-1β", "TNF-a"="TNF-α")
  )
arrow_mult <- as.numeric(quantile(sqrt(scores$PC1^2 + scores$PC2^2), 0.90) / max(sqrt(loadings$PC1^2 + loadings$PC2^2)))
loadings <- mutate(loadings, x1 = PC1 * arrow_mult, y1 = PC2 * arrow_mult)

plot_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey75") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey75") +
  stat_ellipse(data = scores, aes(PC1, PC2, fill = Group, color = Group, alpha = Timepoint, group = interaction(Group, Timepoint)),
               geom = "polygon", level = 0.68, linewidth = 0.35, show.legend = c(color=FALSE, fill=FALSE, alpha=TRUE)) +
  geom_point(data = scores, aes(PC1, PC2, color = Group), size = 1.4, alpha = 0.75, shape = 16) +
  geom_segment(data = loadings, aes(x=0, y=0, xend=x1, yend=y1, color=I(arrow_col)), arrow = arrow(length = unit(0.09, "cm")), linewidth = 0.35, show.legend = FALSE) +
  geom_text_repel(data = loadings, aes(x=x1, y=y1, label=lab, color=I(arrow_col)), size=2.8, fontface="bold", bg.color="white", bg.r=0.15, segment.color="grey80", segment.size=0.2, box.padding=0.3, point.padding=0.12, force=1.5, max.overlaps=Inf, min.segment.length=0, seed=42, show.legend=FALSE) +
  scale_color_manual(values = GROUP_COLORS, breaks = CORE_GRP, labels = GROUP_DISPLAY[CORE_GRP], name = "Group") +
  scale_fill_manual(values = GROUP_COLORS, breaks = CORE_GRP, labels = GROUP_DISPLAY[CORE_GRP], name = "Group") +
  scale_alpha_manual(values = TIME_ALPHA, name = "Timepoint") +
  scale_shape_manual(values = c("4h" = 21, "8h" = 21, "24h" = 21)) +
  guides(color = guide_legend(order = 1, override.aes = list(size=1, alpha=1, shape=16, linetype=0)), fill = "none", alpha = guide_legend(order = 2, override.aes = list(fill="grey30", color="grey30"))) +
  labs(x = "PC1 — Activation Magnitude", y = "PC2 — Opposition Axis") + theme_paper() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    axis.text.x = ggtext::element_markdown(size = 8.5, face = "bold"),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 6, face = "plain"),
    legend.key.size = unit(0.35, "lines"),
    legend.margin = margin(t = 2, b = 0)
  )

# ── 7. Panel E: Radar (formerly Panel d) ─────────────────────────────────────
coord_radar <- function(theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggplot2::ggproto("CoordRadar", ggplot2::CoordPolar, theta = theta, r = r, start = start, direction = sign(direction), is_linear = function(coord) TRUE)
}
FLOW_ORDER <- c("CD40", "CD86", "HLA-I", "HLA-II", "CD206", "CD68")
ELISA_ORDER <- c("TNF-a", "IL-1b", "IL-6", "IL-12p70", "IL-10", "CCL2")
MK_LEVELS <- c(FLOW_ORDER, ELISA_ORDER)
greek <- function(m) recode(m, "TNF-a"="TNF-α", "IL-1b"="IL-1β", .default=m)
star_markers <- read_csv(file.path(OUT1, "L1_contrasts.csv"), show_col_types = FALSE) %>% filter(Contrast == "V1-V2", Timepoint == 24, p_holm < 0.05) %>% pull(Marker)
star_html <- "<span style='color:#D35400;font-size:13pt'><b>*</b></span>"
MK_DISPLAY <- setNames(paste0(greek(MK_LEVELS), ifelse(MK_LEVELS %in% star_markers, star_html, "")), MK_LEVELS)
ST_LEVELS <- c(paste("A/PR/8/34 (V1)", c("4h","8h","24h"), sep=" · "), paste("A/NWS/33 (V2)", c("4h","8h","24h"), sep=" · "))

rad <- read_csv(file.path(OUT1, "tidy_factorial.csv"), show_col_types=FALSE) %>%
  filter(Group %in% c("C_PR8", "C_S33"), Marker != "NP") %>%
  group_by(Marker) %>% mutate(scaled = (Value - min(Value)) / (max(Value) - min(Value))) %>% ungroup() %>%
  mutate(strain = factor(if_else(Group == "C_PR8", "A/PR/8/34 (V1)", "A/NWS/33 (V2)"), levels = c("A/PR/8/34 (V1)", "A/NWS/33 (V2)")),
         Virus = factor(Virus, levels = c("None", "V1", "V2")), tlab = paste0(Timepoint, "h"), Timepoint_fac = tp_fac(Timepoint),
         ST = factor(paste(strain, tlab, sep=" · "), levels = ST_LEVELS), Marker = factor(Marker, levels = MK_LEVELS),
         gcol = scales::alpha(VIRUS_COLORS[as.character(Virus)], 0.06)) %>% arrange(ST, Bio_Rep, Marker)
rad_mean <- rad %>% group_by(ST, Virus, Timepoint_fac, Marker) %>% summarise(scaled = mean(scaled), .groups="drop") %>% arrange(ST, Marker)
rings <- tidyr::expand_grid(level = c(0.25, 0.5, 0.75, 1.0), Marker = factor(MK_LEVELS, levels=MK_LEVELS))
ring_lab <- tibble(Marker = factor("CD40", levels=MK_LEVELS), level = c(0.25, 0.5, 0.75, 1.0), lab = sprintf("%.2f", c(0.25, 0.5, 0.75, 1.0)))

plot_radar <- ggplot() +
  geom_polygon(data=rings, aes(Marker, level, group=level), fill=NA, color="grey82", linewidth=0.25) +
  geom_polygon(data=rad, aes(Marker, scaled, group=interaction(ST, Bio_Rep), color=I(gcol)), fill=NA, linewidth=0.15, show.legend=FALSE) +
  geom_polygon(data=rad_mean, aes(Marker, scaled, group=ST, fill=Virus, alpha=Timepoint_fac), color=NA) +
  geom_polygon(data=rad_mean, aes(Marker, scaled, group=ST, color=Virus, linetype=Timepoint_fac, linewidth=Timepoint_fac), fill=NA, key_glyph="path") +
  shadowtext::geom_shadowtext(data=ring_lab, aes(Marker, level, label=lab), inherit.aes=FALSE, size=2.3, fontface="bold", color="black", bg.colour="white", bg.r=0.15) +
  coord_radar(theta="x", start=-75*pi/180) + scale_x_discrete(labels=MK_DISPLAY) + scale_y_continuous(limits=c(0,1), expand=expansion(mult=c(0, 0.05))) +
  scale_color_manual(values=VIRUS_COLORS, breaks=c("V1","V2"), labels=VIRUS_DISPLAY[c("V1","V2")], name="Strain") +
  scale_fill_manual(values=VIRUS_COLORS, guide="none") + scale_alpha_manual(values=c("4h"=0, "8h"=0, "24h"=0.08), guide="none") +
  scale_linetype_manual(values=c("4h"="dotted", "8h"="solid", "24h"="solid"), name="Timepoint") +
  scale_linewidth_manual(values=c("4h"=0.25, "8h"=0.2, "24h"=0.6), name="Timepoint") +
  labs(x=NULL, y=NULL) + theme_paper() +
  theme(
    panel.grid.major.x = element_line(color="grey90", linewidth=0.15),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = ggtext::element_markdown(size=13.5, face="bold"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "vertical",
    legend.spacing.y = unit(2, "pt"),
    legend.margin = margin(t=-2, b=0),
    legend.title = element_text(face="bold", size=7.5),
    legend.text = element_text(size=6.5, face="plain"),
    legend.key.size = unit(0.35, "lines"),
    legend.key.width = unit(1.1, "lines"),
    legend.key = element_rect(fill=NA, color=NA),
    legend.background = element_rect(fill=alpha("white",0.7), color=NA),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(linetype=guide_legend(order=1, override.aes=list(color="grey25")), linewidth=guide_legend(order=1), color=guide_legend(order=2, override.aes=list(linetype="solid", linewidth=0.6)))

# Safe Render Radar Plot (Trim whitespace & separate legend)
plot_radar_noleg <- plot_radar + theme(legend.position = "none")
temp_radar_raw <- tempfile(fileext = ".png")
ggsave(temp_radar_raw, patchwork::wrap_plots(plot_radar_noleg + theme(plot.margin = margin(16,16,16,16))), width = 5.5, height = 5.5, dpi = 1000, bg = "white")

img <- magick::image_read(temp_radar_raw)
img_cropped <- magick::image_trim(img)
temp_radar_cropped <- tempfile(fileext = ".png")
magick::image_write(img_cropped, temp_radar_cropped)

radar_img <- png::readPNG(temp_radar_cropped)
radar_grob <- grid::rasterGrob(radar_img, interpolate = TRUE)
radar_image_plot <- ggplot() + 
  annotation_custom(radar_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

plot_radar_dummy <- plot_radar + coord_cartesian() + 
  scale_color_manual(values=VIRUS_COLORS, guide="none") +
  scale_fill_manual(values=VIRUS_COLORS, breaks=c("V1","V2"), labels=VIRUS_DISPLAY[c("V1","V2")], name="Strain") +
  guides(fill = guide_legend(order=2, override.aes=list(alpha=1)), 
         linetype=guide_legend(order=1, override.aes=list(color="grey25")), 
         linewidth=guide_legend(order=1),
         color = "none")

radar_leg <- cowplot::get_legend(plot_radar_dummy)
plot_radar_safe <- cowplot::plot_grid(NULL, radar_image_plot, NULL, radar_leg, ncol = 1, rel_heights = c(0.02, 1.25, 0.04, 0.16))

# ── 8. Panel F: Heatmap (formerly Panel e) ───────────────────────────────────
MK_ORDER_HM <- c("CD40","CD86","HLA-I","HLA-II","CD68","CD206")
VIR_FULL <- c(V1 = "A/PR/8/34", V2 = "A/NWS/33")
hm <- read.csv(file.path(OUT2, "M2_contrasts.csv"), check.names = FALSE) %>% filter(LPS == "No") %>%
  mutate(Marker = factor(Marker, levels=rev(MK_ORDER_HM)), Virus = factor(VIR_FULL[Virus], levels=VIR_FULL), Tp = tp_fac(Tp))
lim <- max(abs(hm$estimate))
plot_heatmap <- ggplot(hm, aes(Tp, Marker, fill=estimate)) +
  geom_tile(colour="white", linewidth=0.35) + geom_text(aes(label=sig), size=2.8, vjust=0.78, colour="grey15") +
  facet_wrap(~ Virus) + scale_fill_gradient2(low=HM_LOW, mid=HM_MID, high=HM_HIGH, midpoint=0, limits=c(-lim, lim), name = expression("log(OR)"["Infected vs. Bystander"])) +
  guides(fill = guide_colourbar(title.position="top", title.hjust=0.5, barwidth=unit(5.5, "lines"), barheight=unit(0.35, "lines"), reverse = TRUE)) +
  labs(x=NULL, y=NULL) + theme_paper() +
  theme(
    legend.position="bottom",
    panel.grid=element_blank(),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    axis.text.x = ggtext::element_markdown(size = 8.5, face = "bold"),
    legend.title = element_text(size=8, face="bold"),
    legend.text = element_text(size=7),
    legend.margin = margin(t=1, b=0)
  )

bpath <- file.path(OUT2, "bayes_4h_sensing.csv")
if(file.exists(bpath)) {
  bdf <- read.csv(bpath, check.names=FALSE); bsub <- bdf[bdf$Marker %in% c("HLA-I","CD86"), ]
  hl <- expand.grid(Marker=factor(c("HLA-I","CD86"), levels=rev(MK_ORDER_HM)), Tp=factor("4h", levels=c("4h","8h","24h")), Virus=factor(VIR_FULL, levels=VIR_FULL))
  plot_heatmap <- plot_heatmap + geom_tile(data=hl, fill=NA, colour="grey5", linewidth=0.3)
}

# ── 9. Panel G: Paracrine Lines (formerly Panel f) ───────────────────────────
MK_UP3 <- c("CD40","CD68","HLA-I")
em <- read.csv(file.path(OUT2, "M3A_emmeans.csv"), check.names=FALSE) %>% filter(LPS=="No", Marker %in% MK_UP3) %>%
  mutate(Virus=factor(Virus, levels=c("None","V1","V2")), Marker=factor(Marker, levels=MK_UP3), Tp=tp_fac(Tp), CI_low=100*asymp.LCL, CI_high=100*asymp.UCL)
raw <- read.csv(file.path(OUT2, "counts_overall_long.csv"), check.names=FALSE) %>% filter(LPS=="No", Marker %in% MK_UP3) %>%
  mutate(Virus=factor(Virus, levels=c("None","V1","V2")), Marker=factor(Marker, levels=MK_UP3), Tp=tp_fac(Timepoint))
lab_df <- read.csv(file.path(OUT2, "M3A_contrasts.csv"), check.names=FALSE) %>% filter(LPS=="No", Marker %in% MK_UP3) %>%
  mutate(Virus=factor(gsub(" - None","",contrast), levels=c("V1","V2")), Marker=factor(Marker, levels=MK_UP3), Tp=tp_fac(Tp), txt=ifelse(Tp=="24h", paste0("log(OR): ", sprintf("%+.1f", estimate), sig), sig)) %>%
  left_join(em %>% select(Virus, Marker, Tp, y_pct=pct), by=c("Virus","Marker","Tp")) %>%
  filter(txt != "") %>% mutate(xlab=as.integer(Tp)-0.16, ylab=y_pct)
GRP_LAB_VS <- c(None="Naive", V1="A/PR/8/34 Bystander", V2="A/NWS/33 Bystander")

plot_paracrine <- ggplot(em, aes(Tp, pct, colour=Virus, group=Virus)) +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high, fill=Virus), alpha=0.15, colour=NA) +
  geom_jitter(data=raw, aes(Tp, pct, colour=Virus), inherit.aes=FALSE, width=0.14, height=0, size=0.75, alpha=0.4) +
  geom_line(linewidth=0.3) + geom_point(shape=21, fill="white", size=1.1, stroke=0.55) +
  geom_shadowtext(data=lab_df, aes(xlab, ylab, label=txt, colour=Virus), inherit.aes=FALSE, size=2.35, fontface="bold", hjust=1, vjust=0.5, bg.colour="white", bg.r=0.15, show.legend=FALSE) +
  facet_wrap(~ Marker, nrow=1, scales="free_y") +
  scale_colour_manual(values=VIRUS_COLORS, labels=GRP_LAB_VS, name=NULL) +
  scale_fill_manual(values=VIRUS_COLORS, guide="none") +
  scale_x_discrete(expand=expansion(add=c(0.6, 0.25))) + scale_y_continuous(expand=expansion(mult=c(0.10, 0.14))) +
  labs(x=NULL, y="% positive") + theme_paper() +
  theme(
    legend.position="bottom",
    axis.text.x = ggtext::element_markdown(size = 8.5, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    legend.text = element_text(size=7.5, face="plain"),
    legend.key.size = unit(0.4, "lines"),
    legend.margin = margin(t=4, b=0)
  )

# ── 10. Final Layout Assembly ──────────────────────────────────────────────────
# Dimensions adjusted so the 3 data rows fit comfortably below the ~50 mm scheme.

pA <- plot_scheme
pB <- plot_forest + theme(plot.margin = margin(4, 6, 4, 6))
pC <- plot_traj + theme(plot.margin = margin(4, 6, 4, 18))
pD <- plot_biplot + theme(plot.margin = margin(4, 6, 4, 6))
pE <- plot_radar_safe + theme(plot.margin = margin(2, 2, 2, 2))
pF <- plot_heatmap + theme(plot.margin = margin(8, 6, 4, 6))
pG <- plot_paracrine + theme(plot.margin = margin(8, 6, 4, 6))

row_scheme <- cowplot::plot_grid(
  pA,
  ncol = 1,
  labels = "a", label_size = 16, label_fontface = "bold"
)

row1 <- cowplot::plot_grid(
  pB, pC, 
  ncol = 2, align = "h", axis = "b", 
  labels = c("b", "c"), label_size = 16, label_fontface = "bold"
)

row2 <- cowplot::plot_grid(
  pD, pE, 
  ncol = 2, align = "none", 
  labels = c("d", "e"), label_size = 16, label_fontface = "bold"
)

row3 <- cowplot::plot_grid(
  pF, pG, 
  ncol = 2, align = "h", axis = "b", 
  labels = c("f", "g"), label_size = 16, label_fontface = "bold"
)

# Relative heights: Scheme reduced to ~70% area (~0.40);
# the saved ~0.15 height is distributed equally into 3 spacer gaps (0.05 each)
gap_spacer <- NULL

final_layout <- cowplot::plot_grid(
  row_scheme,
  gap_spacer,
  row1,
  gap_spacer,
  row2,
  gap_spacer,
  row3,
  ncol = 1,
  rel_heights = c(0.40, 0.05, 1.0, 0.05, 1.0, 0.05, 0.85)
)

# High-resolution export (600 DPI & 1200 DPI publication quality)
out_pdf <- file.path(OUTDIR, "Figure1.pdf")
out_png_600 <- file.path(OUTDIR, "Figure1.png")
out_png_1200 <- file.path(OUTDIR, "Figure1_1200dpi.png")

message("Saving vector PDF (cairo_pdf)...")
tmp_pdf <- tempfile(fileext = ".pdf")
ggsave(tmp_pdf, final_layout, width = 8.27, height = 11.69, device = cairo_pdf)
file.copy(tmp_pdf, out_pdf, overwrite = TRUE)
unlink(tmp_pdf)

message("Saving high-resolution PNG (600 dpi)...")
tmp_png_600 <- tempfile(fileext = ".png")
ggsave(tmp_png_600, final_layout, width = 8.27, height = 11.69, dpi = 600, bg = "white")
file.copy(tmp_png_600, out_png_600, overwrite = TRUE)
unlink(tmp_png_600)

message("Saving ultra high-resolution PNG (1200 dpi)...")
tmp_png_1200 <- tempfile(fileext = ".png")
ggsave(tmp_png_1200, final_layout, width = 8.27, height = 11.69, dpi = 1200, bg = "white")
file.copy(tmp_png_1200, out_png_1200, overwrite = TRUE)
unlink(tmp_png_1200)

message("DONE: Figure1 created successfully at:\n  ",
        out_pdf, "\n  ", out_png_600, "\n  ", out_png_1200)

