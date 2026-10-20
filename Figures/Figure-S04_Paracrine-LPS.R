OUT <- Sys.getenv("BYSTANDER_OUT", unset = file.path("Flow-ELISA", "R-analysis-figure_VS-bystander", "output"))
FIG <- Sys.getenv("FIGURE_OUT",    unset = file.path("Figures", "output"))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
source(file.path("01_Flow_ELISA", "infected_vs_bystander", "_Style-Bystander.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2);
  library(patchwork); library(forcats); library(shadowtext) })

rd <- function(f) read.csv(file.path(OUT, f), check.names = FALSE)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 320, bg = "white")
  cat("saved", name, "\n")
}
MK_ORDER <- c("CD40","CD86","HLA-I","HLA-II","CD68","CD206")
VIR_FULL <- c(V1 = "A/PR/8/34", V2 = "A/NWS/33")
GRP_COL  <- c(None = "#6E7376", V1 = "#4285F4", V2 = "#C0392B")
GRP_LAB  <- c(None = "Naive", V1 = "A/PR/8/34 Bystander", V2 = "A/NWS/33 Bystander")
MK_UP3   <- c("CD40","CD68","HLA-I")     # main panel B (clear monotonic ↑)
MK_DOWN  <- c("CD86","CD206")            # supplementary (regulatory ↓)

make_heatmap <- function(lps, bayes = FALSE) {
  hm <- rd("M2_contrasts.csv") %>% filter(LPS == lps) %>%
    mutate(Marker = factor(Marker, levels = rev(MK_ORDER)),
           Virus  = factor(VIR_FULL[Virus], levels = VIR_FULL),
           Tp     = factor(paste0(Tp,"h"), levels = TP_LV))
  lim <- max(abs(hm$estimate))
  g <- ggplot(hm, aes(Tp, Marker, fill = estimate)) +
    geom_tile(colour = "white", linewidth = .6) +
    geom_text(aes(label = sig), size = 3.4, vjust = .78, colour = "grey15") +
    facet_wrap(~ Virus) +
    scale_fill_gradient2(low = HM_LOW, mid = HM_MID, high = HM_HIGH, midpoint = 0,
        limits = c(-lim, lim), name = expression(log(OR)["Infected vs. Bystander"])) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = .5,
                                  barwidth = unit(9, "lines"), barheight = unit(.6, "lines"))) +
    labs(x = NULL, y = NULL) +
    theme_paper() + theme(legend.position = "bottom", panel.grid = element_blank())

  bpath <- file.path(OUT, "bayes_4h_sensing.csv")
  if (bayes && file.exists(bpath)) {
    bdf  <- read.csv(bpath, check.names = FALSE)
    bsub <- bdf[bdf$Marker %in% c("HLA-I","CD86"), ]
    hl <- expand.grid(Marker = factor(c("HLA-I","CD86"), levels = rev(MK_ORDER)),
                      Tp = factor("4h", levels = TP_LV),
                      Virus = factor(VIR_FULL, levels = VIR_FULL))
    cap <- paste0(
      "Bayesian 4h sensing (NP+ vs NP−, beta-binomial MCMC; boxed cells):\n",
      paste(sprintf("%s: log-OR %.2f  ·  95%% HDI [%.2f, %.2f]  ·  BF₁₀ %s  ·  BF₀₁ %s",
                    bsub$Marker, bsub$med_logOR, bsub$HDI_low, bsub$HDI_high,
                    bsub$BF10_fmt, bsub$BF01_fmt), collapse = "\n"),
      "\nHDI excludes 0 → robust to small NP+ gate (no asymptotics); BF prior N(0,1.5).")
    g <- g +
      geom_tile(data = hl, fill = NA, colour = "grey5", linewidth = .45) +
      labs(caption = cap) +
      theme(plot.caption = element_text(hjust = 0, size = 7, colour = "grey25",
                                        lineheight = 1.15, margin = margin(t = 7)))
  }
  g
}

make_paracrine_lines <- function(lps, markers, nrow = 1) {
  em <- rd("M3A_emmeans.csv") %>% filter(LPS == lps, Marker %in% markers) %>%
    mutate(Virus = factor(Virus, levels = c("None","V1","V2")),
           Marker = factor(Marker, levels = markers),
           Tp = factor(paste0(Tp,"h"), levels = TP_LV),
           CI_low = 100*asymp.LCL, CI_high = 100*asymp.UCL)
  raw <- rd("counts_overall_long.csv") %>% filter(LPS == lps, Marker %in% markers) %>%
    mutate(Virus = factor(Virus, levels = c("None","V1","V2")),
           Marker = factor(Marker, levels = markers),
           Tp = factor(paste0(Timepoint,"h"), levels = TP_LV))
  # Labels sit immediately to the LEFT of each strain's own point (so the colour
  # follows the line, not a fixed top/bottom) and clear of the raw-point cloud.
  # 24 h: full "log(OR): ±x.x" + stars · 4 h & 8 h: stars only.
  lab <- rd("M3A_contrasts.csv") %>% filter(LPS == lps, Marker %in% markers) %>%
    mutate(Virus = factor(gsub(" - None","",contrast), levels = c("V1","V2")),
           Marker = factor(Marker, levels = markers),
           Tp = factor(paste0(Tp,"h"), levels = TP_LV),
           txt = ifelse(Tp == "24h",
                        paste0("log(OR): ", sprintf("%+.1f", estimate), sig), sig)) %>%
    left_join(em %>% select(Virus, Marker, Tp, y_pct = pct),
              by = c("Virus","Marker","Tp")) %>%
    filter(txt != "") %>%
    mutate(xlab = as.integer(Tp) - 0.14, ylab = y_pct)

  ggplot(em, aes(Tp, pct, colour = Virus, group = Virus)) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = Virus), alpha = .15, colour = NA) +
    geom_jitter(data = raw, aes(Tp, pct, colour = Virus), inherit.aes = FALSE,
                width = .13, height = 0, size = .85, alpha = .35) +
    geom_line(linewidth = .45) +
    geom_point(shape = 21, fill = "white", size = 1.9, stroke = .9) +
    geom_shadowtext(data = lab, aes(xlab, ylab, label = txt, colour = Virus),
                    inherit.aes = FALSE, size = 2.6, fontface = "bold",
                    hjust = 1, vjust = 0.5, bg.colour = "white", bg.r = .18,
                    show.legend = FALSE) +
    facet_wrap(~ Marker, nrow = nrow, scales = "free_y") +
    scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
    scale_fill_manual(values = GRP_COL, guide = "none") +
    scale_x_discrete(expand = expansion(add = c(0.6, 0.25))) +
    scale_y_continuous(expand = expansion(mult = c(.10, .14))) +
    labs(x = NULL, y = "% positive") + theme_paper()
}

build_2panel <- function(lps) {
  (make_heatmap(lps, bayes = (lps == "No")) | make_paracrine_lines(lps, MK_UP3, nrow = 1)) +
    plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "a")
}
save_fig(build_2panel("No"),  "Fig1fg_noLPS_standalone", 14, 5.9)
save_fig(build_2panel("Yes"), "FigureS4_withLPS",         14, 5.2)

# individual halves (both arms)
save_fig(make_heatmap("No"),                       "PanelLeft_heatmap_noLPS",   7.2, 4.6)
save_fig(make_heatmap("Yes"),                      "PanelLeft_heatmap_LPS",     7.2, 4.6)
save_fig(make_paracrine_lines("No",  MK_UP3, 1),   "PanelRight_paracrine_noLPS", 7.5, 4.2)
save_fig(make_paracrine_lines("Yes", MK_UP3, 1),   "PanelRight_paracrine_LPS",   7.5, 4.2)

pDown <- make_paracrine_lines("No", MK_DOWN, nrow = 1)

G4 <- c("C_PR8","LPS_PR8","C_S33","LPS_S33")
kf <- rd("kinetics_emmeans.csv") %>%
  mutate(Group = case_when(Virus=="V1"&LPS=="No"~"C_PR8", Virus=="V1"&LPS=="Yes"~"LPS_PR8",
                           Virus=="V2"&LPS=="No"~"C_S33", Virus=="V2"&LPS=="Yes"~"LPS_S33"),
         Group = factor(Group, levels = G4), Tp = factor(paste0(Tp,"h"), levels = TP_LV))
pF1 <- ggplot(kf, aes(Tp, pctNP, colour = Group, group = Group)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = Group), alpha = .12, colour = NA) +
  geom_line(linewidth = .7) + geom_point(shape = 21, fill = "white", size = 2, stroke = .9) +
  scale_colour_manual(values = GROUP_COLORS[G4], labels = lab_for(G4), name = NULL) +
  scale_fill_manual(values = GROUP_COLORS[G4], guide = "none") +
  labs(x = NULL, y = "% NP+") + theme_paper()
vf <- rd("viability_summary.csv") %>%
  filter(Group %in% c("C","C_PR8","C_S33","LPS_PR8","LPS_S33")) %>%
  mutate(Group = factor(Group, levels = GROUP_LV), Tp = factor(paste0(Timepoint,"h"), levels = TP_LV))
pF2 <- ggplot(vf, aes(Tp, viability, colour = Group, group = Group)) +
  geom_hline(yintercept = 80, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = .7) +
  geom_pointrange(aes(ymin = viability-se, ymax = viability+se), size = .3) +
  scale_colour_manual(values = GROUP_COLORS, labels = lab_for(levels(vf$Group)), name = NULL) +
  coord_cartesian(ylim = c(75, 100)) +
  labs(x = NULL, y = "% viable") + theme_paper()

figComp <- (pDown / (pF1 | pF2)) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a")
save_fig(figComp, "FigureS4_companion_panels", 11, 8)

cat("\n═══ Figure-S04_Paracrine-LPS COMPLETE ═══\n")
