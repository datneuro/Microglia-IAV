# Figure-03-04_Western-Blot.R 

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggsignif); library(grid); library(ggtext)
})
if (!exists("ingest_phase1")) source(file.path("04_Western_Blot", "00_Config-and-Ingestion.R"))
dig <- readRDS(file.path(OUT_DIR, "WB_results_digest.rds"))
sens_data <- read_csv(file.path(OUT_DIR, "WB_sensitivity.csv"), show_col_types = FALSE) %>%
  filter(branch == "A_model_membrane", prior_set == "default")
FIG_OUT <- Sys.getenv("FIGURE_OUT", unset = file.path("Figures", "output"))
dir.create(FIG_OUT, showWarnings = FALSE, recursive = TRUE)

# ---- split-square legend key (light = -LPS, dark = +LPS) ---------------------
PAIR_OF <- function(fc) {
  if (is.na(fc))                          c("grey90","grey60")
  else if (fc %in% c("#D5D8DC","#6E7376")) c("#D5D8DC","#6E7376")
  else if (fc %in% c("#A8C7F0","#4285F4")) c("#A8C7F0","#4285F4")
  else if (fc %in% c("#F2B4B4","#C0392B")) c("#F2B4B4","#C0392B")
  else c(fc, fc)
}
draw_key_split <- function(data, params, size) {
  pr <- PAIR_OF(data$fill)
  grid::grobTree(
    grid::polygonGrob(c(0,0,1), c(0,1,1), gp = grid::gpar(fill = pr[1], col = NA)),
    grid::polygonGrob(c(0,1,1), c(0,0,1), gp = grid::gpar(fill = pr[2], col = NA)),
    grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey50", lwd = 0.7)))
}

p1 <- ingest_phase1()
GRP1 <- c("C","PR8","S33","LPS","LPS/PR8","LPS/S33")
p1 <- p1 %>% mutate(group = factor(group, levels = GRP1))
xidx <- setNames(seq_along(GRP1), GRP1)

PAIR1 <- tribble(
  ~contrast,                              ~x,        ~xend,
  "LPS baseline: C vs LPS",               "C",       "LPS",
  "Infection(-LPS): C vs PR8",            "C",       "PR8",
  "Infection(-LPS): C vs S33",            "C",       "S33",
  "Infection(-LPS): PR8 vs S33",          "PR8",     "S33",
  "Co-stim(+LPS): LPS vs LPS/PR8",        "LPS",     "LPS/PR8",
  "Co-stim(+LPS): LPS vs LPS/S33",        "LPS",     "LPS/S33",
  "Co-stim(+LPS): LPS/PR8 vs LPS/S33",    "LPS/PR8", "LPS/S33")

sig1 <- sens_data %>% filter(phase == "P1") %>%
  rename(target = response) %>%
  left_join(PAIR1, by = "contrast") %>%
  filter(!is.na(x) & pd >= 0.975 & (HDI_low * HDI_high > 0) & prior_stable == TRUE & rope_stable == TRUE) %>%
  filter(target != "syn1.mono") %>%
  mutate(label = sprintf("italic(p)[d]=='%.1f%%'", pd * 100), span = abs(xidx[xend] - xidx[x]))

# Formal per-panel Bayesian coefficient annotation (filtered for pd >= 0.975)
eff_annot <- function(tg) {
  if (tg == "syn1.mono") return("")
  e <- dig$phase1_effects %>% filter(target == tg, pd >= 0.975)
  if (!nrow(e)) return("")
  fmt <- function(lbl, name) {
    r <- e %>% filter(effect == lbl)
    if (!nrow(r)) return(NA_character_)
    sprintf("%s: *β* = %+.2f, *p*<sub>d</sub> = %.1f%%", name, r$median[1], r$pd[1]*100)
  }
  mains <- na.omit(c(fmt("LPS main", "Main effect of LPS"),
                     fmt("Virus main: PR8", "Main effect of A/PR/8/34"),
                     fmt("Virus main: S33", "Main effect of A/NWS/33")))
  ints  <- na.omit(c(fmt("Virus×LPS: PR8", "LPS × A/PR/8/34"),
                     fmt("Virus×LPS: S33", "LPS × A/NWS/33")))
  out <- c(mains, ints)
  if (length(out) == 0) return("")
  paste(out, collapse = "<br>")
}

# targets shown as phospho/total or LC3 ratio -> drop the raw components
TARGETS_DROP <- c("pP65", "P65", "pikb", "ikb", "lc3.1", "lc3.2")
# y-range expansion multiplier: heavy (×3) for tight-spread targets, else ×1.8
Y3X <- c("syn1.mono", "syn211", "psyn", "psmb5", "psmb5.pre", "psmb8",
         "pre.CTSD", "CTSD", "ratio.pP65.P65", "ratio.pIkB.IkB", "k48", "asc")
ymult <- function(tg) if (tg %in% Y3X) 3.0 else 1.8
RATIO_YLAB <- c(ratio.pIkB.IkB = "log₂(p-IκBα / IκBα)",
                ratio.pP65.P65 = "log₂(p-p65 / p65)",
                ratio.LC3      = "log₂(LC3-II / LC3-I)")
y_axis_lab <- function(tg) {
  if (tg %in% names(RATIO_YLAB)) RATIO_YLAB[[tg]]
  else paste0("log₂(", lab_target(tg), " / β-actin)")
}

build_target_panel <- function(tg) {
  d <- filter(p1, target == tg)
  mult <- ymult(tg)
  yr <- range(d$y, na.rm = TRUE); ysp <- max(diff(yr), 1.0)
  # significance brackets WELL above violins (raised to clear kernel tails)
  sg <- sig1 %>% filter(target == tg) %>% arrange(span) %>%
    mutate(yp = yr[2] + ysp * (0.52 + 0.34 * (row_number() - 1)))
  n_brk <- nrow(sg)
  # EXPANDED y-range: head for brackets, foot for +LPS bar; ×mult more whitespace
  head_pad <- ysp * (if (n_brk) 0.30 else 0.55) * mult
  y_top <- (if (n_brk) max(sg$yp) else yr[2]) + head_pad
  y_bot <- yr[1] - ysp * 0.38 * mult
  total <- y_top - y_bot
  # +LPS bar FIXED relative to panel height, set BELOW axis with a clear gap
  bar_y <- y_bot - total * 0.105             # gap from x-axis (no longer touching)
  txt_y <- bar_y - total * 0.030             # label kept snug under the bar

  ggplot(d, aes(group, y)) +
    geom_vline(xintercept = 3.5, linetype = "dotted", color = "grey45", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.3) +
    geom_violin(aes(fill = group), color = NA, alpha = 0.40,
                scale = "width", width = 0.95, trim = FALSE, adjust = 2.2) +
    geom_boxplot(aes(fill = group, color = group), alpha = 0.55, width = 0.4,
                 outlier.shape = NA, lwd = 0.55) +
    geom_jitter(aes(color = group), shape = 21, fill = "white",
                size = 2.4, stroke = 0.7, width = 0.1, alpha = 0.85) +
    { if (n_brk) geom_signif(data = sg, inherit.aes = FALSE,
        aes(xmin = x, xmax = xend, annotations = label, y_position = yp),
        manual = TRUE, tip_length = 0.01, textsize = 2.5, parse = TRUE, color = "grey20")
      else list() } +
    # +LPS bar drawn BELOW the x-axis (outside ylim, clip off), label snug beneath
    annotate("segment", x = 3.6, xend = 6.4, y = bar_y, yend = bar_y,
             color = "#1A1A1A", linewidth = 1.1) +
    annotate("text", x = 5.0, y = txt_y, label = "+ LPS", size = 5.2,
             fontface = "bold", color = "#1A1A1A", vjust = 1) +
    scale_fill_manual(values = PAL_TREAT, guide = "none") +
    scale_color_manual(values = PAL_TREAT, guide = "none") +
    scale_x_discrete(labels = NULL, expand = expansion(add = 0.6)) +
    coord_cartesian(ylim = c(y_bot, y_top), clip = "off") +
    labs(title = lab_target(tg), subtitle = eff_annot(tg), x = NULL, y = y_axis_lab(tg)) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 12.5),
          plot.title = ggtext::element_textbox_simple(
            face = "bold", size = 17, halign = 0.5, fill = "grey85",
            padding = margin(2.5, 3, 2.5, 3), margin = margin(b = 3),
            linetype = 1, r = unit(2, "pt"), box.color = "grey60"),
          plot.subtitle = ggtext::element_markdown(size = 10, color = "grey25",
                                       lineheight = 1.25, hjust = 0, margin = margin(b = 4)),
          plot.margin = margin(t = 4, r = 30, b = 38, l = 30))
}

# standalone split-square legend (3 colors × light/dark)
leg_df <- tibble(g = factor(c("Naive","PR8","S33"), levels = c("Naive","PR8","S33")),
                 fill = c("#D5D8DC","#A8C7F0","#F2B4B4"), y = 1)
p_leg_src <- ggplot(leg_df, aes(g, y, fill = fill)) +
  geom_tile(width = 0.7, height = 0.5, key_glyph = draw_key_split) +
  scale_fill_identity(guide = "legend", name = "Treatment",  # ±LPS note commented out
                      labels = c("#D5D8DC"="Naive","#A8C7F0"="A/PR/8/34",
                                 "#F2B4B4"="A/NWS/33"),
                      breaks = c("#D5D8DC","#A8C7F0","#F2B4B4")) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(face="bold", size=9),
        legend.text = element_text(size=8.5), legend.key.size = unit(0.7,"cm"))
gt <- ggplot_gtable(ggplot_build(p_leg_src))
li <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
p_legend <- if (length(li)) patchwork::wrap_elements(gt$grobs[[li[1]]]) else ggplot() + theme_void()

# Panel order = manuscript §6 (Fig 5): proteasome/proteostasis -> NF-kB -> autophagy -> a-synuclein.
# Precursors placed immediately BEFORE their mature partner; a-syn block mono->multi->211->psyn->25kDa
# (25-kDa last = highlight).
TARGET_ORDER <- c(
  "psmb5.pre", "psmb5", "psmb8", "k48", "pre.CTSD", "CTSD", "asc",   # proteasome / proteostasis
  "ratio.pIkB.IkB", "ratio.pP65.P65",                                # NF-kB signalling
  "ratio.LC3",                                                       # autophagy
  "syn1.mono", "syn1.multi", "syn211", "psyn", "syn1.25"            # a-synuclein (25-kDa last)
)
targets1 <- intersect(TARGET_ORDER, setdiff(unique(p1$target), TARGETS_DROP))  # keep declared order; drop raw p/total & LC3-I/II
panels1 <- lapply(targets1, build_target_panel)
fig1 <- wrap_plots(c(panels1, list(p_legend)), ncol = 3)   # no overall title/subtitle (kept in report md)
ggsave(file.path(FIG_OUT, "WB_Phase1_panels.pdf"), fig1, width = 15, height = 22, device = cairo_pdf)
ggsave(file.path(FIG_OUT, "WB_Phase1_panels.png"), fig1, width = 15, height = 22, dpi = 300)
message("  saved WB_Phase1_panels (", length(targets1), " targets)")

GRP2 <- c("C","PO","PR8","S33")
p2 <- ingest_phase2(); tc <- phase2_timecourse(p2)
STRIP2 <- c(y_main = "α-Synuclein (Syn1)",
            y_aggr = "α-Synuclein 25 kDa (Syn1)",
            y_pocc = "phospho-α-Synuclein")
ROWLAB <- c(y_main = "log₂(α-Syn)", #y_main = "log₂(α-Synuclein (Syn1) / β-actin)",
            y_aggr = "log₂(α-Syn 25 kDa)", #y_main = "log₂(α-Synuclein (Syn1) / β-actin)",
            y_pocc = "log₂(phospho-α-Syn)")     # y_pocc = "log₂(phospho-α-Synuclein / β-actin)")
treat_lab <- c(C="Naive", PO="Poly(I:C)", PR8="A/PR/8/34", S33="A/NWS/33")
COL2 <- setNames(unname(PAL_TREAT[GRP2]), treat_lab[GRP2])

raw_long <- tc %>%
  pivot_longer(all_of(names(STRIP2)), names_to="response", values_to="val") %>%
  filter(is.finite(val), group %in% GRP2) %>%
  mutate(tpx = as.integer(factor(timepoint, levels=c("4h","8h","24h"))),
         response  = factor(response, levels=names(STRIP2), labels=ROWLAB[names(STRIP2)]),
         Treatment = factor(group, levels=GRP2, labels=treat_lab[GRP2]))
sum_long <- raw_long %>% group_by(response, Treatment, tpx) %>%
  summarize(m=mean(val), .groups="drop")

# per-CELL polynomial-trend stats — filter for pd >= 0.975
trend_cell <- dig$phase2_temporal %>%
  filter(branch=="A_model_membrane", kind=="Polynomial_time",
         contrast %in% c("linear","quadratic"), Treatment %in% GRP2, pd >= 0.975) %>%
  mutate(term=sprintf("*β*<sub>%s</sub> = %+.2f (*p*<sub>d</sub> = %.1f%%)", contrast, median, pd * 100)) %>%
  arrange(response, Treatment, desc(contrast=="quadratic")) %>%
  group_by(response, Treatment) %>%
  summarize(lab=paste(term, collapse="<br>"), .groups="drop") %>%
  mutate(response  = factor(response, levels=names(STRIP2), labels=ROWLAB[names(STRIP2)]),
         Treatment = factor(Treatment, levels=GRP2, labels=treat_lab[GRP2]))

# treatment-vs-Naive significance text per timepoint (pd>=0.975 + ROPE & Prior stable)
sig2 <- sens_data %>% filter(phase == "P2" & grepl(" - C$", contrast)) %>%
  filter(pd >= 0.975 & (HDI_low * HDI_high > 0) & prior_stable == TRUE & rope_stable == TRUE) %>%
  mutate(group=sub(" - C$","",contrast)) %>% filter(group %in% GRP2) %>%
  mutate(tpx = as.integer(factor(timepoint, levels=c("4h","8h","24h"))),
         response  = factor(response, levels=names(STRIP2), labels=ROWLAB[names(STRIP2)]),
         Treatment = factor(group, levels=GRP2, labels=treat_lab[GRP2]),
         label = sprintf("italic(p)[d]=='%.1f%%'", pd * 100)) %>%
  left_join(raw_long %>% group_by(response, Treatment, tpx) %>%
              summarize(ytop=max(val), .groups="drop"),
            by=c("response","Treatment","tpx"))

# colored timepoint background bands (applied to every cell)
fig2 <- ggplot(raw_long, aes(tpx, val, fill=Treatment, color=Treatment)) +
  annotate("rect", xmin=0.5, xmax=1.5, ymin=-Inf, ymax=Inf, fill="#FFF3E0", alpha=0.55) +
  annotate("rect", xmin=1.5, xmax=2.5, ymin=-Inf, ymax=Inf, fill="#E8F5E9", alpha=0.55) +
  annotate("rect", xmin=2.5, xmax=3.5, ymin=-Inf, ymax=Inf, fill="#E3F2FD", alpha=0.55) +
  geom_violin(aes(group=tpx), color=NA, alpha=0.35, scale="width", width=0.75,
              trim=FALSE, adjust=1.6) +
  geom_boxplot(aes(group=tpx), width=0.24, alpha=0.55, outlier.shape=NA, lwd=0.45) +
  geom_jitter(width=0.06, height=0, shape=21, fill="white", size=2.5, stroke=0.5, alpha=0.85) + #1.9 -> 4.5 
  geom_line(data=sum_long, aes(tpx, m, group=1), linewidth=0.8) +
  geom_point(data=sum_long, aes(tpx, m), size=2.1) +
  geom_text(data=sig2, aes(tpx, ytop, label=label), inherit.aes=FALSE,
            vjust=-1.2,hjust = 0.65, size=4.5, parse=TRUE, color="grey20") +
  ggtext::geom_richtext(data=trend_cell, aes(x=0.45, y=Inf, label=lab), inherit.aes=FALSE,
       hjust=0, vjust=1.05, size=4.5, fill=NA, label.color=NA, color="grey15", lineheight=1.2) +
  facet_grid(response ~ Treatment, scales="free_y", switch="y") +
  scale_x_continuous(breaks=1:3, labels=c("4h","8h","24h"), expand=expansion(add=0)) +
  scale_y_continuous(expand=expansion(mult=c(0.10, 0.42))) +
  scale_fill_manual(values=COL2, guide="none") +
  scale_color_manual(values=COL2, guide="none") +
  labs(x=NULL, y=NULL) +                       # no overall title/subtitle (kept in report md)
  theme_bw(base_size=11) +
  theme(strip.background.x=element_rect(fill="grey85", color="grey55"),
        strip.background.y=element_blank(),    # row label = plain y-axis title, no grey box
        strip.text.x=element_text(face="bold", size=15),
        strip.text.y.left=element_text(size=15,face = "bold", angle=90),
        strip.placement="outside",
        axis.text.x=element_text(face="bold", size=15),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())
ggsave(file.path(FIG_OUT, "WB_Phase2_aSyn_kinetics.pdf"), fig2, width=12, height=9, device=cairo_pdf)
ggsave(file.path(FIG_OUT, "WB_Phase2_aSyn_kinetics.png"), fig2, width=12, height=9, dpi=300)
message("  saved WB_Phase2_aSyn_kinetics")

# ==============================================================================
# §3  PHASE 2 — 24h transfection-sensing panel {C,T,PO,TPO}
# ==============================================================================
GRP24 <- c("C","T","PO","TPO")
pn <- phase2_panel24(p2) %>%
  pivot_longer(all_of(names(STRIP2)), names_to="response", values_to="val") %>%
  filter(is.finite(val)) %>%
  mutate(group=factor(group, levels=GRP24),
         response=factor(response, levels=names(STRIP2), labels=STRIP2))
pn_sum <- pn %>% group_by(response, group) %>%
  summarize(m=mean(val), sem=sd(val)/sqrt(n()), .groups="drop")
PAL24 <- c(C="#D5D8DC", T="#95A5A6", PO="#2ECC71", TPO="#16A085")

fig3 <- ggplot(pn_sum, aes(group, m, fill=group)) +
  geom_col(width=0.7, color="grey25", linewidth=0.3) +
  geom_errorbar(aes(ymin=m-sem, ymax=m+sem), width=0.25) +
  geom_jitter(data=pn, aes(group, val), inherit.aes=FALSE, shape=21, fill="white",
              width=0.12, size=1.6, stroke=0.5) +
  facet_wrap(~response, scales="free_y") +
  scale_fill_manual(values=PAL24, guide="none") +
  labs(title="Western Blot Phase 2 — 24h transfection-sensing panel",
       subtitle="C=Naive · T=transfection vehicle · PO=Poly I:C (extracellular) · TPO=Poly I:C (transfected, intracellular) · mean ± SEM",
       x=NULL, y="log₂(Target / β-actin)") +
  theme_classic(base_size=11) +
  theme(strip.background=element_rect(fill="grey85", color="grey55"),
        strip.text=element_text(face="bold", size=8),
        plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey35", size=8))
ggsave(file.path(FIG_OUT, "WB_Phase2_panel24.pdf"), fig3, width=12, height=5, device=cairo_pdf)
ggsave(file.path(FIG_OUT, "WB_Phase2_panel24.png"), fig3, width=12, height=5, dpi=300)
message("  saved WB_Phase2_panel24\n[Figure-03-04_Western-Blot] DONE.")
