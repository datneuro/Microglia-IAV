
# ── 0. Dependencies ──────────────────────────────────────────────────────────
.req <- c("tidyverse", "patchwork", "car", "emmeans", "sandwich", "ggsignif")
.missing <- .req[!vapply(.req, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing)) stop("Missing packages: ", paste(.missing, collapse = ", "))
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
  library(car); library(emmeans); library(sandwich); library(ggsignif)
})

set.seed(2026)

# ── 1. Directories & global style ────────────────────────────────────────────
OUT    <- Sys.getenv("FLOW_ELISA_OUT", unset = file.path("Flow-ELISA", "R_analysis", "output"))
OUTDIR <- Sys.getenv("FIGURE_OUT",     unset = file.path("Figures", "output"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
source(file.path("Figures", "_Style.R"))   # theme_paper(), GROUP_COLORS, GROUP_LV ...

# =============================================================================
# 2. PANEL GEOMETRY MAPS
# =============================================================================
# x-axis = 21 bars = 3 timepoints × 7 groups, laid out in three timepoint blocks.
GRP_ORDER <- c("C", "C_PR8", "C_S33", "LPS", "LPS_PR8", "LPS_S33", "PolyIC")
GRP_POS   <- setNames(seq_along(GRP_ORDER), GRP_ORDER)        # within-block 1..7
TP_ORDER  <- c("4", "8", "24")
BLOCK_OFF <- setNames(c(0, 7, 14), TP_ORDER)                  # block offsets
TP_LAB    <- c("4" = "4 h", "8" = "8 h", "24" = "24 h")
xpos <- function(tp, grp) unname(BLOCK_OFF[as.character(tp)] + GRP_POS[grp])

# timepoint block background tints + matching header-text colours
TP_FILL     <- c("4" = "#4285F4", "8" = "#F39C12", "24" = "#C0392B")
TP_TEXTCOL  <- c("4" = "#2C5AA0", "8" = "#B9770E", "24" = "#922B21")
BLOCK_XMIN  <- c("4" = 0.4,  "8" = 7.5,  "24" = 14.5)
BLOCK_XMAX  <- c("4" = 7.5,  "8" = 14.5, "24" = 21.6)
BLOCK_XC    <- c("4" = 4,    "8" = 11,   "24" = 18)

# Split-square legend glyph (upper-left = light −LPS · lower-right = dark +LPS),
# inherited house style from the qPCR master figure.
draw_key_split_grp <- function(data, params, size) {
  fc <- data$fill
  pair <- if (is.na(fc))                       c("grey90", "grey60") else
          if (fc %in% c("#D5D8DC", "#6E7376")) c("#D5D8DC", "#6E7376") else
          if (fc %in% c("#A8C7F0", "#4285F4")) c("#A8C7F0", "#4285F4") else
          if (fc %in% c("#F2B4B4", "#C0392B")) c("#F2B4B4", "#C0392B") else
                                               c(fc, fc)
  grid::grobTree(
    grid::polygonGrob(c(0, 0, 1), c(0, 1, 1), gp = grid::gpar(fill = pair[1], col = NA)),
    grid::polygonGrob(c(0, 1, 1), c(0, 0, 1), gp = grid::gpar(fill = pair[2], col = NA)),
    grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey50", lwd = 0.7)))
}
LEG_BREAKS <- c("C", "C_PR8", "C_S33", "PolyIC")
LEG_LABELS <- c("Naive", "A/PR/8/34", "A/NWS/33", "PolyI:C")

# Marker layout — Flow (rows 1–2) then ELISA (rows 3–4); Anti → Pro within each.
#   Flow : CD206 (M2/anti) · CD68 (pan/other) · CD40 CD86 HLA-I HLA-II (M1/pro)
#   ELISA: IL-10 (anti)    · CCL2 (chemokine) · TNF-a IL-1b IL-6 IL-12p70 (pro)
MARKER_ORDER <- c("CD206", "CD68", "CD40", "CD86", "HLA-I", "HLA-II",
                  "IL-10", "CCL2", "TNF-a", "IL-1b", "IL-6", "IL-12p70")

# Greek display + y-axis title per platform
greek <- function(m) dplyr::recode(m, "TNF-a" = "TNF-α", "IL-1b" = "IL-1β", .default = m)
ylab_for <- function(m, platform)
  if (platform == "flow") paste0("%", greek(m)) else paste0(greek(m), " (pg/mL)")

# Planned-contrast → the two groups it spans (for bracket placement)
CONTRAST_PAIR <- list(
  "Ctrl-LPS"      = c("C",       "LPS"),
  "Ctrl-V1"       = c("C",       "C_PR8"),
  "Ctrl-V2"       = c("C",       "C_S33"),
  "V1-V2"         = c("C_PR8",   "C_S33"),
  "LPS-LPS.V1"    = c("LPS",     "LPS_PR8"),
  "LPS-LPS.V2"    = c("LPS",     "LPS_S33"),
  "LPS.V1-LPS.V2" = c("LPS_PR8", "LPS_S33"),
  "PolyIC-Ctrl"   = c("PolyIC",  "C"),
  "PolyIC-LPS"    = c("PolyIC",  "LPS"),
  "PolyIC-V1"     = c("PolyIC",  "C_PR8"),
  "PolyIC-V2"     = c("PolyIC",  "C_S33")
)

# ANOVA terms, canonical order + compact display names
TERM_ORDER <- c("Timepoint", "LPS", "Virus",
                "Timepoint:LPS", "Timepoint:Virus", "LPS:Virus",
                "Timepoint:LPS:Virus")
TERM_DISP  <- c("Timepoint" = "Time", "LPS" = "LPS", "Virus" = "Virus",
                "Timepoint:LPS" = "Time×LPS", "Timepoint:Virus" = "Time×Virus",
                "LPS:Virus" = "LPS×Virus", "Timepoint:LPS:Virus" = "Time×LPS×Virus")

# =============================================================================
# 3. STAT FORMATTERS
# =============================================================================
p_to_stars <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}
fmt_p    <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
fmt_f    <- function(f) ifelse(is.na(f), "NA", sprintf("%.2f", f))
fmt_eta2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))

# =============================================================================
# 4. DATA
# =============================================================================
fac <- read_csv(file.path(OUT, "tidy_factorial.csv"), show_col_types = FALSE) %>%
  filter(Marker != "NP")
poly <- read_csv(file.path(OUT, "tidy_polyic.csv"), show_col_types = FALSE) %>%
  filter(Marker != "NP")

platform_of <- bind_rows(fac, poly) %>% distinct(Marker, Platform) %>% deframe()

# factor levels for the model (as in 01_Univariate-Factorial.R)
mk_factors <- function(d) d %>%
  mutate(Timepoint = factor(Timepoint, levels = c(4, 8, 24)),
         LPS       = factor(LPS,       levels = c("No", "Yes")),
         Virus     = factor(Virus,     levels = c("None", "V1", "V2")),
         Bio_Rep   = factor(Bio_Rep))
df_fact <- fac  %>% mk_factors() %>%
  mutate(Value_t = ifelse(Platform == "flow", log(Value / (100 - Value)), log2(Value)))
df_poly <- poly %>%
  mutate(Timepoint = factor(Timepoint, levels = c(4, 8, 24)), Bio_Rep = factor(Bio_Rep),
         Value_t = ifelse(Platform == "flow", log(Value / (100 - Value)), log2(Value)))

# raw (untransformed) pool for plotting bars + dots (all 7 groups)
raw_all <- bind_rows(fac, poly) %>%
  mutate(Group     = factor(Group, levels = GRP_ORDER),
         Timepoint = as.character(Timepoint))

# =============================================================================
# 5. STAGE-1 STATISTICS, RE-COMPUTED PER MARKER  (identical to stage1_L1_factorial.R)
# =============================================================================
old_opts <- options(contrasts = c("contr.sum", "contr.poly"))

build_contrasts_from_emm <- function(emm) {
  g   <- as.data.frame(summary(emm))
  tp1 <- as.character(g$Timepoint)[1]
  sub <- g[as.character(g$Timepoint) == tp1, c("LPS", "Virus")]
  key <- paste(as.character(sub$LPS), as.character(sub$Virus), sep = ".")
  pos <- function(lps, vir) which(key == paste(lps, vir, sep = "."))
  mk  <- function(plus, minus) { v <- numeric(length(key)); v[plus] <- 1; v[minus] <- -1; v }
  list(
    "Ctrl-LPS"      = mk(pos("No",  "None"), pos("Yes", "None")),
    "Ctrl-V1"       = mk(pos("No",  "None"), pos("No",  "V1")),
    "Ctrl-V2"       = mk(pos("No",  "None"), pos("No",  "V2")),
    "V1-V2"         = mk(pos("No",  "V1"),   pos("No",  "V2")),
    "LPS-LPS.V1"    = mk(pos("Yes", "None"), pos("Yes", "V1")),
    "LPS-LPS.V2"    = mk(pos("Yes", "None"), pos("Yes", "V2")),
    "LPS.V1-LPS.V2" = mk(pos("Yes", "V1"),   pos("Yes", "V2"))
  )
}
ref_groups <- list(
  "PolyIC-Ctrl" = list(LPS = "No",  Virus = "None"),
  "PolyIC-LPS"  = list(LPS = "Yes", Virus = "None"),
  "PolyIC-V1"   = list(LPS = "No",  Virus = "V1"),
  "PolyIC-V2"   = list(LPS = "No",  Virus = "V2")
)

analyze_marker <- function(m) {
  dat <- df_fact %>% filter(Marker == m)
  mod <- lm(Value_t ~ Timepoint * LPS * Virus + Timepoint:Bio_Rep, data = dat)

  anova_std <- car::Anova(mod, type = 3)
  anova_hc3 <- car::Anova(mod, type = 3, white.adjust = "hc3")
  ss_resid  <- anova_std["Residuals", "Sum Sq"]
  df_resid  <- df.residual(mod)

  anova_out <- data.frame(
    Term         = TERM_ORDER,
    Df           = anova_std[TERM_ORDER, "Df"],
    Df_resid     = df_resid,
    partial_eta2 = anova_std[TERM_ORDER, "Sum Sq"] /
                     (anova_std[TERM_ORDER, "Sum Sq"] + ss_resid),
    F_hc3        = anova_hc3[TERM_ORDER, "F"],
    p_hc3        = anova_hc3[TERM_ORDER, "Pr(>F)"],
    row.names = NULL, stringsAsFactors = FALSE
  )

  # planned factorial contrasts (HC3 vcov)
  vcov_hc3 <- sandwich::vcovHC(mod, type = "HC3")
  emm <- emmeans::emmeans(mod, ~ LPS:Virus | Timepoint, vcov. = vcov_hc3)
  con <- as.data.frame(summary(
    emmeans::contrast(emm, method = build_contrasts_from_emm(emm))))
  con_out <- data.frame(Contrast = con$contrast,
                        Timepoint = as.character(con$Timepoint),
                        p_raw = con$p.value, family = "factorial",
                        stringsAsFactors = FALSE)

  # Poly(I:C) Welch t-tests
  dat_p <- df_poly %>% filter(Marker == m)
  poly_rows <- list(); k <- 1
  for (tp_val in levels(dat$Timepoint)) {
    pv <- dat_p %>% filter(Timepoint == tp_val) %>% pull(Value_t)
    for (cn in names(ref_groups)) {
      r  <- ref_groups[[cn]]
      rv <- dat %>% filter(Timepoint == tp_val, LPS == r$LPS, Virus == r$Virus) %>% pull(Value_t)
      tt <- t.test(pv, rv)
      poly_rows[[k]] <- data.frame(Contrast = cn, Timepoint = tp_val,
                                   p_raw = tt$p.value, family = "polyic",
                                   stringsAsFactors = FALSE); k <- k + 1
    }
  }
  con_all <- bind_rows(con_out, bind_rows(poly_rows)) %>%
    group_by(family) %>% mutate(p_holm = p.adjust(p_raw, method = "holm")) %>% ungroup()

  list(anova = anova_out, contrasts = con_all)
}

message("Computing factorial statistics for ", length(MARKER_ORDER), " markers ...")
STATS <- setNames(lapply(MARKER_ORDER, analyze_marker), MARKER_ORDER)
options(old_opts)

# =============================================================================
# 6. ANNOTATION BUILDERS
# =============================================================================
# (a) Omnibus stat block — one line per term, F / p / η²ₚ + stars
omnibus_block <- function(m) {
  a <- STATS[[m]]$anova
  lines <- vapply(TERM_ORDER, function(tm) {
    r <- a[a$Term == tm, ]
    stars <- p_to_stars(r$p_hc3)
    sfx   <- if (!stars %in% c("ns", "")) paste0(" ", stars) else ""
    p_str <- if (!is.na(r$p_hc3) && r$p_hc3 < 0.001) "p<0.001"
             else paste0("p=", fmt_p(r$p_hc3))
    sprintf("%-14s F(%d,%d)=%s, %s, η²p=%s%s",
            TERM_DISP[[tm]], as.integer(r$Df), as.integer(r$Df_resid),
            fmt_f(r$F_hc3), p_str, fmt_eta2(r$partial_eta2), sfx)
  }, character(1))
  paste(lines, collapse = "\n")
}

# (b) Significance brackets — significant Holm contrasts (both families) per panel.
#     Returns full geometry (horizontal bar + tips + centred star position).
build_brackets <- function(m, y_base, y_step) {
  cdf <- STATS[[m]]$contrasts %>% filter(!is.na(p_holm), p_holm < 0.05)
  if (nrow(cdf) == 0)
    return(tibble(xmin = numeric(), xmax = numeric(), xc = numeric(),
                  y_position = numeric(), tip = numeric(), ty = numeric(),
                  label = character()))
  cdf %>%
    rowwise() %>%
    mutate(pair = list(CONTRAST_PAIR[[Contrast]]),
           x1   = xpos(Timepoint, pair[[1]]),
           x2   = xpos(Timepoint, pair[[2]]),
           xmin = min(x1, x2), xmax = max(x1, x2),
           span = abs(x2 - x1),
           label = p_to_stars(p_holm)) %>%
    ungroup() %>%
    group_by(Timepoint) %>%
    arrange(span, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    mutate(y_position = y_base + y_step * rank,
           xc  = (xmin + xmax) / 2,
           tip = y_position - y_step * 0.28,
           ty  = y_position + y_step * 0.10) %>%
    select(xmin, xmax, xc, y_position, tip, ty, label)
}

# =============================================================================
# 7. PANEL BUILDER
# =============================================================================
build_panel <- function(m) {
  platform <- platform_of[[m]]
  rd <- raw_all %>% filter(Marker == m) %>% mutate(x = xpos(Timepoint, Group))

  summ <- rd %>% group_by(Timepoint, Group, x) %>%
    summarise(mean = mean(Value), sem = sd(Value) / sqrt(n()), .groups = "drop")

  # data ceiling = tallest of (bar+SEM) or raw dot, per panel
  y_top_data <- max(c(summ$mean + summ$sem, rd$Value), na.rm = TRUE)
  y_step     <- y_top_data * 0.07

  brk     <- build_brackets(m, y_top_data * 1.02, y_step)
  brk_top <- if (nrow(brk)) max(brk$y_position) else y_top_data * 1.02
  time_y  <- brk_top + y_step * 1.7          # timepoint header sits above brackets
  y_top   <- time_y + y_step * 0.9

  # bottom indicator geometry (negative y → drawn into the bottom margin, clip off)
  rng   <- y_top
  bar_y <- -rng * 0.040
  txt_y <- -rng * 0.090
  off   <- c(0, 7, 14)
  lps_seg  <- data.frame(x = off + 3.6, xend = off + 6.4, y = bar_y)
  poly_seg <- data.frame(x = off + 6.6, xend = off + 7.4, y = bar_y)
  lps_txt  <- data.frame(x = off + 5.0, y = txt_y, label = "+ LPS")
  poly_txt <- data.frame(x = off + 8.0, y = txt_y, label = "PolyI:C") #Origincal + 7

  p <- ggplot() +
    # ── timepoint block background tints (one colour per timepoint) + separators ──
    annotate("rect", xmin = unname(BLOCK_XMIN[TP_ORDER]),
             xmax = unname(BLOCK_XMAX[TP_ORDER]), ymin = -Inf, ymax = Inf,
             fill = unname(TP_FILL[TP_ORDER]), alpha = 0.08) +
    geom_vline(xintercept = c(7.5, 14.5), linetype = "dotted",
               color = "grey55", linewidth = 0.4) +

    # ── bars: group-coloured fill (light) + group-coloured border ──
    geom_col(data = summ, aes(x = x, y = mean, fill = Group, color = Group),
             width = 0.62, alpha = 0.55, linewidth = 0.45,
             key_glyph = draw_key_split_grp) +
    # SEM whisker in the group colour
    geom_errorbar(data = summ, aes(x = x, ymin = pmax(mean - sem, 0),
                                   ymax = mean + sem, color = Group),
                  width = 0.26, linewidth = 0.5) +
    # raw replicate dots
    geom_jitter(data = rd, aes(x = x, y = Value, color = Group),
                shape = 21, fill = "white", stroke = 0.55, size = 2.3,
                width = 0.12, height = 0, alpha = 0.9) +

    { if (nrow(brk)) list(
        geom_segment(data = brk, aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
                     color = "grey20", linewidth = 0.4, inherit.aes = FALSE),
        geom_segment(data = brk, aes(x = xmin, xend = xmin, y = y_position, yend = tip),
                     color = "grey20", linewidth = 0.4, inherit.aes = FALSE),
        geom_segment(data = brk, aes(x = xmax, xend = xmax, y = y_position, yend = tip),
                     color = "grey20", linewidth = 0.4, inherit.aes = FALSE),
        geom_text(data = brk, aes(x = xc, y = ty, label = label),
                  color = "grey20", size = 3.4, vjust = 0, inherit.aes = FALSE)
      ) else list() } +

    # ── timepoint headers, on top of the bars, coloured per block ──
    annotate("text", x = unname(BLOCK_XC[TP_ORDER]), y = time_y,
             label = unname(TP_LAB[TP_ORDER]), color = unname(TP_TEXTCOL[TP_ORDER]),
             fontface = "bold", size = 4.3) +

    # ── x-axis indicator bars: +LPS (groups 4–6) and PolyI:C (group 7) per block ──
    geom_segment(data = lps_seg,  aes(x = x, xend = xend, y = y, yend = y),
                 color = "#1A1A1A", linewidth = 1.0, inherit.aes = FALSE) +
    geom_text(data = lps_txt,  aes(x = x, y = y, label = label), vjust = 1,
              color = "#1A1A1A", fontface = "bold", size = 4.5, inherit.aes = FALSE) +
    geom_segment(data = poly_seg, aes(x = x, xend = xend, y = y, yend = y),
                 color = "#2ECC71", linewidth = 1.0, inherit.aes = FALSE) +
    geom_text(data = poly_txt, aes(x = x, y = y, label = label), vjust = 1,
              color = "#2ECC71", fontface = "bold", size = 4.5, inherit.aes = FALSE) +

    scale_fill_manual(values = GROUP_COLORS, breaks = LEG_BREAKS,
                      labels = LEG_LABELS, drop = FALSE, name = NULL) +
    scale_color_manual(values = GROUP_COLORS, breaks = GRP_ORDER,
                       drop = FALSE, guide = "none") +
    scale_x_continuous(expand = expansion(add = 0.6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, y_top), clip = "off") +
    labs(x = NULL, y = ylab_for(m, platform)) +
    theme_paper(base = 11) + 
    theme(
      plot.subtitle = element_text(size = 7.2, family = "mono", color = "grey25",
                                   hjust = 0, lineheight = 1.2,
                                   margin = margin(t = 0, b = 10)),
      axis.title.y  = element_text(face = "bold", size = 20),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      plot.margin   = margin(t = 45, r = 24, b = 65, l = 22), 
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)))

  p
}

# =============================================================================
# 8. ASSEMBLE — 4 × 3 grid, shared bottom legend
# =============================================================================
message("Building panels ...")
panels <- lapply(MARKER_ORDER, build_panel)

fig <- wrap_plots(panels, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5,
                                   margin = margin(b = 6)))) &
  theme(legend.position   = "bottom",
        legend.title      = element_blank(),
        legend.text       = element_text(size = 12),
        legend.key.size   = unit(0.9, "cm"),
        legend.key.width  = unit(1.1, "cm"),
        legend.box.margin = margin(t = 85, b = 6),   # push legend clear of panels above
        legend.box.spacing = unit(0.6, "cm"))

fig_w <- 21.0
fig_h <- 29.7 # Original 26 
ggsave(file.path(OUTDIR, "FigureS2_marker_barplots.pdf"), fig,
       width = fig_w, height = fig_h, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUTDIR, "FigureS2_marker_barplots.png"), fig,
       width = fig_w, height = fig_h, dpi = 400, limitsize = FALSE)

# ── caption ───────────────────────────────────────────────────────────────────
caption <- paste(
  "Supplementary Fig. S2 · Univariate Factorial Response of Microglial Surface and Secreted Markers.",
  "",
  "Twelve markers (NP excluded) arranged Flow cytometry (top two rows) then ELISA",
  "(bottom two rows); within each platform panels run anti-inflammatory →",
  "pro-inflammatory. Each panel shows 21 bars: 3 timepoints (4/8/24 h) × 7 groups",
  "(Naive, A/PR/8/34, A/NWS/33, LPS, LPS+A/PR/8/34, LPS+A/NWS/33, PolyI:C).",
  "Bars are mean ± SEM on the raw measurement scale (Flow: % of live cells, y-axis",
  "'%<Marker>'; ELISA: pg/mL, y-axis '<Marker> (pg/mL)'); bar border and SEM whisker",
  "carry the group colour and individual biological replicates are overlaid as dots.",
  "",
  "Statistics were computed with the model of 01_Univariate-Factorial.R, lm(transform(Value) ~",
  "Timepoint × LPS × Virus + Timepoint:Bio_Rep) under sum-to-zero contrasts with",
  "HC3 robust standard errors (flow logit-, ELISA log2-transformed). The stat block",
  "above each panel reports the Type-III F, p and partial eta-squared for every main",
  "effect and interaction. Brackets mark Holm-adjusted planned contrasts within each",
  "timepoint (factorial and PolyI:C families corrected separately): * p<0.05,",
  "** p<0.01, *** p<0.001.",
  sep = "\n")
writeLines(caption, file.path(OUTDIR, "FigureS2_caption.txt"))

message("✓ Supplementary Fig. S2 written to: ", OUTDIR,
        "  (", fig_w, " × ", fig_h, " in)  | caption -> FigureS2_caption.txt")
