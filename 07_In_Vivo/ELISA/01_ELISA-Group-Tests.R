#!/usr/bin/env Rscript
# =============================================================================
# ELISA
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(cowplot); library(ggsignif); library(scales)
})
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)))

xlsx_name  <- "20260901 ELISA Lung and Brain 5day post infection.xlsx"
sheet_name <- "R_data"
out_fig    <- file.path("ELISA", "figures", "bar_figures")   
out_tab    <- file.path("ELISA", "tables")

group_levels <- c("Control", "PR8", "S33")
COL <- c(Control = "#9E9E9E", PR8 = "#4C72B0", S33 = "#C44E52")
XLB <- c(Control = "Naive",   PR8 = "A/PR/8/34", S33 = "A/NWS/33")

tissue_levels <- c("Lung", "Brain")
unit_default  <- "pg/mL"
unit_override <- c()       

star_size    <- 5.2           
dot_size     <- 3.2         
dot_stroke   <- 1.0
bold_ylab    <- TRUE          
panel_labels <- FALSE       
save_panels  <- TRUE          

fig_ncol   <- 5             
panel_w    <- 2.9            
panel_h    <- 3.3
dpi        <- 1200           
# -----------------------------------------------------------------------------


## ---- read data ---------------------------------------------------------------
find_file <- function(f) {
  cand <- c(f, file.path("ELISA", f), file.path("..", f))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("File not found: ", f, call. = FALSE)
  hit[1]
}
read_excel_safe <- function(path, ...) {
  out <- try(read_excel(path, ...), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  tmp <- file.path(tempdir(), "elisa_tmp.xlsx"); file.copy(path, tmp, overwrite = TRUE)
  read_excel(tmp, ...)
}
# ggplot2 4.x defaults to the ragg device, which cannot write to some paths
save_png <- function(filename, plot, ...) ggsave(filename, plot, ..., device = grDevices::png)

dat <- read_excel_safe(find_file(xlsx_name), sheet = sheet_name) %>%
  select(group = Group, Tissue, unit = Replicate, Cytokine, v = starts_with("Conc")) %>%
  mutate(across(c(group, Tissue, Cytokine), ~ str_squish(as.character(.x))),
         v = suppressWarnings(as.numeric(v))) %>%
  filter(!is.na(v), !is.na(group), !is.na(Cytokine))

cyto_levels <- unique(dat$Cytokine)
dat <- dat %>% mutate(group    = factor(group, levels = group_levels),
                      Tissue   = factor(Tissue, levels = tissue_levels),
                      Cytokine = factor(Cytokine, levels = cyto_levels))

y_title <- function(cyt) {
  u <- if (cyt %in% names(unit_override)) unit_override[[cyt]] else unit_default
  paste0(cyt, " (", u, ")")
}


pairwise_lm <- function(r) {
  r <- droplevels(r)
  m  <- lm(v ~ group, data = r)
  df <- m$df.residual
  mse <- sum(residuals(m)^2) / df
  mu <- tapply(r$v, r$group, mean); nn <- tapply(r$v, r$group, length)
  lv <- levels(r$group); cmb <- combn(lv, 2)
  pv <- setNames(vapply(seq_len(ncol(cmb)), function(i) {
    a <- cmb[1, i]; b <- cmb[2, i]
    se <- sqrt(mse * (1 / nn[[a]] + 1 / nn[[b]]))
    2 * pt(-abs((mu[[a]] - mu[[b]]) / se), df)
  }, numeric(1)), paste(cmb[1, ], "-", cmb[2, ]))
  list(p = pv, emmean = mu, p_omni = anova(m)$`Pr(>F)`[1])
}

sig_lab <- function(p) if (is.null(p) || is.na(p)) "" else
  if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else ""


bar_panel <- function(r, ylab) {
  st <- pairwise_lm(r)
  pv <- st$p
  su <- r %>% group_by(group) %>%
    summarise(M = mean(v), SE = sd(v) / sqrt(n()), .groups = "drop") %>%
    mutate(group = factor(group, levels = group_levels))

  top  <- max(c(su$M + su$SE, r$v), na.rm = TRUE)
  nmc  <- c("Control - PR8", "Control - S33", "PR8 - S33")
  ann  <- vapply(nmc, function(k) sig_lab(pv[[k]]), character(1))
  ymax <- top * (1.10 + .115 * max(1, sum(ann != "")))
  cmp  <- list(c("Control", "PR8"), c("Control", "S33"), c("PR8", "S33"))
  keepi <- which(ann != "")

  g <- ggplot(su, aes(group, M)) +
    geom_col(aes(fill = group), alpha = .6, width = .7, colour = NA) +
    geom_errorbar(aes(ymin = pmax(0, M - SE), ymax = M + SE, colour = group),
                  width = .22, linewidth = .6) +
    geom_jitter(data = r, aes(group, v, colour = group), inherit.aes = FALSE,
                shape = 21, fill = "white", size = dot_size, stroke = dot_stroke,
                width = .07, height = 0) +
    scale_fill_manual(values = COL, guide = "none") +
    scale_colour_manual(values = COL, guide = "none") +
    scale_x_discrete(labels = XLB) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), oob = scales::oob_keep) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_text(size = 9, face = "bold"),
          axis.title.y = element_text(size = 9.8, margin = margin(r = 4),
                                      face = if (bold_ylab) "bold" else "plain"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(8, 5, 3, 5))

  if (length(keepi)) {
    ys <- top * (1.05 + .105 * (seq_along(keepi) - 1))
    g <- g + geom_signif(
                         xmin = vapply(cmp[keepi], function(z) match(z[1], group_levels), numeric(1)),
                         xmax = vapply(cmp[keepi], function(z) match(z[2], group_levels), numeric(1)),
                         annotations = ann[keepi], y_position = ys, manual = FALSE,
                         tip_length = .012, textsize = star_size, vjust = .45,
                         colour = "black")
  }
  list(plot = g, stats = data.frame(
    Control = unname(st$emmean[["Control"]]), PR8 = unname(st$emmean[["PR8"]]),
    S33 = unname(st$emmean[["S33"]]), p_omni = st$p_omni,
    p_C_PR8 = unname(pv[["Control - PR8"]]), p_C_S33 = unname(pv[["Control - S33"]]),
    p_PR8_S33 = unname(pv[["PR8 - S33"]])))
}


pand <- file.path(out_fig, "panels")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
if (save_panels) dir.create(pand, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

safe <- function(x) x %>% str_replace_all("α", "a") %>% str_replace_all("β", "b") %>%
  str_replace_all("γ", "g") %>% str_replace_all("[^A-Za-z0-9._-]+", "_") %>% str_remove("_$")

ST <- list()
for (tis in levels(dat$Tissue)) {
  dt  <- droplevels(filter(dat, Tissue == tis))
  cys <- levels(dt$Cytokine)
  P <- lapply(cys, function(cy) {
    r <- filter(dt, Cytokine == cy)
    o <- bar_panel(r, y_title(cy))
    s <- o$stats; s$panel <- paste(tis, cy, sep = " | ")
    ST[[s$panel]] <<- s[, c("panel", setdiff(names(s), "panel"))]
    if (save_panels)
      save_png(file.path(pand, sprintf("%s_%s.png", safe(tis), safe(cy))),
               o$plot, width = panel_w, height = panel_h, dpi = dpi, bg = "white")
    o$plot
  })
  nr <- ceiling(length(P) / fig_ncol)
  f <- plot_grid(plotlist = P, nrow = nr, ncol = fig_ncol,
                 labels = if (panel_labels) letters[seq_along(P)] else NULL,
                 label_size = 13) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save_png(file.path(out_fig, paste0("FIG_", tis, ".png")), f,
           width = fig_ncol * panel_w, height = nr * panel_h, dpi = dpi, bg = "white")
  message("  + ", out_fig, "/FIG_", tis, ".png  (", length(P), " panel)")
}
if (save_panels) message("  + ", pand, "/*.png")

stats <- bind_rows(ST)
write.csv(stats, file.path(out_tab, "ELISA_bar_figures_stats.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")
message("  + ", out_tab, "/ELISA_bar_figures_stats.csv")

cat("\n== Comparisons with p < 0.05 ==\n")
stats %>%
  pivot_longer(c(p_C_PR8, p_C_S33, p_PR8_S33), names_to = "comparison", values_to = "p") %>%
  filter(p < .05) %>% arrange(p) %>%
  mutate(comparison = recode(comparison, p_C_PR8 = "Naive vs A/PR/8/34",
                             p_C_S33 = "Naive vs A/NWS/33", p_PR8_S33 = "A/PR/8/34 vs A/NWS/33"),
         p = signif(p, 3)) %>%
  select(panel, comparison, p) %>% as.data.frame() %>% print(row.names = FALSE)
message("\nXong.")
