
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(patchwork)
})

BASE   <- "Flow-ELISA" 
FIGDIR <- Sys.getenv("FIGURE_OUT", unset = file.path("Figures", "output"))
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
source(file.path("Figures", "_Style.R"))

map_grp <- tibble(
  group = GROUP_LV,
  LPS   = c("No LPS","No LPS","No LPS","LPS","LPS","LPS","PolyIC"),
  Virus = c("None","V1","V2","None","V1","V2","PolyIC")
)

facs <- read_excel(file.path(BASE, "DATA", "FACS_ForR.xlsx"), sheet = "FACS_merged") %>%
  mutate(timepoint = factor(timepoint, levels = c("4h","8h","24h")),
         group     = factor(group, levels = GROUP_LV)) %>%
  left_join(map_grp, by = "group")

summ <- function(df, col) {
  df %>% group_by(group, Virus, LPS, timepoint) %>%
    summarise(m  = mean(.data[[col]], na.rm = TRUE),
              se = sd(.data[[col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[col]]))),
              .groups = "drop") %>%
    mutate(Virus = factor(Virus, levels = c("None","V1","V2")),
           LPS   = factor(LPS,   levels = c("No LPS","LPS")))
}

kinetics_panel <- function(df_sum, ycol_lab, ttl, sub) {
  ggplot(df_sum, aes(timepoint, m, color = Virus, fill = Virus, group = Virus)) +
    geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(shape = 21, fill = "white", size = 2.6, stroke = 1.1) +
    facet_wrap(~ LPS) +
    scale_color_manual(values = VIRUS_COLORS,
                       labels = unname(VIRUS_LAB[c("None","V1","V2")]), name = NULL) +
    scale_fill_manual(values = VIRUS_COLORS, guide = "none") +
    labs(title = ttl, subtitle = sub, x = NULL, y = ycol_lab) +
    theme_paper()
}

facs_fac <- filter(facs, LPS != "PolyIC")

p1a <- kinetics_panel(summ(facs_fac, "NP"),
                      "% NP⁺ (of live cells)",
                      "NP⁺ infection kinetics",
                      "Replication readout · V1 (A/PR/8/34) > V2 (A/NWS/33)")

p1b <- kinetics_panel(summ(facs_fac, "Live_cell"),
                      "% Live cells",
                      "Cell viability kinetics",
                      "Viability / cell death across time")

fig1 <- (p1a / p1b) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "IAV replication and viability in HMC3 microglia",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
  ) & theme(legend.position = "bottom")

ggsave(file.path(FIGDIR, "FigureS3_viability.pdf"), fig1, width = 8.5, height = 7.2, device = cairo_pdf)
ggsave(file.path(FIGDIR, "FigureS3_viability.png"), fig1, width = 8.5, height = 7.2, dpi = 400)
message("Saved FigureS3_viability.pdf / .png  → ", FIGDIR)
