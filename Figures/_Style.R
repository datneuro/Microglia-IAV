
GROUP_LV  <- c("C","C_PR8","C_S33","LPS","LPS_PR8","LPS_S33","PolyIC")
GROUP_LV6 <- GROUP_LV[1:6]

GROUP_LAB <- c(C="Naive", C_PR8="A/PR/8/34", C_S33="A/NWS/33",
               LPS="LPS", LPS_PR8="LPS+A/PR/8/34", LPS_S33="LPS+A/NWS/33",
               PolyIC="PolyI:C")

GROUP_COLORS <- c(C="#D5D8DC", C_PR8="#A8C7F0", C_S33="#F2B4B4",
                  LPS="#6E7376", LPS_PR8="#4285F4", LPS_S33="#C0392B",
                  PolyIC="#2ECC71")

VIRUS_COLORS <- c(None="#6E7376", V1="#4285F4", V2="#C0392B")
VIRUS_LAB    <- c(None="Naive / Ctrl", V1="A/PR/8/34 (V1)", V2="A/NWS/33 (V2)")

# Marker classes (composite sets of 01_Flow_ELISA/02_Composite-Indices.R)
SURFACE_M1  <- c("CD40","CD86","HLA-I","HLA-II")
SECRETED_M1 <- c("TNF-a","IL-1b","IL-6","IL-12p70")
M2_MARK     <- c("CD206","IL-10")
OTHER_MARK  <- c("CD68","NP","CCL2")

# Diverging heatmap anchors (navy <-> white <-> deep red), centered at 0
HM_LOW  <- "#114B9A"
HM_MID  <- "#FFFFFF"
HM_HIGH <- "#C0392B"

theme_paper <- function(base = 12) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(
      plot.title     = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5,
                                              margin = ggplot2::margin(t = 5, b = 4)),
      plot.subtitle  = ggplot2::element_text(size = 8, color = "grey40", hjust = 0.5,
                                              lineheight = 1.2, margin = ggplot2::margin(b = 8)),
      axis.title     = ggplot2::element_text(face = "bold", size = 11),
      axis.text      = ggplot2::element_text(size = 9),
      axis.text.x    = ggplot2::element_text(size = 9, face = "bold"),
      strip.text     = ggplot2::element_text(face = "bold", size = 10),
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold", size = 9),
      legend.text      = ggplot2::element_text(size = 8.5)
    )
}

lab_for <- function(codes) unname(GROUP_LAB[codes])
