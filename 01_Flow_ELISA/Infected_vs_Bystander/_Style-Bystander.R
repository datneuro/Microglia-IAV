# ==============================================================================
# _Style-Bystander.R  —  Shared house style for the Infected-vs-Bystander module
# Palette identical to the main Flow/ELISA figures (Figures/_Style.R):
#   Naive/Ctrl = gray · V1 (A/PR/8/34) = blue · V2 (A/NWS/33) = red · PolyI:C = green
#   light shade = -LPS · dark shade = +LPS
# Adds NP-gate (State) colours used only in this module.
# Sourced by the analysis scripts (marker sets, sig_star) and by
# Figures/Figure-01_Microglial-Activation.R.
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

GROUP_LV  <- c("C","C_PR8","C_S33","LPS","LPS_PR8","LPS_S33","PolyIC")
GROUP_LAB <- c(C="Naive", C_PR8="A/PR/8/34", C_S33="A/NWS/33",
               LPS="LPS", LPS_PR8="LPS+A/PR/8/34", LPS_S33="LPS+A/NWS/33",
               PolyIC="PolyI:C")
GROUP_COLORS <- c(C="#D5D8DC", C_PR8="#A8C7F0", C_S33="#F2B4B4",
                  LPS="#6E7376", LPS_PR8="#4285F4", LPS_S33="#C0392B",
                  PolyIC="#2ECC71")

# Virus strain (averaged over LPS) — the headline axis for this module
VIRUS_COLORS <- c(None="#6E7376", V1="#4285F4", V2="#C0392B")
VIRUS_LAB    <- c(None="Naive / Ctrl", V1="A/PR/8/34 (V1)", V2="A/NWS/33 (V2)")

# Infection state (NP gate) — used for the cell-autonomous vs paracrine story
STATE_COLORS <- c("NP+" = "#B03A2E",   # infected  = deep red
                  "NP-" = "#7F8C8D",   # bystander = slate gray
                  "Ctrl"= "#BDC3C7")   # naive cells (M3A reference)
STATE_LAB    <- c("NP+" = "Infected (NP+)",
                  "NP-" = "Bystander (NP−)",
                  "Ctrl"= "Naive (Ctrl)")

# Markers measured inside the NP+ / NP- gates (6)
GATE_MARKERS <- c("CD40","CD68","CD86","CD206","HLA-I","HLA-II")
SURFACE_M1   <- c("CD40","CD86","HLA-I","HLA-II")   # pro-inflammatory surface
M2_MARK      <- c("CD206")                          # regulatory
NEUTRAL_MARK <- c("CD68")                           # phagocytic / neutral

MARKER_CLASS <- c(CD40="M1 surface", CD86="M1 surface",
                  `HLA-I`="M1 surface", `HLA-II`="M1 surface",
                  CD206="M2 regulatory", CD68="Phagocytic")

# Diverging heatmap anchors (navy ↔ white ↔ deep red), centered at 0
HM_LOW  <- "#114B9A"
HM_MID  <- "#FFFFFF"
HM_HIGH <- "#C0392B"

TP_LV <- c("4h","8h","24h")
tp_fac <- function(x) factor(paste0(gsub("h$","",x), "h"), levels = TP_LV)

theme_paper <- function(base = 12) {
  theme_classic(base_size = base) +
    theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0.5,
                                   margin = margin(t = 5, b = 4)),
      plot.subtitle = element_text(size = 8, colour = "grey40", hjust = 0.5,
                                   lineheight = 1.2, margin = margin(b = 8)),
      axis.title    = element_text(face = "bold", size = 11),
      axis.text     = element_text(size = 9),
      axis.text.x   = element_text(size = 9, face = "bold"),
      strip.text    = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 8.5)
    )
}

lab_for <- function(codes) unname(GROUP_LAB[codes])

# significance star helper (Holm-adjusted p)
sig_star <- function(p) ifelse(is.na(p), "",
                        ifelse(p < 0.001, "***",
                        ifelse(p < 0.01,  "**",
                        ifelse(p < 0.05,  "*",
                        ifelse(p < 0.10,  ".", "")))))
