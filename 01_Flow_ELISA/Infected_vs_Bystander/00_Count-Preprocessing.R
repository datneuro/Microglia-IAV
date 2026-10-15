# ==============================================================================
# 00_Count-Preprocessing.R  —  Ingest raw FlowJo event COUNTS for the NP+/NP- modules
#  Source: DATA/FACS_Count_ForR.xlsx  (numerators = gated event counts,
#           denominators = parent-gate counts → enables beta-binomial GLMM)
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(purrr)
})
set.seed(2026)

BASE <- "./Flow-ELISA"                                     # data folder (contains DATA/)
DATA <- file.path(BASE, "DATA")
HERE <- file.path(BASE, "R-analysis-figure_VS-bystander")  # module output folder
OUT  <- file.path(HERE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── Load raw counts ───────────────────────────────────────────────────────────
cnt <- read_excel(file.path(DATA, "FACS_Count_ForR.xlsx"), sheet = 1)
cat("Raw count file:", nrow(cnt), "wells x", ncol(cnt), "cols\n")

# ── Design factors (as in 01_Flow_ELISA/00_Data-Preprocessing.R) ─────────────────
add_design <- function(df) {
  df %>% mutate(
    Timepoint = as.integer(str_extract(timepoint, "\\d+")),
    LPS = if_else(group %in% c("LPS","LPS_PR8","LPS_S33"), "Yes", "No"),
    Virus = case_when(
      group %in% c("C_PR8","LPS_PR8") ~ "V1",
      group %in% c("C_S33","LPS_S33") ~ "V2",
      group == "PolyIC"               ~ "PolyIC",
      TRUE                            ~ "None"),
    Bio_Rep = as.integer(rep),
    Group = group)
}
cnt <- add_design(cnt)

markers6 <- c("CD40","CD68","CD86","CD206","HLA_I","HLA_II")
nice <- c(CD40="CD40", CD68="CD68", CD86="CD86", CD206="CD206",
          HLA_I="HLA-I", HLA_II="HLA-II")

# helper: clamp numerator to its gate denominator (guards FlowJo rounding)
clamp <- function(num, den) pmin(num, den)

# ══════════════════════════════════════════════════════════════════════════════
# 1. GATED long table  (NP+ and NP- compartments) — infected wells only
#    one row per well × marker × State{NP+, NP-}
# ══════════════════════════════════════════════════════════════════════════════
infected_groups <- c("C_PR8","C_S33","LPS_PR8","LPS_S33")

gated <- map_dfr(markers6, function(m) {
  npp <- paste0(m, "/NP+"); npn <- paste0(m, "/NP-")
  cnt %>%
    filter(Group %in% infected_groups) %>%
    transmute(
      sample_id, Group, Timepoint, LPS, Virus, Bio_Rep,
      Marker = nice[[m]],
      NPp_pos = clamp(.data[[npp]], NP),     NPp_den = NP,
      NPn_pos = clamp(.data[[npn]], NP_neg), NPn_den = NP_neg)
}) %>%
  pivot_longer(c(NPp_pos, NPp_den, NPn_pos, NPn_den),
               names_to = c("State","quant"), names_sep = "_") %>%
  pivot_wider(names_from = quant, values_from = value) %>%
  mutate(State = recode(State, NPp = "NP+", NPn = "NP-"),
         State = factor(State, levels = c("NP-","NP+")),   # NP- = reference
         frac  = pos / den,
         pct   = 100 * frac)

write.csv(gated, file.path(OUT, "counts_gated_long.csv"), row.names = FALSE)
cat("counts_gated_long.csv :", nrow(gated), "rows (",
    n_distinct(gated$sample_id), "wells x", length(markers6), "markers x 2 states )\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. OVERALL / paracrine long table 
#    Naive & LPS wells contribute total-live marker fraction (Virus = None);
#    infected wells contribute their NP- (bystander) fraction.
# ══════════════════════════════════════════════════════════════════════════════
ref_overall <- map_dfr(markers6, function(m) {
  cnt %>%
    filter(Group %in% c("C","LPS")) %>%
    transmute(sample_id, Group, Timepoint, LPS, Virus, Bio_Rep,
              Marker = nice[[m]], Compartment = "Ctrl",
              pos = clamp(.data[[m]], Live_cell), den = Live_cell)
})
bystander <- map_dfr(markers6, function(m) {
  npn <- paste0(m, "/NP-")
  cnt %>%
    filter(Group %in% infected_groups) %>%
    transmute(sample_id, Group, Timepoint, LPS, Virus, Bio_Rep,
              Marker = nice[[m]], Compartment = "Bystander",
              pos = clamp(.data[[npn]], NP_neg), den = NP_neg)
})
overall <- bind_rows(ref_overall, bystander) %>%
  mutate(Virus = factor(Virus, levels = c("None","V1","V2")),
         frac = pos / den, pct = 100 * frac)

write.csv(overall, file.path(OUT, "counts_overall_long.csv"), row.names = FALSE)
cat("counts_overall_long.csv:", nrow(overall), "rows (Ctrl + bystander reference set)\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. INFECTION KINETICS per well
# ══════════════════════════════════════════════════════════════════════════════
kin <- cnt %>%
  transmute(sample_id, Group, Timepoint, LPS, Virus, Bio_Rep,
            Live_cell, Single_cell,
            NP_pos = NP, NP_den = Live_cell,
            pctNP  = 100 * NP / Live_cell,
            viability = 100 * Live_cell / Single_cell,
            # whole-well M1-surface activation (% of live cells) for M4
            surfM1_pos = clamp(CD40, Live_cell) + clamp(CD86, Live_cell) +
                         clamp(HLA_I, Live_cell) + clamp(HLA_II, Live_cell),
            surfM1_den = 4 * Live_cell,
            surfM1_pct = 100 * surfM1_pos / surfM1_den)
write.csv(kin, file.path(OUT, "infection_kinetics.csv"), row.names = FALSE)
cat("infection_kinetics.csv :", nrow(kin), "wells\n")

# ── Diagnostics ───────────────────────────────────────────────────────────────
cat("\n── Gate-size diagnostics (infected wells) ──\n")
gated %>% filter(State == "NP+") %>%
  group_by(Group, Timepoint) %>%
  summarise(min_NPgate = min(den), med_NPgate = median(den), .groups="drop") %>%
  as.data.frame() %>% print()
cat("\nWells with NP+ gate < 100 events (interpret with caution):",
    sum(unique(gated[gated$State=="NP+", c("sample_id","den")])$den < 100), "\n")
cat("\n═══ 00_Count-Preprocessing COMPLETE ═══\n")
