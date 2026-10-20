# ==============================================================================
# 00_Build-Matched-Dataset.R 
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(readxl) })

PROJ <- "."
OUT  <- file.path(PROJ, "CrossPlatform_aSyn", "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TREATS <- c("Control", "V1", "V2", "PolyIC")
TPS    <- c("4h", "8h", "24h")
logit  <- function(p) { p <- pmin(pmax(p, 1e-4), 1 - 1e-4); log(p / (1 - p)) }
key_ok <- function(df) df %>% mutate(treatment = factor(treatment, levels = TREATS),
                                     timepoint = factor(timepoint, levels = TPS),
                                     rep = as.integer(rep)) %>% filter(rep <= 4)

# ---- WB Phase2 (anchor) ------------------------------------------------------
wb <- read_excel(file.path(PROJ, "Western blot/DATA/WB_FORR.xlsx"), sheet = "Phase2_Timecourse") %>%
  filter(group %in% c("C","PO","PR8","S33")) %>%
  transmute(timepoint, rep,
            treatment = recode(group, C="Control", PR8="V1", S33="V2", PO="PolyIC"),
            # naming + formula matched to Phase-1 WB targets (each = log2(band/actin))
            WB_aSyn_25kDa    = log2(as.numeric(syn1.25) / as.numeric(actin)),
            WB_aSyn_multimer = log2(as.numeric(syn1)    / as.numeric(actin)),  # = Phase-1 multimer
            WB_aSyn_pSyn     = log2(as.numeric(psyn)    / as.numeric(actin))) %>% # = Phase-1 phospho-aSyn
  key_ok()

# ---- ELISA -------------------------------------------------------------------
elisa <- read_excel(file.path(PROJ, "Flow-ELISA/DATA/ELISA_ForR.xlsx")) %>%
  filter(group %in% c("C","C_PR8","C_S33","PolyIC")) %>%
  mutate(treatment = recode(group, C="Control", C_PR8="V1", C_S33="V2", PolyIC="PolyIC")) %>%
  transmute(timepoint, rep, treatment,
            across(c(CCL2, IL_1b, IL_6, IL_10, TNFa, IL_12p70),
                   ~ log2(as.numeric(.x)), .names = "ELISA_{.col}")) %>%
  key_ok()

# ---- Flow (FACS, % total live → logit); include NP as viral-load covariate ---
facs <- read_excel(file.path(PROJ, "Flow-ELISA/DATA/FACS_ForR.xlsx")) %>%
  filter(group %in% c("C","C_PR8","C_S33","PolyIC")) %>%
  mutate(treatment = recode(group, C="Control", C_PR8="V1", C_S33="V2", PolyIC="PolyIC")) %>%
  transmute(timepoint, rep, treatment,
            across(c(CD40, CD68, CD86, CD206, HLA_I, HLA_II, NP),
                   ~ logit(as.numeric(.x) / 100), .names = "FACS_{.col}")) %>%
  key_ok()

# ---- qPCR (per-rep FC, 4h/8h only) → log2 ; pivot genes wide -----------------
qmap <- tibble(Group_ID = paste0("G", c(1,2,3,7,8,9,10,14)),
               treatment = rep(c("Control","V1","V2","PolyIC"), 2),
               timepoint = rep(c("4h","8h"), each = 4))
qpcr <- read_csv(file.path(PROJ, "2606_PCR_VirusHMC3/analysis_master/04_per_rep_fc_sdm.csv"),
                 show_col_types = FALSE) %>%
  filter(Group_ID %in% qmap$Group_ID, is.finite(FoldChange), FoldChange > 0,
         Gene %in% c("IL1B","IL6","IRF7","RSAD2","MX1","OAS1","PSMB5","PSMB8")) %>%
  mutate(rep = as.integer(sub("Rep", "", Bio_Rep))) %>%
  left_join(qmap, by = "Group_ID") %>%
  group_by(timepoint, treatment, rep, Gene) %>%        # collapse technical dups
  summarize(v = log2(mean(FoldChange)), .groups = "drop") %>%
  pivot_wider(names_from = Gene, values_from = v, names_prefix = "qPCR_") %>%
  key_ok()

# ---- JOIN (left on WB anchor) ------------------------------------------------
master <- wb %>%
  left_join(elisa, by = c("timepoint","treatment","rep")) %>%
  left_join(facs,  by = c("timepoint","treatment","rep")) %>%
  left_join(qpcr,  by = c("timepoint","treatment","rep")) %>%
  arrange(timepoint, treatment, rep)

write_csv(master, file.path(OUT, "master_matched.csv"))

# ---- diagnostics -------------------------------------------------------------
cat("== master_matched ==\n"); cat("rows:", nrow(master), " cols:", ncol(master), "\n\n")
cat("rows per timepoint × treatment (anchor WB):\n")
print(master %>% count(timepoint, treatment) %>% pivot_wider(names_from=treatment, values_from=n))
cat("\nmarker coverage (non-NA count / ", nrow(master), "):\n", sep="")
cov <- master %>% summarize(across(starts_with(c("WB_","ELISA_","FACS_","qPCR_")), ~sum(!is.na(.x)))) %>%
  pivot_longer(everything(), names_to="marker", values_to="n_obs")
print(as.data.frame(cov), row.names = FALSE)
cat("\nanchor WB_aSyn_25kDa non-NA:", sum(is.finite(master$WB_aSyn_25kDa)), "\n")
message("\n[00_Build-Matched-Dataset] wrote output/master_matched.csv")
