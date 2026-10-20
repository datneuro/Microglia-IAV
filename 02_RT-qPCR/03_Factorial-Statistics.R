# ==============================================================================
# 03_Factorial-Statistics.R  —  RT-qPCR factorial statistics
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(car)
  library(emmeans)
  library(ggsignif)
  library(patchwork)
  library(openxlsx)
  library(grid)
})

# ==============================================================================
# §1  GLOBAL CONFIGURATION
# ==============================================================================

PROJ_DIR <- "./2606_PCR_VirusHMC3"
OUT_DIR  <- file.path(PROJ_DIR, "analysis_master")

INPUT_PER_REP <- file.path(OUT_DIR, "04_per_rep_fc_sdm.csv")
FIG_OUT       <- file.path(OUT_DIR, "04_master_expression_log2_violin")
XLSX_OUT      <- file.path(OUT_DIR, "04_master_qpcr_stats_report.xlsx")

ADJUST_METHOD   <- "holm"
SHOW_ONLY_SIG   <- TRUE
REFERENCE_GENES <- c("GAPDH", "RPS11")

STRAIN_V1 <- "A/PR/8/34"
STRAIN_V2 <- "A/NWS/33"

# ==============================================================================
# §2  GROUP DESIGN MAP
# ==============================================================================

GROUP_MAP <- tribble(
  ~Group_ID, ~Timepoint, ~LPS,      ~Virus,  ~ShortLabel,                    ~Label,                         ~in_factorial,
  "G1",  "4h",  "No LPS",  "None", "Ctrl",                            "4h-Ctrl",                            TRUE,
  "G2",  "4h",  "No LPS",  "V1",   STRAIN_V1,                         paste0("4h-", STRAIN_V1),             TRUE,
  "G3",  "4h",  "No LPS",  "V2",   STRAIN_V2,                         paste0("4h-", STRAIN_V2),             TRUE,
  "G4",  "4h",  "LPS",     "None", "LPS",                             "4h-LPS",                             TRUE,
  "G5",  "4h",  "LPS",     "V1",   paste0("LPS+", STRAIN_V1),         paste0("4h-LPS+", STRAIN_V1),         TRUE,
  "G6",  "4h",  "LPS",     "V2",   paste0("LPS+", STRAIN_V2),         paste0("4h-LPS+", STRAIN_V2),         TRUE,
  "G7",  "4h",  NA,        NA,     "PolyI:C",                         "4h-PolyI:C",                         FALSE,
  "G8",  "8h",  "No LPS",  "None", "Ctrl",                            "8h-Ctrl",                            TRUE,
  "G9",  "8h",  "No LPS",  "V1",   STRAIN_V1,                         paste0("8h-", STRAIN_V1),             TRUE,
  "G10", "8h",  "No LPS",  "V2",   STRAIN_V2,                         paste0("8h-", STRAIN_V2),             TRUE,
  "G11", "8h",  "LPS",     "None", "LPS",                             "8h-LPS",                             TRUE,
  "G12", "8h",  "LPS",     "V1",   paste0("LPS+", STRAIN_V1),         paste0("8h-LPS+", STRAIN_V1),         TRUE,
  "G13", "8h",  "LPS",     "V2",   paste0("LPS+", STRAIN_V2),         paste0("8h-LPS+", STRAIN_V2),         TRUE,
  "G14", "8h",  NA,        NA,     "PolyI:C",                         "8h-PolyI:C",                         FALSE
)

GRP_LEVELS <- paste0("G", 1:14)
SHORT_MAP  <- setNames(GROUP_MAP$ShortLabel, GROUP_MAP$Group_ID)

# ==============================================================================
# §3  4-FAMILY COLOR PALETTE + SPLIT-KEY LEGEND GLYPH
# ==============================================================================
# Light shade = −LPS state.  Dark shade = +LPS state.
# PolyI:C: solid green.

GROUP_COLORS <- c(
  "G1"  = "#D5D8DC",  "G8"  = "#D5D8DC",   # light gray  — Ctrl
  "G2"  = "#A8C7F0",  "G9"  = "#A8C7F0",   # light blue  — V1
  "G3"  = "#F2B4B4",  "G10" = "#F2B4B4",   # light red   — V2
  "G4"  = "#6E7376",  "G11" = "#6E7376",   # dark gray   — LPS
  "G5"  = "#4285F4",  "G12" = "#4285F4",   # dark blue   — LPS+V1
  "G6"  = "#C0392B",  "G13" = "#C0392B",   # dark red    — LPS+V2
  "G7"  = "#2ECC71",  "G14" = "#2ECC71"    # stable green — PolyI:C
)

# 4-item legend: one per color family; split key shows light (−LPS) / dark (+LPS)
LEGEND_BREAKS <- c("G1",  "G2",  "G3",  "G7")
LEGEND_LABELS <- c("Naive", STRAIN_V1, STRAIN_V2, "PolyI:C")

# Split-triangle legend key: upper-left triangle = light shade, lower-right = dark shade
# Note: key glyph functions use plain [0,1] normalised coords, NOT unit() objects.
draw_key_split_grp <- function(data, params, size) {
  fc <- data$fill
  pair <- if (is.na(fc))                               c("grey90",  "grey60") else
          if (fc %in% c("#D5D8DC", "#6E7376"))         c("#D5D8DC", "#6E7376") else
          if (fc %in% c("#A8C7F0", "#4285F4"))         c("#A8C7F0", "#4285F4") else
          if (fc %in% c("#F2B4B4", "#C0392B"))         c("#F2B4B4", "#C0392B") else
                                                       c(fc, fc)
  grid::grobTree(
    grid::polygonGrob(x = c(0, 0, 1), y = c(0, 1, 1),
                      gp = grid::gpar(fill = pair[1], col = NA)),
    grid::polygonGrob(x = c(0, 1, 1), y = c(0, 0, 1),
                      gp = grid::gpar(fill = pair[2], col = NA)),
    grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey50", lwd = 0.7))
  )
}

# Per-gene minimum y-floor and y-ceiling for genes with small effects or bottom clips.
# Values are added symmetrically (or asymmetrically where commented) to the natural range.
Y_GENE_FLOOR <- c(
  "IL1B"  = -6.0,
  "IL6"   = -3.0,
  "IRF7"  = -4.0,
  "MX1"   = -4.0,
  "OAS1"  = -4.0,
  "PSMB5" = -4.0,
  "PSMB8" = -4.0,
  "RSAD2" = -4.0
)
Y_GENE_CEIL <- c(
  "IL1B"  = 6.0,
  "IL6"   = 3.0,
  "IRF7"  = 12.0,
  "MX1"   = 12.0,
  "OAS1"  = 12.0,
  "PSMB5" = 5.0,
  "PSMB8" = 6.0,
  "RSAD2" = 14.0
)

# ==============================================================================
# §4  EXPANDED CONTRAST SPECIFICATION (CORRECTED EMMEANS GRID ORDER)
# ==============================================================================
# True emmeans grid order for ~ LPS * Virus | Timepoint (expand.grid alignment):
#   pos 1: (No LPS, None) [Ctrl]       pos 2: (LPS, None)    [LPS]
#   pos 3: (No LPS, V1)   [A/PR/8/34]  pos 4: (LPS, V1)      [LPS+A/PR/8/34]
#   pos 5: (No LPS, V2)   [A/NWS/33]   pos 6: (LPS, V2)      [LPS+A/NWS/33]

FACTORIAL_CONTRASTS <- list(
  list(
    label  = "LPS baseline: Ctrl vs LPS",
    family = "LPS Baseline",
    vec    = c(-1, +1,  0,  0,  0,  0), # pos 2 - pos 1
    grp_4h = c("G1", "G4"),  grp_8h = c("G8",  "G11")
  ),
  list(
    label  = paste0("Infection (−LPS): Ctrl vs ", STRAIN_V1),
    family = "Infection",
    vec    = c(-1,  0, +1,  0,  0,  0), # pos 3 - pos 1
    grp_4h = c("G1", "G2"),  grp_8h = c("G8",  "G9")
  ),
  list(
    label  = paste0("Infection (−LPS): Ctrl vs ", STRAIN_V2),
    family = "Infection",
    vec    = c(-1,  0,  0,  0, +1,  0), # pos 5 - pos 1
    grp_4h = c("G1", "G3"),  grp_8h = c("G8",  "G10")
  ),
  list(
    label  = paste0("Infection (−LPS): ", STRAIN_V1, " vs ", STRAIN_V2),
    family = "Infection",
    vec    = c( 0,  0, -1,  0, +1,  0), # pos 5 - pos 3
    grp_4h = c("G2", "G3"),  grp_8h = c("G9",  "G10")
  ),
  list(
    label  = paste0("Co-stimulation (+LPS): LPS vs LPS+", STRAIN_V1),
    family = "Co-stimulation",
    vec    = c( 0, -1,  0, +1,  0,  0), # pos 4 - pos 2
    grp_4h = c("G4", "G5"),  grp_8h = c("G11", "G12")
  ),
  list(
    label  = paste0("Co-stimulation (+LPS): LPS vs LPS+", STRAIN_V2),
    family = "Co-stimulation",
    vec    = c( 0, -1,  0,  0,  0, +1), # pos 6 - pos 2
    grp_4h = c("G4", "G6"),  grp_8h = c("G11", "G13")
  ),
  list(
    label  = paste0("Co-stimulation (+LPS): LPS+", STRAIN_V1, " vs LPS+", STRAIN_V2),
    family = "Co-stimulation",
    vec    = c( 0,  0,  0, -1,  0, +1), # pos 6 - pos 4
    grp_4h = c("G5", "G6"),  grp_8h = c("G12", "G13")
  )
)

POLYIC_CONTRASTS <- list(
  list(
    label  = "PolyI:C cross-ref: Ctrl vs PolyI:C",
    family = "PolyI:C Reference",
    grp_4h = c("G1",  "G7"),  grp_8h = c("G8",  "G14")
  ),
  list(
    label  = "PolyI:C cross-ref: LPS vs PolyI:C",
    family = "PolyI:C Reference",
    grp_4h = c("G4",  "G7"),  grp_8h = c("G11", "G14")
  ),
  list(
    label  = paste0("PolyI:C cross-ref: ", STRAIN_V1, " vs PolyI:C"),
    family = "PolyI:C Reference",
    grp_4h = c("G2",  "G7"),  grp_8h = c("G9",  "G14")
  ),
  list(
    label  = paste0("PolyI:C cross-ref: ", STRAIN_V2, " vs PolyI:C"),
    family = "PolyI:C Reference",
    grp_4h = c("G3",  "G7"),  grp_8h = c("G10", "G14")
  )
)

# Pre-build lookup for joining emmeans results to group pairs
FACTORIAL_SPEC_TBL <- tibble(
  contrast        = sapply(FACTORIAL_CONTRASTS, `[[`, "label"),
  Contrast_Family = sapply(FACTORIAL_CONTRASTS, `[[`, "family"),
  grp_4h_ref      = sapply(FACTORIAL_CONTRASTS, function(x) x$grp_4h[1]),
  grp_4h_test     = sapply(FACTORIAL_CONTRASTS, function(x) x$grp_4h[2]),
  grp_8h_ref      = sapply(FACTORIAL_CONTRASTS, function(x) x$grp_8h[1]),
  grp_8h_test     = sapply(FACTORIAL_CONTRASTS, function(x) x$grp_8h[2])
)

# ==============================================================================
# §5  UTILITY FUNCTIONS
# ==============================================================================

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

fmt_p <- function(p) {
  if (is.na(p))  return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_f    <- function(f)  ifelse(is.na(f), "NA", sprintf("%.2f", f))
fmt_eta2 <- function(x)  ifelse(is.na(x), "NA", sprintf("%.3f", x))

# ==============================================================================
# §6  DATA INGESTION 
# ==============================================================================
message("\n=== STAGE 1: DATA INGESTION ===")

if (!file.exists(INPUT_PER_REP))
  stop("Per-rep CSV not found: ", INPUT_PER_REP,
       "\nRun 02_Preprocessing-Pipeline.R first.")

raw <- read_csv(INPUT_PER_REP, show_col_types = FALSE)

# Column name normalisation
if ("Gene" %in% colnames(raw) && !"Target_Name" %in% colnames(raw))
  raw <- rename(raw, Target_Name = Gene)
if ("Target Name" %in% colnames(raw) && !"Target_Name" %in% colnames(raw))
  raw <- rename(raw, Target_Name = `Target Name`)
if ("Method" %in% colnames(raw))
  raw <- filter(raw, grepl("SDM", Method))
if (!"Timepoint" %in% colnames(raw))
  raw <- mutate(raw, Timepoint = if_else(Group_ID %in% paste0("G", 1:7), "4h", "8h"))

# No additional Exclude.amplicon filter
df_clean <- raw %>%
  filter(!is.na(FoldChange), is.finite(FoldChange), FoldChange > 0,
         !Target_Name %in% REFERENCE_GENES) %>%
  left_join(
    GROUP_MAP %>% select(Group_ID, LPS, Virus, ShortLabel, Label, in_factorial),
    by = "Group_ID"
  ) %>%
  mutate(
    Group_ID  = factor(Group_ID,  levels = GRP_LEVELS),
    Timepoint = factor(Timepoint, levels = c("4h", "8h")),
    LPS       = factor(LPS,       levels = c("No LPS", "LPS")),
    Virus     = factor(Virus,     levels = c("None", "V1", "V2")),
    Bio_Rep   = factor(Bio_Rep)
  )

df_factorial <- filter(df_clean, in_factorial == TRUE)
target_genes <- c("IL1B", "IL6", "IRF7", "RSAD2", "MX1", "OAS1", "PSMB5", "PSMB8")
n_genes      <- length(target_genes)
n_bioreps    <- n_distinct(df_clean$Bio_Rep)

message("  Target genes  : ", paste(target_genes, collapse = ", "))
message("  Bio replicates: ", n_bioreps)
message("  Total rows    : ", nrow(df_clean), "  (no Exclude.amplicon filter)")
message("  Factorial rows: ", nrow(df_factorial))

# ==============================================================================
# §7  GROUP SUMMARY STATS
# ==============================================================================

summary_df <- df_clean %>%
  group_by(Group_ID, Target_Name) %>%
  summarize(
    N_Bio_Reps  = n_distinct(Bio_Rep),
    N_obs       = n(),
    Mean_FC     = mean(FoldChange),
    SEM_FC      = sd(FoldChange) / sqrt(n()),
    Mean_log2FC = mean(log2(FoldChange)),
    SD_log2FC   = sd(log2(FoldChange)),
    SEM_log2FC  = sd(log2(FoldChange)) / sqrt(n()),
    .groups     = "drop"
  ) %>%
  left_join(
    GROUP_MAP %>% select(Group_ID, ShortLabel, Label, Timepoint, LPS, Virus, in_factorial),
    by = "Group_ID"
  ) %>%
  mutate(Group_ID = factor(Group_ID, levels = GRP_LEVELS))

# ==============================================================================
# §8  RCBD FIXED-BLOCK LINEAR MODEL
# ==============================================================================
message("\n=== STAGE 2: RCBD LINEAR MODEL (log2 scale) ===")

old_contrasts <- getOption("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))

omnibus_list <- list()
model_store  <- list()

for (gene in target_genes) {
  message("  Fitting: ", gene)
  gdf <- filter(df_factorial, Target_Name == gene,
                !is.na(FoldChange), FoldChange > 0)
  if (nrow(gdf) < 6) { message("    SKIP: n=", nrow(gdf)); next }

  fit <- tryCatch(
    lm(log2(FoldChange) ~ Timepoint * LPS * Virus + Bio_Rep, data = gdf),
    error = function(e) { message("    lm failed: ", e$message); NULL }
  )
  if (is.null(fit)) next
  model_store[[gene]] <- list(fit = fit, data = gdf)

  aov_tab <- tryCatch(
    as.data.frame(car::Anova(fit, type = 3)),
    error = function(e) { message("    Anova error: ", e$message); NULL }
  )
  if (is.null(aov_tab)) next

  aov_tab$Term <- rownames(aov_tab)
  rownames(aov_tab) <- NULL

  # Partial eta-squared = SS_term / (SS_term + SS_residual)
  ss_col   <- intersect(c("Sum Sq", "Sum.Sq"), colnames(aov_tab))[1]
  ss_resid <- if (!is.na(ss_col)) {
    r <- aov_tab[[ss_col]][aov_tab$Term == "Residuals"]
    if (length(r) == 1 && !is.na(r) && r > 0) r else NA_real_
  } else NA_real_

  aov_tab$eta2p <- if (!is.na(ss_col) && !is.na(ss_resid))
    aov_tab[[ss_col]] / (aov_tab[[ss_col]] + ss_resid)
  else NA_real_

  p_col <- intersect(c("Pr(>F)", "Pr(>Chisq)"), colnames(aov_tab))[1]
  if (!is.na(p_col))
    aov_tab$Significance <- sapply(aov_tab[[p_col]], p_to_stars)

  omnibus_list[[gene]] <- mutate(aov_tab, Gene = gene, .before = 1)
}

options(contrasts = old_contrasts)

omnibus_df <- bind_rows(omnibus_list) %>%
  filter(!Term %in% c("(Intercept)", "Intercept"))

message("  Omnibus table: ", nrow(omnibus_df), " rows, ",
        n_distinct(omnibus_df$Gene), " genes.")

# ==============================================================================
# §9  EXPANDED PLANNED CONTRASTS  +  GLOBAL HOLM PER GENE
# ==============================================================================
message("\n=== STAGE 3: EXPANDED PLANNED CONTRASTS ===")

CMAT <- setNames(
  lapply(FACTORIAL_CONTRASTS, `[[`, "vec"),
  sapply(FACTORIAL_CONTRASTS, `[[`, "label")
)

all_contrasts_list <- list()

for (gene in target_genes) {
  if (is.null(model_store[[gene]])) next
  fit <- model_store[[gene]]$fit
  message("  Contrasts: ", gene)

  # ── Factorial contrasts via emmeans ────────────────────────────────────────
  emm_lv <- tryCatch(
    emmeans(fit, ~ LPS * Virus | Timepoint),
    error = function(e) { message("    emmeans failed: ", e$message); NULL }
  )

  fac_rows <- NULL
  if (!is.null(emm_lv)) {
    ctr_raw <- tryCatch(
      as.data.frame(summary(contrast(emm_lv, method = CMAT, adjust = "none"))),
      error = function(e) { message("    contrast() failed: ", e$message); NULL }
    )
    if (!is.null(ctr_raw) && nrow(ctr_raw) > 0) {
      fac_rows <- ctr_raw %>%
        left_join(FACTORIAL_SPEC_TBL, by = "contrast") %>%
        mutate(
          Gene           = gene,
          Contrast_Label = contrast,
          Group_ref      = if_else(as.character(Timepoint) == "4h", grp_4h_ref,  grp_8h_ref),
          Group_test     = if_else(as.character(Timepoint) == "4h", grp_4h_test, grp_8h_test),
          P_raw          = p.value
        ) %>%
        select(Gene, Timepoint, Contrast_Label, Contrast_Family,
               Group_ref, Group_test,
               estimate, SE, df, t_stat = t.ratio, P_raw)
    }
  }

  # ── PolyI:C cross-comparisons via Welch t-test (log2 scale) ───────────────
  pic_rows_list <- list()
  for (tp in c("4h", "8h")) {
    for (ps in POLYIC_CONTRASTS) {
      grp_pair <- if (tp == "4h") ps$grp_4h else ps$grp_8h
      g_ref <- grp_pair[1];  g_tst <- grp_pair[2]

      v_ref <- df_clean %>%
        filter(Target_Name == gene, Group_ID == g_ref, FoldChange > 0) %>%
        pull(FoldChange) %>% log2()
      v_tst <- df_clean %>%
        filter(Target_Name == gene, Group_ID == g_tst, FoldChange > 0) %>%
        pull(FoldChange) %>% log2()

      tt <- if (length(v_ref) >= 2 && length(v_tst) >= 2)
        tryCatch(t.test(v_tst, v_ref, var.equal = FALSE), error = function(e) NULL)
      else NULL

      pic_rows_list[[length(pic_rows_list) + 1]] <- tibble(
        Gene           = gene,
        Timepoint      = tp,
        Contrast_Label = ps$label,
        Contrast_Family = ps$family,
        Group_ref      = g_ref,
        Group_test     = g_tst,
        estimate = if (!is.null(tt)) as.numeric(tt$estimate[1] - tt$estimate[2]) else NA_real_,
        SE             = NA_real_,
        df             = if (!is.null(tt)) as.numeric(tt$parameter) else NA_real_,
        t_stat         = if (!is.null(tt)) as.numeric(tt$statistic)  else NA_real_,
        P_raw          = if (!is.null(tt)) tt$p.value               else NA_real_
      )
    }
  }
  pic_rows <- bind_rows(pic_rows_list)

  # ── Pool + global Holm per gene ────────────────────────────────────────────
  gene_rows <- bind_rows(fac_rows, pic_rows)
  if (nrow(gene_rows) == 0) next
  gene_rows$P_adjusted_Holm <- p.adjust(gene_rows$P_raw, method = ADJUST_METHOD)
  gene_rows$Significance    <- sapply(gene_rows$P_adjusted_Holm, p_to_stars)
  all_contrasts_list[[gene]] <- gene_rows
}

contrasts_df <- bind_rows(all_contrasts_list)
n_sig <- sum(!is.na(contrasts_df$P_adjusted_Holm) & contrasts_df$P_adjusted_Holm < 0.05,
             na.rm = TRUE)
message("  Total contrast rows : ", nrow(contrasts_df))
message("  Holm-significant    : ", n_sig)

# ==============================================================================
# §10  SIGNIFICANCE BRACKET BUILDER
# ==============================================================================

grp_order <- setNames(seq_along(GRP_LEVELS), GRP_LEVELS)

build_signif_annot <- function(cdf, only_sig = SHOW_ONLY_SIG) {
  df <- cdf
  if (only_sig)
    df <- filter(df, !is.na(P_adjusted_Holm), P_adjusted_Holm < 0.05,
                 !Significance %in% c("ns", ""))
  if (nrow(df) == 0) return(tibble())

  df %>%
    select(Target_Name = Gene,
           Timepoint,
           x    = Group_ref,
           xend = Group_test,
           label = Significance,
           Contrast_Family) %>%
    mutate(span = grp_order[xend] - grp_order[x]) %>%
    group_by(Target_Name, Timepoint) %>%          # rank within each gene × timepoint
    arrange(abs(span), .by_group = TRUE) %>%
    mutate(bracket_rank = row_number()) %>%
    ungroup() %>%
    select(-span)
}

signif_annot_base <- build_signif_annot(contrasts_df)
message("  Significant brackets: ", nrow(signif_annot_base))

# ==============================================================================
# §11  OMNIBUS ANNOTATION BUILDER
#      Full term names + partial eta-squared
#      Format: "Term: F(df1, df2) = X.XX, p = X.XXX, η²p = X.XXX ***"
# ==============================================================================

get_omnibus_annot <- function(gene, omnibus_df) {
  if (nrow(omnibus_df) == 0) return("")
  gtbl <- filter(omnibus_df, Gene == gene)
  if (nrow(gtbl) == 0) return("")

  f_col  <- intersect(c("F value", "F"), colnames(gtbl))[1]
  p_col  <- intersect(c("Pr(>F)", "Pr(>Chisq)"), colnames(gtbl))[1]
  if (is.na(f_col) || is.na(p_col)) return("")

  df2_resid <- gtbl$Df[gtbl$Term == "Residuals"]
  df2_s     <- if (length(df2_resid) == 1 && !is.na(df2_resid))
    as.character(as.integer(df2_resid)) else "?"

  fmt_row <- function(term_re, full_name) {
    row <- gtbl %>% filter(grepl(term_re, Term, perl = TRUE)) %>% slice(1)
    if (nrow(row) == 0) return(NULL)
    f_val <- row[[f_col]];  p_val <- row[[p_col]]
    df1_s <- if ("Df" %in% colnames(row) && !is.na(row$Df[1]))
      as.character(as.integer(row$Df[1])) else "?"
    eta_v <- if ("eta2p" %in% colnames(row)) row$eta2p[1] else NA_real_
    stars <- p_to_stars(p_val)
    sfx   <- if (!stars %in% c("ns", "")) paste0(" ", stars) else ""
    sprintf("%s: F(%s,%s)=%s, p=%s, η²p=%s%s",
            full_name, df1_s, df2_s, fmt_f(f_val), fmt_p(p_val), fmt_eta2(eta_v), sfx)
  }

  lines <- Filter(Negate(is.null), list(
    fmt_row("Timepoint.*LPS.*Virus|Virus.*LPS.*Timepoint", "Timepoint:LPS:Virus"),
    fmt_row("^LPS:Virus$|^Virus:LPS$",                     "LPS:Virus"),
    fmt_row("^Timepoint:LPS$|^LPS:Timepoint$",             "Timepoint:LPS")
  ))

  paste(lines, collapse = "\n")
}

# ==============================================================================
# §12  VIOLIN-BOXPLOT PLOT BUILDER  (log2 scale)
# ==============================================================================
# Visual design:
#   • Violin: trim=FALSE, adjust=2.0 (smooth over N=3 outlier bumps), wider (width=1.0)
#   • Boxplot: group-colored border, wider (width=0.30), dominant layer
#   • Jitter: small (size=1.8) and transparent (alpha=0.50), subordinate
#   • Timeline "← 4 hours →" / "← 8 hours →" at top of each panel (bold)
#   • No text labels under x-axis for timepoints; only +LPS and PolyI:C bars
#   • Omnibus ANOVA stats in the subtitle (below gene name, above plot area)
#   • No legend in individual panels (legend lives in 9th patchwork slot)

build_gene_plot <- function(gene) {

  jitter_g <- df_clean %>%
    filter(Target_Name == gene, !is.na(FoldChange), FoldChange > 0) %>%
    mutate(.y_val = log2(FoldChange))

  if (nrow(jitter_g) == 0) return(ggplot() + labs(title = gene))

  # ── Y-range: minimum 5 units; always show baseline; expand for no-sig genes ──
  y_vals    <- jitter_g$.y_val
  y_raw_min <- min(y_vals, na.rm = TRUE)
  y_raw_max <- max(y_vals, na.rm = TRUE)
  y_span    <- max(y_raw_max - y_raw_min, 0.5)

  floor_val <- { tmp <- Y_GENE_FLOOR[gene]; if (!is.na(tmp)) tmp else -0.5 }
  ceil_val  <- { tmp <- Y_GENE_CEIL[gene];  if (!is.na(tmp)) tmp else  0.5 }
  y_bot0 <- min(y_raw_min - y_span * 0.28, floor_val)
  y_top0 <- max(y_raw_max + y_span * 0.28, ceil_val)

  # Enforce minimum total range of 5 log2 units
  if (y_top0 - y_bot0 < 5.0) {
    y_mid  <- (y_top0 + y_bot0) / 2
    y_bot0 <- y_mid - 2.5
    y_top0 <- y_mid + 2.5
  }

  # ── Significance brackets: smaller step to reduce crowding ────────────────
  sig_base <- filter(signif_annot_base, Target_Name == gene)
  if (nrow(sig_base) > 0) {
    brk_base  <- y_raw_max + max(y_span * 0.22, 2.5)
    n_brk_max <- max(sig_base$bracket_rank)
    # Distribute brackets evenly across available space up to y_top0;
    # minimum gap 0.40 log2 units keeps labels readable.
    brk_step  <- max((y_top0 - brk_base) / (n_brk_max + 1), 0.40)
    sig_g <- sig_base %>%
      mutate(y_position = brk_base + brk_step * bracket_rank)
    y_top <- max(y_top0, max(sig_g$y_position) + brk_step)
  } else {
    sig_g <- sig_base
    y_top <- y_top0
  }
  y_bot <- y_bot0

  # ── Bottom stratification: fixed physical offset (scales proportionally) ──
  rng     <- y_top - y_bot
  bar_off <- max(rng * 0.060, 0.28)   # +LPS / PolyI:C indicator bar
  txt_off <- max(rng * 0.105, 0.50)   # annotation text below bar
  bar_y   <- y_bot - bar_off
  txt_y   <- y_bot - txt_off

  # ── Omnibus ANOVA text in subtitle ────────────────────────────────────────
  annot_txt <- get_omnibus_annot(gene, omnibus_df)

  ggplot(jitter_g, aes(x = Group_ID, y = .y_val)) +

    # ── Block shading ──────────────────────────────────────────────────────
    annotate("rect", xmin = 0.45, xmax = 7.5,   ymin = -Inf, ymax = Inf,
             fill = "#3B82F6", alpha = 0.04) +
    annotate("rect", xmin = 7.5,  xmax = 14.55, ymin = -Inf, ymax = Inf,
             fill = "#EF4444", alpha = 0.04) +

    # ── 4h | 8h separator ─────────────────────────────────────────────────
    geom_vline(xintercept = 7.5, linetype = "dotted",
               color = "grey35", linewidth = 0.45) +

    # ── Baseline log2 = 0 ─────────────────────────────────────────────────
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.35) +

    # ── Timeline at top of panel (full "hours", bold) ──────────────────────
    annotate("text", x = 4.0,  y = Inf,
             label = "← 4 hours →",
             vjust = -0.35, size = 3.5, fontface = "bold", color = "grey25") +
    annotate("text", x = 11.0, y = Inf,
             label = "← 8 hours →",
             vjust = -0.35, size = 3.5, fontface = "bold", color = "grey25") +

    # ── Layer 1: Violin — no trim; adjust=4 merges isolated N=3 bumps into one body ─
    geom_violin(aes(fill = Group_ID),
                alpha = 0.40, color = NA,
                scale = "width", width = 1.0, trim = FALSE, adjust = 4.0) +

    # ── Layer 2: Boxplot — group-colored border ─────────────────────────────
    geom_boxplot(aes(fill = Group_ID, color = Group_ID),
                 alpha = 0.55, width = 0.30,
                 outlier.shape = NA, lwd = 0.55) +

    # ── Layer 3: Jitter — tight width keeps replicates inside violin body ──
    geom_jitter(aes(color = Group_ID),
                shape = 21, fill = "white",
                size = 2.8, stroke = 0.70, width = 0.05, alpha = 0.85) +

    # ── Significance brackets ──────────────────────────────────────────────
    { if (nrow(sig_g) > 0)
        geom_signif(
          data        = sig_g,
          aes(xmin        = x,
              xmax        = xend,
              annotations = label,
              y_position  = y_position),
          manual      = TRUE,
          tip_length  = 0.01,
          textsize    = 3.0,
          color       = "grey20",
          inherit.aes = FALSE
        )
      else list()
    } +

    # ── Bottom stratification: +LPS indicator bars ─────────────────────────
    # (annotate() draws exactly one element — avoids N-row duplication from geom_*)
    # 4h +LPS (G4–G6)
    annotate("segment", x = 3.6,  xend = 6.4,  y = bar_y, yend = bar_y,
             color = "#1A1A1A", linewidth = 1.1) +
    annotate("text",   x = 5.0,  y = txt_y, label = "+ LPS",
             size = 2.7, fontface = "bold", color = "#1A1A1A", vjust = 1.0) +
    # 8h +LPS (G11–G13)
    annotate("segment", x = 10.6, xend = 13.4, y = bar_y, yend = bar_y,
             color = "#1A1A1A", linewidth = 1.1) +
    annotate("text",   x = 12.0, y = txt_y, label = "+ LPS",
             size = 2.7, fontface = "bold", color = "#1A1A1A", vjust = 1.0) +

    # ── Bottom stratification: PolyI:C indicator bars ──────────────────────
    # 4h PolyI:C (G7)
    annotate("segment", x = 6.6,  xend = 7.4,  y = bar_y, yend = bar_y,
             color = "#2ECC71", linewidth = 1.1) +
    annotate("text",   x = 7.0,  y = txt_y, label = "PolyI:C",
             size = 2.3, color = "#2ECC71", vjust = 1.0) +
    # 8h PolyI:C (G14)
    annotate("segment", x = 13.6, xend = 14.4, y = bar_y, yend = bar_y,
             color = "#2ECC71", linewidth = 1.1) +
    annotate("text",   x = 14.0, y = txt_y, label = "PolyI:C",
             size = 2.3, color = "#2ECC71", vjust = 1.0) +

    # ── Scales ────────────────────────────────────────────────────────────
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_color_manual(values = GROUP_COLORS, guide = "none") +
    scale_x_discrete(labels = NULL, expand = expansion(add = 0.55)) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.0))) +
    coord_cartesian(ylim = c(y_bot, y_top), clip = "off") +

    labs(
      title    = gene,
      subtitle = if (nzchar(annot_txt)) annot_txt else NULL,
      x        = NULL,
      y        = expression(log[2](Fold~Change))
    ) +

    theme_classic(base_size = 10) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.text.y     = element_text(size = 8),
      axis.title.y    = element_text(size = 8.5, face = "bold"),
      plot.title      = element_text(face = "bold", size = 10.5,
                                     hjust = 0.5, margin = margin(b = 1)),
      plot.subtitle   = element_text(size = 6.3, family = "mono", color = "grey30",
                                     lineheight = 1.20, hjust = 0,
                                     margin = margin(t = 0, b = 22)),
      plot.margin     = margin(t = 28, r = 10, b = 62, l = 10),
      legend.position = "none"
    )
}

# ==============================================================================
# §13  STANDALONE LEGEND PANEL  (placed in 9th patchwork slot)
# ==============================================================================

# Extract the legend grob from a fully rendered ggplot — the only reliable way to
# get key_glyph to render when alpha or other aesthetics would suppress it in a
# standalone panel.
legend_df <- tibble(
  x        = factor(LEGEND_BREAKS, levels = GRP_LEVELS),
  y        = 1,
  Group_ID = factor(LEGEND_BREAKS, levels = GRP_LEVELS)
)

p_legend_src <- ggplot(legend_df, aes(x = x, y = y, fill = Group_ID)) +
  geom_tile(width = 0.7, height = 0.5, key_glyph = draw_key_split_grp) +
  scale_fill_manual(
    values = GROUP_COLORS,
    breaks = LEGEND_BREAKS,
    labels = LEGEND_LABELS,
    name   = "Treatment"
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9.5,
                                     lineheight = 1.3, margin = margin(b = 6)),
    legend.text       = element_text(size = 9.0),
    legend.key.size   = unit(0.72, "cm"),
    legend.key.height = unit(0.72, "cm"),
    legend.key.width  = unit(0.90, "cm"),
    legend.spacing.y  = unit(0.25, "cm"),
    legend.key        = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Pull the legend grob out of the rendered gtable and wrap it for patchwork
{
  gt      <- ggplot_gtable(ggplot_build(p_legend_src))
  leg_idx <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  p_legend <- if (length(leg_idx) > 0)
    patchwork::wrap_elements(gt$grobs[[leg_idx[1]]])
  else
    ggplot() + theme_void()
}

# ==============================================================================
# §14  ASSEMBLE: 8 GENE PANELS + 1 LEGEND PANEL  (3 × 3 grid)
# ==============================================================================
message("\n=== STAGE 4: BUILDING & SAVING FIGURE ===")

n_cols_fig <- 3L
n_rows_fig <- ceiling((n_genes + 1L) / n_cols_fig)   # +1 for legend slot
fig_w      <- n_cols_fig * 5.8 + 0.6
fig_h      <- n_rows_fig * 5.5 + 2.5

plot_list <- lapply(target_genes, build_gene_plot)
names(plot_list) <- target_genes

p_combined <-
  wrap_plots(c(plot_list, list(p_legend)), ncol = n_cols_fig) +
  plot_annotation(
    title    = "HMC3-IAV / αSyn — RT-qPCR Three-Way RCBD Analysis",
    subtitle = paste0(
      "Model: lm(log₂(FC) ~ Timepoint×LPS×Virus + Bio_Rep)  |  ",
      "Strains: ", STRAIN_V1, " / ", STRAIN_V2, "  |  ",
      "n = ", n_bioreps, " biological replicates  |  ",
      "Global Holm FWER per gene (7 factorial + 4 PolyI:C contrasts per timepoint)  |  ",
      "All data inclusive"
    ),
    caption  = paste0(
      "Violin: kernel density (no trim, adjust=2)  ·  Box: IQR + median  ·  ",
      "Dots: individual replicates  ·  ",
      "Dashed = log₂(FC) = 0  ·  Dotted vline = 4h | 8h  ·  ",
      "Brackets = Holm-adjusted p  ·  ",
      "ANOVA stats (Timepoint:LPS:Virus / LPS:Virus / Timepoint:LPS) in subtitle of each panel"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12.5, hjust = 0.5),
      plot.subtitle = element_text(size  = 8,  color = "grey40", hjust = 0.5,
                                   margin = margin(b = 6)),
      plot.caption  = element_text(size  = 7,  color = "grey50", hjust = 0)
    )
  )

ggsave(paste0(FIG_OUT, ".pdf"), p_combined,
       width = fig_w, height = fig_h, device = cairo_pdf)
ggsave(paste0(FIG_OUT, ".png"), p_combined,
       width = fig_w, height = fig_h, dpi = 1000)
message("  Saved: ", basename(paste0(FIG_OUT, ".pdf")),
        "  (", round(fig_w, 1), " × ", round(fig_h, 1), " in)")

# ==============================================================================
# §15  EXCEL EXPORT  (4 sheets)
# ==============================================================================
message("\n=== STAGE 5: EXCEL EXPORT ===")

wb     <- createWorkbook()
hl_hdr <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2")
hl_sig <- createStyle(fgFill = "#E2EFDA")
hl_qc  <- createStyle(fgFill = "#FCE4D6")

write_sheet <- function(wb, name, df, sig_col = NULL,
                        qc_col = NULL, qc_thresh = NULL, qc_dir = ">") {
  addWorksheet(wb, name)
  writeData(wb, name, df)
  addStyle(wb, name, hl_hdr, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
  if (!is.null(sig_col) && sig_col %in% colnames(df)) {
    sig_rows <- which(!is.na(df[[sig_col]]) & df[[sig_col]] < 0.05) + 1L
    if (length(sig_rows) > 0)
      addStyle(wb, name, hl_sig, rows = sig_rows,
               cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
  }
  if (!is.null(qc_col) && !is.null(qc_thresh) && qc_col %in% colnames(df)) {
    qc_rows <- (if (qc_dir == ">") which(df[[qc_col]] > qc_thresh)
                else which(df[[qc_col]] < qc_thresh)) + 1L
    if (length(qc_rows) > 0)
      addStyle(wb, name, hl_qc, rows = qc_rows,
               cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
  }
}

s0 <- data.frame(
  Section = c("Design","Input","Strains","Amplicon QC","Normalization",
              "Target genes","Model","Blocking","Type III SS",
              "Factorial contrasts","PolyI:C contrasts","Multiple testing",
              "Output","Software"),
  Description = c(
    "2(Timepoint:4h/8h) x 2(LPS:No/Yes) x 3(Virus:None/A/PR/8/34/A/NWS/33); G7/G14 = PolyI:C",
    basename(INPUT_PER_REP),
    paste0("V1=", STRAIN_V1, " (H1N1);  V2=", STRAIN_V2, " (H1N1)"),
    "No additional filter — input contains melt-QC-clean wells only (Exclude.amplicon == 0)",
    "geNORM geometric mean of GAPDH + RPS11",
    paste(target_genes, collapse = ", "),
    "lm(log2(FC) ~ Timepoint*LPS*Virus + Bio_Rep)  [RCBD fixed-block]",
    "Bio_Rep fixed factor: absorbs plate-to-plate variation (100% stable, no singularity)",
    "car::Anova(type=3); contr.sum/contr.poly",
    paste0(
      "7 planned contrasts per timepoint: (1) Ctrl vs LPS; ",
      "(2) Ctrl vs A/PR/8/34; (3) Ctrl vs A/NWS/33; (4) A/PR/8/34 vs A/NWS/33 (-LPS); ",
      "(5) LPS vs LPS+A/PR/8/34; (6) LPS vs LPS+A/NWS/33; (7) LPS+strains vs each (+LPS)"
    ),
    "4 Welch t-tests per timepoint (log2): Ctrl/LPS/A/PR/8/34/A/NWS/33 vs PolyI:C",
    "Global Holm FWER per gene: 14 factorial + 8 PolyI:C = 22 p-values pooled per gene",
    "04_master_expression_log2_violin.pdf/png — log2(FC) violin-boxplot; 3x3 grid with legend in slot 9",
    paste0("R ", paste(R.version$major, R.version$minor, sep = "."),
           "; car, emmeans, ggsignif, patchwork, openxlsx, grid")
  ),
  stringsAsFactors = FALSE
)

s1 <- omnibus_df %>%
  select(Gene, Term,
         any_of(c("Df","Sum Sq","Mean Sq","F value","F","Pr(>F)","Pr(>Chisq)","eta2p")),
         Significance) %>%
  arrange(Gene, Term)

s2 <- contrasts_df %>%
  select(Gene, Timepoint, Contrast_Family, Contrast_Label,
         Group_ref, Group_test, estimate, SE, df, t_stat,
         P_raw, P_adjusted_Holm, Significance) %>%
  arrange(Gene, Timepoint, Contrast_Family, Contrast_Label)

s3 <- summary_df %>%
  select(Target_Name, Group_ID, Label, Timepoint, LPS, Virus,
         in_factorial, N_Bio_Reps, N_obs,
         Mean_FC, SEM_FC, Mean_log2FC, SD_log2FC, SEM_log2FC) %>%
  arrange(Target_Name, Group_ID)

addWorksheet(wb, "0_Methods")
writeData(wb, "0_Methods", s0)
addStyle(wb, "0_Methods", hl_hdr, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "0_Methods", cols = 1:2, widths = c(26, 130))

write_sheet(wb, "1_Omnibus_ANOVA",      s1, sig_col = "Pr(>F)")
write_sheet(wb, "2_All_Contrasts",      s2, sig_col = "P_adjusted_Holm")
write_sheet(wb, "3_Group_Summaries_QC", s3)

worksheetOrder(wb) <- c(
  which(names(wb) == "0_Methods"),
  which(names(wb) == "1_Omnibus_ANOVA"),
  which(names(wb) == "2_All_Contrasts"),
  which(names(wb) == "3_Group_Summaries_QC")
)

saveWorkbook(wb, XLSX_OUT, overwrite = TRUE)
message("  Saved: ", basename(XLSX_OUT))

# ==============================================================================
# §16  CONSOLE SUMMARY
# ==============================================================================
message("\n=== COMPLETE ===")
message("  Figure : ", basename(paste0(FIG_OUT, ".pdf")))
message("  XLSX   : ", basename(XLSX_OUT))

p_col_om <- intersect(c("Pr(>F)", "Pr(>Chisq)"), colnames(omnibus_df))[1]
if (!is.na(p_col_om)) {
  sig_om <- omnibus_df %>%
    filter(.data[[p_col_om]] < 0.05,
           !Term %in% c("(Intercept)", "Intercept", "Residuals")) %>%
    mutate(p_fmt = sprintf("%.4f", .data[[p_col_om]])) %>%
    select(Gene, Term, p_fmt, Significance, any_of("eta2p"))
  message("\nSignificant omnibus effects (p < 0.05):")
  if (nrow(sig_om) > 0) print(as_tibble(sig_om)) else message("  (none)")
}

if ("P_adjusted_Holm" %in% colnames(contrasts_df)) {
  sig_c <- contrasts_df %>%
    filter(!is.na(P_adjusted_Holm), P_adjusted_Holm < 0.05) %>%
    mutate(p_fmt = sprintf("%.4f", P_adjusted_Holm)) %>%
    select(Gene, Timepoint, Contrast_Family, Contrast_Label, p_fmt, Significance)
  message("\nSignificant planned contrasts (Holm p < 0.05):")
  if (nrow(sig_c) > 0) print(as_tibble(sig_c)) else message("  (none)")
}

message("\nDone.")
