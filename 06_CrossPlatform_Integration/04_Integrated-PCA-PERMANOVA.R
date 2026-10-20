# ==============================================================================
# 04_Integrated-PCA-PERMANOVA.R  —  Integrated multivariate analysis
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(vegan); library(ggrepel); library(patchwork) })
PROJ <- "."
DIR  <- file.path(PROJ,"CrossPlatform_aSyn"); OUT <- file.path(DIR,"output"); FIG <- file.path(DIR,"figures")
CODE <- "06_CrossPlatform_Integration" 
source(file.path(CODE,"_Labels.R"))
ANCHOR <- "WB_aSyn_25kDa"; set.seed(20260615)
m <- read_csv(file.path(OUT,"master_matched.csv"), show_col_types=FALSE) %>%
  mutate(timepoint=factor(timepoint,c("4h","8h","24h")), treatment=factor(treatment,c("Control","V1","V2","PolyIC")))
pretty_mk <- disp_marker

run_block <- function(df, feats, tag) {
  d <- df %>% select(timepoint, treatment, all_of(feats)) %>% drop_na()
  X <- scale(as.matrix(d[feats]))
  pca <- prcomp(X, center=FALSE, scale.=FALSE)        # already scaled (correlation PCA)
  ve  <- round(100*summary(pca)$importance[2,1:2],1)
  scores <- bind_cols(d %>% select(timepoint,treatment), as_tibble(pca$x[,1:2])) %>%
    mutate(treatment = factor(disp_treat(as.character(treatment)), levels=TRT_LEVELS))
  load   <- as_tibble(pca$rotation[,1:2], rownames="marker") %>%
    mutate(label=pretty_mk(marker), is_asyn = marker==ANCHOR)
  write_csv(scores, file.path(OUT, paste0("tierC_pca_scores_",tag,".csv")))
  write_csv(load,   file.path(OUT, paste0("tierC_pca_loadings_",tag,".csv")))
  write_csv(tibble(PC=c("PC1","PC2"), pct=ve), file.path(OUT, paste0("tierC_var_",tag,".csv")))


  D <- dist(X, method="euclidean")
  pm <- adonis2(D ~ timepoint*treatment, data=d, permutations=9999, by="terms")
  bd <- bind_rows(lapply(c("timepoint","treatment"), function(f){
    h <- betadisper(D, d[[f]]); pt <- permutest(h, permutations=9999)
    tibble(factor=f, betadisper_F=pt$tab$F[1], betadisper_p=pt$tab$`Pr(>F)`[1]) }))
  pmdf <- as.data.frame(pm) %>% rownames_to_column("term") %>%
    mutate(block=tag) %>% left_join(bd, by=c("term"="factor"))
  write_csv(pmdf, file.path(OUT, paste0("tierC_permanova_",tag,".csv")))

  cat(sprintf("\n[%s] PC1=%.1f%% PC2=%.1f%% | PERMANOVA treatment R²=%.2f p=%.4f, timepoint R²=%.2f p=%.4f\n",
              tag, ve[1], ve[2], pm$R2[2], pm$`Pr(>F)`[2], pm$R2[1], pm$`Pr(>F)`[1]))
  cat("αSyn-25kDa nearest-aligned markers (|cos| with αSyn loading):\n")
  av <- as.numeric(load[load$marker==ANCHOR, c("PC1","PC2")])
  load %>% filter(!is_asyn) %>%
    mutate(cos = (PC1*av[1]+PC2*av[2])/(sqrt(PC1^2+PC2^2)*sqrt(sum(av^2)))) %>%
    arrange(desc(cos)) %>% transmute(label, cos=round(cos,2)) %>% slice(c(1:3, (n()-2):n())) %>%
    as.data.frame() %>% print(row.names=FALSE)
  invisible(pmdf)
}

WB    <- grep("^WB_",    names(m), value=TRUE)   
ELISA <- grep("^ELISA_", names(m), value=TRUE)
FACS  <- grep("^FACS_",  names(m), value=TRUE)
QPCR  <- grep("^qPCR_",  names(m), value=TRUE)

cat("==================== PCA-A (all timepoints, no qPCR) ====================\n")
run_block(m, c(WB, ELISA, FACS), "allTP")
cat("\n==================== PCA-B (4h/8h, +qPCR) ====================\n")
run_block(filter(m, timepoint %in% c("4h","8h")) %>% mutate(timepoint=droplevels(timepoint)),
          c(WB, ELISA, FACS, QPCR), "early_qPCR")
message("\n[04_Integrated-PCA-PERMANOVA] PCA + PERMANOVA done")
