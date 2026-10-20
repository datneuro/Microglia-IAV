TISSUE <- "lung"
source("07_In_Vivo/Flow/Config.R"); source("07_In_Vivo/Flow/Common.R")
suppressMessages({library(readxl); library(ggplot2); library(cowplot)})
ROOT <- "."
GRP <- c("Control","PR8","S33")
COL <- c(Control="#9E9E9E", PR8="#4C72B0", S33="#C44E52")
XLB <- c(Control="Naive", PR8="A/PR/8/34", S33="A/NWS/33")
key <- function(g, i) paste(g, i, sep="_")


B <- read.csv(file.path(INVIVO_OUT,"out_brain","_tables","mouse_level_table.csv"))
B$id <- B$unit
brain_vars <- c("Y_uG_aSyn_pct","Y_uG_aSyn_MFI","Y_nUG_aSyn_pct","uG_CD206_pct",
                "uG_CD206_MFI","pct_microglia","uG_pct_CD45hi")
B <- B[, c("id","group", intersect(brain_vars, names(B)))]
names(B)[-(1:2)] <- paste0("brain.", names(B)[-(1:2)])


Lf <- file.path(INVIVO_OUT,"out_lung","_tables","lung_per_mouse_flowjo.csv")
stopifnot("lung_per_mouse_flowjo.csv missing -- run 79 + 81 first" = file.exists(Lf))
L <- read.csv(Lf)
for (p in c("AM","Mono","DC","Eos","Neu","Mac")) L[[paste0("pct_",p)]] <- 100*L[[paste0("n_",p)]]/pmax(L$n_leuk,1)
L$pct_leuk <- 100*L$n_leuk/pmax(L$n_live,1)
L$pct_live <- 100*L$n_live/pmax(L$n_singlet,1)
keepL <- c("unit","group","pct_live","pct_leuk","pct_AM","pct_Mono","pct_DC",
           "pct_Eos","pct_Neu","npct_epi","npct_leuk","npmfi_epi","npmfi_leuk")
L <- L[, intersect(keepL, names(L))]; names(L)[1] <- "id"
names(L)[-(1:2)] <- paste0("lung.", names(L)[-(1:2)])

Ef <- file.path(ROOT,"ELISA/20260901 ELISA Lung and Brain 5day post infection.xlsx")
E <- as.data.frame(read_excel(Ef, sheet="R_data"))
names(E) <- c("Group","Tissue","Rep","Cytokine","Conc")
E$Group[E$Group %in% c("S88")] <- "S33"
E$id <- key(E$Group, E$Rep)
E$var <- paste0("elisa.", E$Tissue, ".", gsub("[^A-Za-z0-9]","", E$Cytokine))
EW <- reshape(E[,c("id","Group","var","Conc")], idvar=c("id","Group"),
              timevar="var", direction="wide")
names(EW) <- sub("^Conc\\.", "", names(EW)); names(EW)[2] <- "group"

W <- read.csv(file.path(ROOT,"Body weight/tables/BW02_per_mouse_summary.csv"),
              fileEncoding="UTF-8-BOM")
# accepts both the corrected strain label and the misspelt "A/PR8/43" of earlier body-weight tables
map <- c("Naïve"="Control","Naive"="Control","A/PR8/43"="PR8","A/PR/8/34"="PR8",
         "A/NWS/33"="S33","S33"="S33")
W$group <- unname(map[W$Group])
W$idx <- as.integer(sub(".*_", "", W$ID)); W$id <- key(W$group, W$idx)
W <- W[, c("id","group","nadir_pct","pct_lastday","AUC_loss")]
names(W)[-(1:2)] <- paste0("bw.", names(W)[-(1:2)])


EXCLUDE <- "S33_4"
L12 <- L; EW12 <- EW; W12 <- W          # keep the 12-animal versions for comparison
drop_ex <- function(d) if (is.null(d)) NULL else d[!(d$id %in% EXCLUDE), , drop=FALSE]
L <- drop_ex(L); EW <- drop_ex(EW); W <- drop_ex(W); B <- drop_ex(B)
cat(sprintf("Excluded %s from all four blocks (reason, Config.R: %s)\n\n",
            paste(EXCLUDE, collapse=", "), paste(exclude_reason[EXCLUDE], collapse="; ")))

cat("########## The four data blocks ##########\n")
for (nm in c("brain","lung","elisa","bw")) {
  d <- list(brain=B, lung=L, elisa=EW, bw=W)[[nm]]
  if (is.null(d)) { cat(sprintf("  %-6s NOT AVAILABLE\n", nm)); next }
  cat(sprintf("  %-6s %2d animals | %2d variables | %s\n", nm, nrow(d), ncol(d)-2,
      paste(sort(d$id), collapse=" "))) }


cat("\n\n########## TIER 1: do the four blocks tell the same story? ##########\n")
cat("(does not use the ID pairing -- valid whether or not the IDs line up)\n\n")
tier1 <- function(d, tag) { if (is.null(d)) return(NULL)
  vs <- setdiff(names(d), c("id","group")); d$group <- factor(d$group, levels=GRP)
  do.call(rbind, lapply(vs, function(v) { z <- d[is.finite(d[[v]]),]
    if (nlevels(droplevels(z$group)) < 3) return(NULL)
    mn <- tapply(z[[v]], z$group, mean)
    data.frame(block=tag, variable=v, Control=mn[1], PR8=mn[2], S33=mn[3],
      p=anova(lm(as.formula(paste(v,"~group")), z))$`Pr(>F)`[1],
      direction=ifelse(mn[3] > mn[1], "S33 higher", "S33 lower")) })) }
T1 <- do.call(rbind, list(tier1(B,"brain"), tier1(L,"lung"), tier1(EW,"elisa"), tier1(W,"body weight")))
T1$FDR <- p.adjust(T1$p, "BH")
T1 <- T1[order(T1$FDR), ]
cat(sprintf("%-9s %-28s %11s %11s %11s %8s %8s  %s\n","block","variable","Control","PR8","S33","p","FDR","direction"))
for (i in seq_len(nrow(T1))) with(T1[i,], cat(sprintf("%-9s %-28s %11.4g %11.4g %11.4g %8.4f %8.4f  %s%s\n",
    block, variable, Control, PR8, S33, p, FDR, direction, ifelse(FDR<0.05," *",""))))
write.csv(T1, file.path(INVIVO_OUT,"out_brain","_tables","INTEGRATED_tier1_group.csv"), row.names=FALSE)

# --- cross-check: same tier 1 but KEEPING S33_4 in the three non-brain blocks ---
T1b <- do.call(rbind, list(tier1(L12,"lung"), tier1(EW12,"elisa"), tier1(W12,"body weight")))
T1b$FDR <- p.adjust(T1b$p, "BH")
cmp <- merge(T1[,c("variable","p")], T1b[,c("variable","p")], by="variable", suffixes=c("_11mice","_12mice"))
cmp$conclusion_changes <- (cmp$p_11mice < .05) != (cmp$p_12mice < .05)
cat(sprintf("\n-- Keeping vs excluding S33_4: %d/%d variables change conclusion at p = 0.05 --\n",
            sum(cmp$conclusion_changes), nrow(cmp)))
if (any(cmp$conclusion_changes)) { z <- cmp[cmp$conclusion_changes,]
  for (i in seq_len(nrow(z))) with(z[i,], cat(sprintf("   %-28s 11 mice p=%.4f | 12 mice p=%.4f\n", variable, p_11mice, p_12mice))) }
write.csv(T1b, file.path(INVIVO_OUT,"out_brain","_tables","INTEGRATED_tier1_group_12mice.csv"), row.names=FALSE)

cat("\n\n########## TIER 2: animal-level correlation ##########\n")
cat("The ID pairing IS VERIFIED (numbers follow the dissection order, shared by all four blocks).\n")
cat("The permutation test below therefore asks 'is there real between-block structure',\n")
cat("not 'do the IDs line up'.\n\n")
blocks <- Filter(Negate(is.null), list(brain=B, lung=L, elisa=EW, bw=W))
J <- Reduce(function(a,b) merge(a, b[, setdiff(names(b),"group")], by="id", all=FALSE), blocks)
J$group <- factor(blocks[[1]]$group[match(J$id, blocks[[1]]$id)], levels=GRP)
cat(sprintf("%d animals have ALL FOUR blocks: %s\n", nrow(J), paste(sort(J$id), collapse=" ")))
vs <- setdiff(names(J), c("id","group"))
num <- vs[sapply(J[vs], is.numeric)]

# RESIDUAL correlation after removing the group effect. Without this, every pair
# would correlate simply because both separate the three groups, not because of
# any individual-level relationship.
resid_of <- function(v) { z <- J[[v]]; ok <- is.finite(z)
  r <- rep(NA_real_, length(z)); if (sum(ok) < 5) return(r)
  r[ok] <- residuals(lm(z[ok] ~ J$group[ok])); r }
RS <- sapply(num, resid_of); rownames(RS) <- J$id
cross <- function(a, b) grepl(paste0("^",a,"\\."), num) & FALSE
pairs_df <- do.call(rbind, lapply(seq_along(num), function(i) {
  if (i == length(num)) return(NULL)
  do.call(rbind, lapply((i+1):length(num), function(j) {
    ba <- sub("\\..*","",num[i]); bb <- sub("\\..*","",num[j])
    if (ba == bb) return(NULL)                       # BETWEEN blocks only
    x <- RS[,i]; y <- RS[,j]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 6) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="pearson"))
    data.frame(v1=num[i], v2=num[j], n=sum(ok), r=unname(ct$estimate), p=ct$p.value) })) }))
pairs_df$FDR <- p.adjust(pairs_df$p, "BH")
pairs_df <- pairs_df[order(pairs_df$p), ]
cat(sprintf("\n%d BETWEEN-block pairs (group effect removed). Top 15 by p:\n", nrow(pairs_df)))
cat(sprintf("%-30s %-30s %3s %7s %8s %8s\n","variable 1","variable 2","n","r","p","FDR"))
for (i in head(seq_len(nrow(pairs_df)), 15)) with(pairs_df[i,],
  cat(sprintf("%-30s %-30s %3d %+7.3f %8.4f %8.4f%s\n", v1, v2, n, r, p, FDR, ifelse(FDR<0.10," .",""))))
cat(sprintf("\nPairs with FDR < 0.10: %d / %d\n", sum(pairs_df$FDR<0.10), nrow(pairs_df)))

cat("\n---- Counterfactual test: animal indices permuted within each group ----\n")
set.seed(20260903); NPERM <- 2000L
obs_stat <- sum(pairs_df$FDR < 0.10)
obs_max  <- max(abs(pairs_df$r), na.rm=TRUE)
nullv <- replicate(NPERM, {
  RSp <- RS
  for (bk in unique(sub("\\..*","",num))) { ci <- which(sub("\\..*","",num)==bk)
    for (g in GRP) { ri <- which(J$group==g); if (length(ri)>1) RSp[ri,ci] <- RSp[sample(ri),ci] } }
  rs <- sapply(seq_len(nrow(pairs_df)), function(k) {
    x <- RSp[, pairs_df$v1[k]]; y <- RSp[, pairs_df$v2[k]]
    ok <- is.finite(x)&is.finite(y); if (sum(ok)<6) NA else abs(cor(x[ok], y[ok])) })
  max(rs, na.rm=TRUE) })
cat(sprintf("  largest observed |r| = %.3f\n", obs_max))
cat(sprintf("  null distribution (2000 permutations): median %.3f | Q95 %.3f | max %.3f\n",
    median(nullv), quantile(nullv,.95), max(nullv)))
cat(sprintf("  p = %.4f  -> %s\n", (sum(nullv >= obs_max)+1)/(NPERM+1),
    ifelse((sum(nullv>=obs_max)+1)/(NPERM+1) < 0.05,
      "between-block structure exceeds chance",
      "does NOT exceed chance: either no real relationship or the IDs do not line up -- the data CANNOT tell the two apart")))
write.csv(pairs_df, file.path(INVIVO_OUT,"out_brain","_tables","INTEGRATED_tier2_pairs.csv"), row.names=FALSE)
write.csv(J, file.path(INVIVO_OUT,"out_brain","_tables","INTEGRATED_mouse_table.csv"), row.names=FALSE)
cat("\n[saved] INTEGRATED_tier1_group.csv, INTEGRATED_tier2_pairs.csv, INTEGRATED_mouse_table.csv\n")
