#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(readxl) })
source("07_In_Vivo/Flow/Config.R")                 # INVIVO_OUT
ROOT <- "."                                        # ELISA/ and Body weight/ data folders
OUTT <- file.path(INVIVO_OUT, "out_brain", "_tables")
GRP  <- c("Control", "PR8", "S33")
EXCLUDE <- "S33_4"
set.seed(20260903)


B <- read.csv(file.path(OUTT, "mouse_level_table.csv")); B$id <- B$unit
L <- read.csv(file.path(INVIVO_OUT, "out_lung", "_tables", "lung_per_mouse_flowjo.csv"))
L$id <- L$unit
L$pct_AM   <- 100 * L$n_AM   / pmax(L$n_leuk, 1)
L$pct_Mono <- 100 * L$n_Mono / pmax(L$n_leuk, 1)
L$pct_DC   <- 100 * L$n_DC   / pmax(L$n_leuk, 1)

E <- as.data.frame(read_excel(file.path(ROOT, "ELISA/20260901 ELISA Lung and Brain 5day post infection.xlsx"),
                              sheet = "R_data"))
names(E) <- c("Group", "Tissue", "Rep", "Cytokine", "Conc")
E$Group[E$Group == "S88"] <- "S33"
E$id  <- paste(E$Group, E$Rep, sep = "_")
E$var <- paste0("elisa.", E$Tissue, ".", gsub("[^A-Za-z0-9]", "", E$Cytokine))
EW <- reshape(E[, c("id", "var", "Conc")], idvar = "id", timevar = "var", direction = "wide")
names(EW) <- sub("^Conc\\.", "", names(EW))

W <- read.csv(file.path(ROOT, "Body weight/tables/BW02_per_mouse_summary.csv"),
              fileEncoding = "UTF-8-BOM")
map <- c("Naïve"="Control","Naive"="Control","A/PR8/43"="PR8","A/PR/8/34"="PR8",
         "A/NWS/33"="S33","S33"="S33")
W$group <- unname(map[W$Group]); W$id <- paste(W$group, sub(".*_", "", W$ID), sep = "_")

D <- merge(B[, c("id","group","Y_uG_aSyn_MFI","Y_nUG_aSyn_pct","uG_CD206_MFI","pct_microglia")],
           L[, c("id","npct_epi","npct_leuk","npmfi_epi","pct_AM","pct_Mono","pct_DC")], by = "id")
D <- merge(D, EW, by = "id")
D <- merge(D, W[, c("id","AUC_loss","nadir_pct")], by = "id")
D <- D[!(D$id %in% EXCLUDE), ]
D$group <- factor(D$group, levels = GRP)

X_LUNG   <- "npct_epi"                       # primary lung variable, fixed in advance
Y_BRAIN  <- c("Y_uG_aSyn_MFI", "Y_nUG_aSyn_pct", "uG_CD206_MFI", "pct_microglia")
NICE <- c(Y_uG_aSyn_MFI  = "a-Syn MFI, microglia",
          Y_nUG_aSyn_pct = "% a-Syn+, non-microglia",
          uG_CD206_MFI   = "CD206 MFI, microglia",
          pct_microglia  = "% microglia")

cat(sprintf("n = %d animals | lung variable = %s | %d brain read-outs (fixed in advance)\n",
            nrow(D), X_LUNG, length(Y_BRAIN)))
cat(sprintf("excluded: %s\n", paste(EXCLUDE, collapse = ", ")))

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
NPERM <- 20000L

rset <- function(x, y_mat, strat = NULL) {
  if (!is.null(strat)) {
    x <- residuals(lm(x ~ strat))
    y_mat <- apply(y_mat, 2L, function(y) residuals(lm(y ~ strat)))
  }
  apply(y_mat, 2L, function(y) suppressWarnings(cor(x, y)))
}

perm_max <- function(x, y_mat, strat = NULL, groups = NULL) {
  obs <- max(abs(rset(x, y_mat, strat)), na.rm = TRUE)
  null <- replicate(NPERM, {
    xp <- if (is.null(groups)) sample(x) else {
      z <- x; for (g in levels(droplevels(groups))) { k <- which(groups == g)
        if (length(k) > 1) z[k] <- z[sample(k)] }; z }
    max(abs(rset(xp, y_mat, strat)), na.rm = TRUE) })
  list(obs = obs, p = (sum(null >= obs) + 1) / (NPERM + 1),
       q95 = unname(quantile(null, .95)), med = median(null))
}

report <- function(tag, d, strat_by_group, perm_within_group, note) {
  x  <- d[[X_LUNG]]
  Y  <- as.matrix(d[, Y_BRAIN])
  st <- if (strat_by_group) droplevels(d$group) else NULL
  gg <- if (perm_within_group) droplevels(d$group) else NULL

  cat(sprintf("\n\n########## %s  (n = %d) ##########\n", tag, nrow(d)))
  cat(note, "\n\n")
  rr <- rset(x, Y, st)
  df <- if (is.null(st)) nrow(d) - 2 else nrow(d) - nlevels(st) - 1
  pv <- sapply(rr, function(r) { t <- r * sqrt(df / (1 - r^2)); 2 * pt(-abs(t), df) })
  q  <- p.adjust(pv, "BH")
  cat(sprintf("%-26s %7s %5s %9s %9s\n", "brain read-out", "r", "df", "p", "q(BH,4)"))
  for (k in seq_along(rr))
    cat(sprintf("%-26s %+7.3f %5d %9.4f %9.4f%s\n", NICE[[Y_BRAIN[k]]], rr[k], df, pv[k], q[k],
                ifelse(q[k] < .05, " *", ifelse(q[k] < .10, " .", ""))))
  pm <- perm_max(x, Y, st, gg)
  cat(sprintf("\n  permutation on the largest |r| of the 4: observed %.3f | null median %.3f Q95 %.3f | p = %.4f%s\n",
              pm$obs, pm$med, pm$q95, pm$p, ifelse(pm$p < .05, "  <= exceeds chance", "")))
  data.frame(level = tag, readout = Y_BRAIN, r = unname(rr), df = df, p = unname(pv),
             q = unname(q), perm_p_family = pm$p, n = nrow(d), stringsAsFactors = FALSE)
}

R1 <- report("(i) 11 animals -- raw correlation", D, FALSE, FALSE,
  "CONFOUNDED BY GROUP. Both variables track the group, so this correlation appears\neven without any individual-level relationship. NOT used as evidence; for completeness only.")

R2 <- report("(ii) WITHIN-GROUP residuals -- group effect removed", D, TRUE, TRUE,
  "Within one group, does an animal with a more severe lung also have a different brain?\nThis is the INDIVIDUAL-level relationship. Permutation preserves the group structure.")

Dinf <- D[D$group != "Control", ]
R3 <- report("(iii) ONLY the 7 infected animals (PR8 + S33)", Dinf, FALSE, FALSE,
  "Naive anchor removed. Among animals that ARE infected, does viral dose set the size of\nthe brain change -- a genuine dose-response -- or is it only 'infected vs not'?")

AX <- rbind(R1, R2, R3)
write.csv(AX, file.path(OUTT, "AXIS_lung_to_brain.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# B. Robustness: swap the lung variable. Does the conclusion depend on npct_epi?
# ---------------------------------------------------------------------------
cat("\n\n########## B. Swapping the lung variable -- is the conclusion robust? ##########\n")
cat("npct_epi is the choice FIXED IN ADVANCE. This table shows whether other lung\n")
cat("severity measures lead to the same conclusion. The best-looking row of this\n")
cat("table must NOT be chosen as the primary result.\n\n")
alt <- c("npct_epi","npct_leuk","npmfi_epi","pct_AM","pct_Mono","pct_DC","AUC_loss")
cat(sprintf("%-12s %s\n", "lung var", paste(sprintf("%22s", NICE[Y_BRAIN]), collapse = "")))
SENS <- do.call(rbind, lapply(alt, function(v) {
  x <- residuals(lm(D[[v]] ~ D$group))
  rr <- sapply(Y_BRAIN, function(y) cor(x, residuals(lm(D[[y]] ~ D$group))))
  df <- nrow(D) - 3 - 1
  pv <- sapply(rr, function(r) 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df))
  cat(sprintf("%-12s %s\n", v,
      paste(sprintf("%16s%6s", sprintf("r=%+.2f", rr), sprintf("p=%.3f", pv)), collapse = "")))
  data.frame(lung_variable = v, readout = Y_BRAIN, r = unname(rr), p = unname(pv)) }))
write.csv(SENS, file.path(OUTT, "AXIS_lung_variable_sensitivity.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# C. GROUP-level comparison: do the lung severity order and the brain order agree?
# ---------------------------------------------------------------------------
cat("\n\n########## C. GROUP level: do the two scales rank the groups in the same order? ##########\n")
cat("Independent of ID matching and of n = 11; the most robust statement these\n")
cat("data support.\n\n")
gm <- function(v) tapply(D[[v]], D$group, mean)
cat(sprintf("%-26s %10s %10s %10s   %s\n", "", "Control", "PR8", "S33", "order"))
for (v in c(X_LUNG, "pct_AM", Y_BRAIN)) { m <- gm(v)
  ord <- if (m[1] < m[2] && m[2] < m[3]) "increasing" else
         if (m[1] > m[2] && m[2] > m[3]) "decreasing" else "NOT monotonic"
  cat(sprintf("%-26s %10.4g %10.4g %10.4g   %s\n",
      ifelse(v %in% names(NICE), NICE[[v]], v), m[1], m[2], m[3], ord)) }

# ---------------------------------------------------------------------------
# D. Does the path run through a BRAIN cytokine?  (limited by n = 11)
# ---------------------------------------------------------------------------
cat("\n\n########## D. Does a brain cytokine lie on the path? ##########\n")
cat("Formal mediation needs a much larger n. What is possible: test the three legs\n")
cat("of the path and report the limitation. IL-1beta is used because it is the ONLY\n")
cat("brain cytokine with a signal (trend test, ELISA/tables/ELISA_trend_test.csv).\n\n")
MED <- "elisa.Brain.IL1"
if (MED %in% names(D)) {
  a <- cor.test(D[[X_LUNG]], D[[MED]])
  cat(sprintf("  a) lung -> brain IL-1b      r = %+.3f  p = %.4f\n", a$estimate, a$p.value))
  for (y in Y_BRAIN) {
    b  <- cor.test(D[[MED]], D[[y]])
    cf <- coef(summary(lm(D[[y]] ~ D[[X_LUNG]] + D[[MED]])))
    c0 <- coef(summary(lm(D[[y]] ~ D[[X_LUNG]])))
    cat(sprintf("  b) IL-1b -> %-24s r = %+.3f p = %.4f | c) lung coefficient: raw %+.4g -> adjusted %+.4g (%+.0f%%)\n",
        NICE[[y]], b$estimate, b$p.value, c0[2,1], cf[2,1], 100*(cf[2,1]-c0[2,1])/abs(c0[2,1]))) }
  cat("\n  Reading: legs (a) and (b) must both be real for the path to exist. If the lung\n")
  cat("  coefficient barely changes when IL-1b is added, IL-1b does NOT mediate. With\n")
  cat("  n = 11 weak mediation cannot be excluded -- only strong mediation can.\n")
}

cat("\n[saved] AXIS_lung_to_brain.csv, AXIS_lung_variable_sensitivity.csv\n")


cat("\n\n########## E. PR8 vs S33 in the periphery -- why does only S33 change the brain? ##########\n")
L2 <- L; L2$pct_Eos <- 100*L2$n_Eos/pmax(L2$n_leuk,1); L2$pct_Neu <- 100*L2$n_Neu/pmax(L2$n_leuk,1)
L2$pct_leuk <- 100*L2$n_leuk/pmax(L2$n_live,1); L2$pct_live <- 100*L2$n_live/pmax(L2$n_singlet,1)
P <- merge(L2[, c("id","group","npct_epi","npct_leuk","npmfi_epi","npmfi_leuk",
                  "pct_AM","pct_Mono","pct_DC","pct_Eos","pct_Neu","pct_leuk","pct_live")],
           EW, by = "id")
P <- merge(P, W[, c("id","AUC_loss","nadir_pct")], by = "id")
P <- P[P$group %in% c("PR8","S33"), ]

pr8_vs_s33 <- function(d, tag) {
  vs <- setdiff(names(d), c("id","group","unit"))
  out <- do.call(rbind, lapply(vs, function(v) {
    a <- d[[v]][d$group=="PR8"]; b <- d[[v]][d$group=="S33"]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 3 || length(b) < 3) return(NULL)
    # exact permutation on the mean difference (choose(n,k) combinations, enumerated)
    z <- c(a,b); n1 <- length(a); cb <- combn(length(z), n1)
    null <- apply(cb, 2, function(k) mean(z[k]) - mean(z[-k]))
    obs <- mean(a) - mean(b)
    data.frame(variable = v, PR8 = mean(a), S33 = mean(b),
               ratio = mean(b)/mean(a),
               p_perm = mean(abs(null) >= abs(obs) - 1e-12),
               n_perm = ncol(cb), stringsAsFactors = FALSE) }))
  out$FDR <- p.adjust(out$p_perm, "BH"); out <- out[order(out$p_perm), ]
  cat(sprintf("\n--- %s (PR8 n=%d vs S33 n=%d) ---\n", tag,
              sum(d$group=="PR8"), sum(d$group=="S33")))
  cat(sprintf("%-24s %11s %11s %8s %9s %8s\n","variable","PR8","S33","S33/PR8","p_perm","FDR"))
  for (i in seq_len(min(14, nrow(out)))) with(out[i,],
    cat(sprintf("%-24s %11.4g %11.4g %8.2f %9.4f %8.4f%s\n", variable, PR8, S33, ratio, p_perm, FDR,
        ifelse(FDR<.05," *", ifelse(FDR<.10," .","")))))
  out
}
E12 <- pr8_vs_s33(P,                        "S33_4 KEPT     (4 vs 4) -- preferred for lung/ELISA/body weight")
E11 <- pr8_vs_s33(P[!(P$id %in% EXCLUDE),], "S33_4 EXCLUDED (4 vs 3) -- matches the 11-animal brain set")
cat("\nsmallest attainable p_perm: 4vs4 -> 2/70 = 0.0286 | 4vs3 -> 2/35 = 0.0571\n")
cat("With 4vs3 NO variable can reach FDR < 0.05 after adjustment, which is why the\n")
cat("4vs4 version is preferred for the three non-brain blocks.\n")
write.csv(cbind(E12, version="4vs4"), file.path(OUTT,"AXIS_PR8_vs_S33_periphery.csv"), row.names=FALSE)
cat("\n[saved] AXIS_PR8_vs_S33_periphery.csv\n")
