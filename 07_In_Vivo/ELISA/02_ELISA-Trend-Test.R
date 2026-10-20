#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(readxl) })

HERE  <- "ELISA"
XLS   <- file.path(HERE, "20260901 ELISA Lung and Brain 5day post infection.xlsx")
OUTT  <- file.path(HERE, "tables"); dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)

GRP   <- c("Control", "PR8", "S33")
SCORE <- c(Control = 0, PR8 = 1, S33 = 2)   

E <- as.data.frame(read_excel(XLS, sheet = "R_data"))
names(E) <- c("Group", "Tissue", "Rep", "Cytokine", "Conc")
E$Group[E$Group == "S88"] <- "S33"
E <- E[E$Group %in% GRP & is.finite(E$Conc), ]
E$Group <- factor(E$Group, levels = GRP)
E$y     <- log10(pmax(E$Conc, 1e-9))

panels <- unique(E[, c("Tissue", "Cytokine")])
panels <- panels[order(panels$Tissue, panels$Cytokine), ]

jt_stat <- function(y, g) {           # g: factor da sap thu tu
  lv <- levels(g); s <- 0
  for (a in seq_len(length(lv) - 1L)) for (b in (a + 1L):length(lv)) {
    xa <- y[g == lv[a]]; xb <- y[g == lv[b]]
    s <- s + sum(outer(xb, xa, ">")) + 0.5 * sum(outer(xb, xa, "=="))
  }
  s
}

label_grid <- function(ns) {
  n <- sum(ns); idx <- seq_len(n)
  c1 <- combn(n, ns[1])
  out <- vector("list", 0L)
  for (k in seq_len(ncol(c1))) {
    rest <- setdiff(idx, c1[, k])
    c2 <- combn(rest, ns[2])
    for (m in seq_len(ncol(c2))) {
      lab <- integer(n)
      lab[c1[, k]] <- 1L; lab[c2[, m]] <- 2L; lab[lab == 0L] <- 3L
      out[[length(out) + 1L]] <- lab
    }
  }
  do.call(rbind, out)
}

GRID_CACHE <- new.env()
get_grid <- function(ns) {
  key <- paste(ns, collapse = "-")
  if (is.null(GRID_CACHE[[key]])) GRID_CACHE[[key]] <- label_grid(ns)
  GRID_CACHE[[key]]
}

jt_exact <- function(y, g) {
  g   <- droplevels(g)
  ns  <- as.integer(table(g))
  ord <- order(g)
  y   <- y[ord]; g <- g[ord]
  obs <- jt_stat(y, g)
  G   <- get_grid(ns)
  lv  <- levels(g)
  null <- apply(G, 1L, function(lab) jt_stat(y, factor(lv[lab], levels = lv)))
  p_inc <- mean(null >= obs)                     # one-sided: increasing
  p_dec <- mean(null <= obs)                     # one-sided: decreasing
  list(stat = obs, n_perm = nrow(G),
       p_inc = p_inc, p_dec = p_dec,
       p_two = min(1, 2 * min(p_inc, p_dec)))
}


res <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  ti <- panels$Tissue[i]; cy <- panels$Cytokine[i]
  d  <- E[E$Tissue == ti & E$Cytokine == cy, ]
  if (nlevels(droplevels(d$Group)) < 3L) return(NULL)
  d$x <- SCORE[as.character(d$Group)]

  # partition the 2-df group sum of squares into linear (1) + departure (1)
  m_full <- lm(y ~ Group, d)
  m_lin  <- lm(y ~ x,     d)
  m_null <- lm(y ~ 1,     d)
  a_lin  <- anova(m_null, m_lin)                 # linear vs constant
  a_dev  <- anova(m_lin,  m_full)                # departure from linearity
  p_omni <- anova(m_full)$`Pr(>F)`[1]

  jt <- jt_exact(d$y, d$Group)

  mn <- tapply(d$Conc, d$Group, mean)
  data.frame(Tissue = ti, Cytokine = cy,
             Control = mn[["Control"]], PR8 = mn[["PR8"]], S33 = mn[["S33"]],
             slope_log10 = unname(coef(m_lin)["x"]),
             fold_per_step = 10^unname(coef(m_lin)["x"]),
             p_omni       = p_omni,
             p_trend      = a_lin$`Pr(>F)`[2],
             p_dev        = a_dev$`Pr(>F)`[2],
             JT           = jt$stat,
             p_JT_inc     = jt$p_inc,
             p_JT_two     = jt$p_two,
             monotone     = (mn[["Control"]] < mn[["PR8"]] & mn[["PR8"]] < mn[["S33"]]) |
                            (mn[["Control"]] > mn[["PR8"]] & mn[["PR8"]] > mn[["S33"]]),
             stringsAsFactors = FALSE)
}))

for (col in c("p_omni", "p_trend", "p_JT_two")) {
  q <- rep(NA_real_, nrow(res))
  for (ti in unique(res$Tissue)) { k <- res$Tissue == ti
    q[k] <- p.adjust(res[[col]][k], "BH") }
  res[[sub("^p_", "q_", col)]] <- q
}


limma_tab <- NULL
if (requireNamespace("limma", quietly = TRUE)) {
  limma_tab <- do.call(rbind, lapply(unique(E$Tissue), function(ti) {
    d  <- E[E$Tissue == ti, ]
    d$id <- paste(d$Group, d$Rep, sep = "_")
    cys <- sort(unique(d$Cytokine))
    ids <- sort(unique(d$id))
    M <- matrix(NA_real_, length(cys), length(ids), dimnames = list(cys, ids))
    M[cbind(match(d$Cytokine, cys), match(d$id, ids))] <- d$y
    keep <- rowSums(is.finite(M)) == ncol(M)
    M <- M[keep, , drop = FALSE]
    grp <- sub("_.*", "", colnames(M))
    des <- model.matrix(~ SCORE[grp])           # 1 df: the same linear contrast
    fit <- limma::eBayes(limma::lmFit(M, des))
    tt  <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    data.frame(Tissue = ti, Cytokine = rownames(tt),
               p_limma_trend = tt$P.Value, q_limma_trend = tt$adj.P.Val,
               stringsAsFactors = FALSE)
  }))
  res <- merge(res, limma_tab, by = c("Tissue", "Cytokine"), all.x = TRUE)
}

res <- res[order(res$Tissue, res$p_trend), ]

fmt <- function(p) ifelse(is.na(p), "     -", formatC(p, format = "f", digits = 4))
star <- function(q) ifelse(is.na(q), "", ifelse(q < 0.05, " *", ifelse(q < 0.10, " .", "")))

for (ti in unique(res$Tissue)) {
  z <- res[res$Tissue == ti, ]
  cat(sprintf("\n\n########## %s -- %d panel, FDR within tissue ##########\n", toupper(ti), nrow(z)))
  cat(sprintf("%-11s %9s %9s %9s %7s | %8s %8s %8s %8s %8s | %8s %s\n",
      "cytokine", "Control", "PR8", "S33", "x/step",
      "p_omni", "p_trend", "q_trend", "p_JT", "q_JT", "p_dev", "monotone"))
  cat(strrep("-", 132), "\n")
  for (i in seq_len(nrow(z))) with(z[i, ],
    cat(sprintf("%-11s %9.3g %9.3g %9.3g %7.2f | %8s %8s %8s%-2s %8s %8s%-2s | %8s %s\n",
        Cytokine, Control, PR8, S33, fold_per_step,
        fmt(p_omni), fmt(p_trend), fmt(q_trend), star(q_trend),
        fmt(p_JT_two), fmt(q_JT_two), star(q_JT_two),
        fmt(p_dev), ifelse(monotone, "yes", "NO"))))
}

cat("\n\n########## How to read the table ##########\n")
cat("q_trend / q_JT  : primary result of this family of tests (FDR within tissue).\n")
cat("p_dev           : departure from linearity. SMALL p_dev => non-linear relation;\n")
cat("                  p_trend must NOT be reported for that panel.\n")
cat("monotone = NO   : the three group means are not ordered; the trend is meaningless.\n")
cat("p_JT            : non-parametric, exact. If p_trend and p_JT disagree strongly\n")
cat("                  the result depends on distributional assumptions -> not reported.\n")

if (!is.null(limma_tab)) {
  cat("\n########## (D) limma/eBayes -- secondary check ##########\n")
  cat("Shrinks variances across the cytokines of ONE tissue. Assumes the cytokine\n")
  cat("variances are exchangeable after log10 -- may fail if assay precision differs\n")
  cat("strongly between cytokines. Used for comparison, NOT as the primary result.\n\n")
  z <- res[order(res$Tissue, res$p_limma_trend), ]
  cat(sprintf("%-7s %-11s %10s %10s %10s\n", "tissue", "cytokine", "p_trend", "p_limma", "q_limma"))
  for (i in seq_len(nrow(z))) with(z[i, ],
    cat(sprintf("%-7s %-11s %10s %10s %10s%s\n", Tissue, Cytokine,
        fmt(p_trend), fmt(p_limma_trend), fmt(q_limma_trend), star(q_limma_trend))))
}

out <- file.path(OUTT, "ELISA_trend_test.csv")
write.csv(res, out, row.names = FALSE)
cat(sprintf("\n[saved] %s\n", out))
