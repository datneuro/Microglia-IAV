
source("07_In_Vivo/Flow/Config.R"); source("07_In_Vivo/Flow/Common.R")
suppressMessages({library(lme4); library(lmerTest); library(emmeans)
                  library(ggplot2); library(cowplot); library(ggsignif)})
COL <- c(Control="#9E9E9E", PR8="#4C72B0", S33="#C44E52")   # nguon: ELISA/elisa_barplot_clean.R (group_cols)
XLB <- c(Control="Naive", PR8="A/PR/8/34", S33="A/NWS/33")
ev <- readRDS(outf("events_microglia.rds")); mc <- marker_cutoff()
fmo <- attr(ev,"fmo")
am <- fmo$cutoff[fmo$marker=="aSyn" & fmo$parent=="microglia"]
an <- fmo$cutoff[fmo$marker=="aSyn" & fmo$parent=="non_microglia"]
ev <- ev[!gsub("_","-",ev$unit) %in% gsub("_","-",exclude_units),]
mg <- ev[ev$pop=="microglia",]; nm <- ev[ev$pop=="non_microglia",]
c45 <- asinh(fixed_cutoffs_raw[["CD45"]]/arcsinh_cofactor)
c11 <- asinh(fixed_cutoffs_raw[["CD11b"]]/arcsinh_cofactor)
key <- unique(ev[,c("file","group","unit")]); key$group <- factor(key$group, levels=group_levels)

getv <- function(d, cut, what) {
  r <- key
  r$v <- switch(what,
    mfi = tapply(sinh(d$aSyn)*arcsinh_cofactor, d$file, mean)[r$file],
    pct = 100*tapply(d$aSyn > cut, d$file, mean)[r$file],
    cd206mfi = tapply(sinh(d$CD206)*arcsinh_cofactor, d$file, mean)[r$file],
    cd206pct = 100*tapply(d$CD206 > mc, d$file, mean)[r$file])
  r[is.finite(r$v),]
}
sig_lab <- function(p) if (is.na(p)) "" else
  if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else
  if (p < .10) sprintf("p = %.3f", p) else ""

bar_panel <- function(r, ylab, title) {
  m  <- suppressMessages(lmer(v ~ group + (1|unit), r))
  em <- suppressMessages(emmeans(m, ~group))
  ct <- as.data.frame(pairs(em, adjust="none"))
  pv <- setNames(ct$p.value, ct$contrast)
  mu <- aggregate(v ~ group + unit, r, mean)          # mean per mouse
  su <- do.call(rbind, lapply(split(mu, mu$group), function(z)
        data.frame(group=z$group[1], M=mean(z$v), SE=sd(z$v)/sqrt(nrow(z)))))
  su$group <- factor(su$group, levels=group_levels)
  top <- max(c(su$M+su$SE, mu$v, r$v), na.rm=TRUE)
  ymax <- top*(1.10 + .115*max(1, sum(sapply(c('Control - PR8','Control - S33','PR8 - S33'), function(k) sig_lab(pv[[k]])!=''))))
  cmp <- list(c("Control","PR8"), c("Control","S33"), c("PR8","S33"))
  nmc <- c("Control - PR8","Control - S33","PR8 - S33")
  ann <- sapply(nmc, function(k) sig_lab(pv[[k]]))
  keepi <- which(ann != "")
  g <- ggplot(su, aes(x=group, y=M)) +
    geom_col(data=su, aes(group, M, fill=group), alpha=.6, width=.7, colour=NA) +
    geom_errorbar(data=su, aes(group, ymin=pmax(0,M-SE), ymax=M+SE, colour=group),
                  width=.22, linewidth=.6) +
    geom_jitter(data=r, aes(group, v, colour=group), shape=16, size=.95,
                alpha=.35, width=.13, height=0) +
    geom_jitter(data=mu, aes(group, v, colour=group), shape=21, fill="white",
                size=2.3, stroke=.9, width=.07, height=0) +
    scale_fill_manual(values=COL, guide="none") +
    scale_colour_manual(values=COL, guide="none") +
    scale_x_discrete(labels=XLB) +
    scale_y_continuous(expand=c(0,0), limits=c(0,ymax), oob=scales::oob_keep) +
    labs(x=NULL, y=ylab) +
    theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=9, face="bold"),
          axis.title.y=element_text(size=9.8, margin=margin(r=4)),
          plot.background=element_rect(fill="white", colour=NA),
          plot.margin=margin(8,5,3,5))
  if (length(keepi)) {
    ys <- top*(1.05 + .105*(seq_along(keepi)-1))
    g <- g + geom_signif(data=su, aes(x=group, y=M), inherit.aes=FALSE,
                         xmin=sapply(cmp[keepi], function(z) match(z[1], group_levels)),
                         xmax=sapply(cmp[keepi], function(z) match(z[2], group_levels)),
                         annotations=ann[keepi], y_position=ys, manual=FALSE,
                         tip_length=.012, textsize=3.0, vjust=-.15, colour="black")
  }
  list(plot=g, stats=data.frame(Control=summary(em)$emmean[1],
       PR8=summary(em)$emmean[2], S33=summary(em)$emmean[3],
       p_omni=anova(m)$`Pr(>F)`[1], p_C_PR8=pv[["Control - PR8"]],
       p_C_S33=pv[["Control - S33"]], p_PR8_S33=pv[["PR8 - S33"]]))
}
outd <- file.path(output_dir,"bar_figures"); dir.create(outd, showWarnings=FALSE)
pand <- file.path(outd,"panels"); dir.create(pand, showWarnings=FALSE)
ST <- list()
mk <- function(k, r, yl, ti) {
  o <- bar_panel(r, yl, ti)
  st <- o$stats; st$panel <- k
  ST[[k]] <<- st[, c("panel", setdiff(names(st), "panel"))]
  save_png(file.path(pand, paste0(k, ".png")), o$plot, width=2.9, height=3.3, dpi=1200, bg="white")
  o$plot
}
P1 <- list(
  mk("uG",     getv(mg, am, "mfi"), expression(alpha*"-Syn MFI in microglia"), NULL),
  mk("cd45hi", getv(mg[mg$cd45_high,], am, "mfi"),
     expression(alpha*"-Syn MFI in CD45"^"high"*" microglia"), NULL),
  mk("cd45lo", getv(mg[!mg$cd45_high,], am, "mfi"),
     expression(alpha*"-Syn MFI in CD45"^"low"*" microglia"), NULL),
  mk("tmemhi", getv(mg[mg$ug_bright,], am, "mfi"),
     expression(alpha*"-Syn MFI in TMEM119"^"bright"*" microglia"), NULL),
  mk("tmemlo", getv(mg[mg$ug_dim,], am, "mfi"),
     expression(alpha*"-Syn MFI in TMEM119"^"dim"*" microglia"), NULL),
  mk("ssched", getv(mg[mg$SSC_A>166867,], am, "mfi"),
     expression(alpha*"-Syn MFI in SSC-A"^"high"*" microglia"), NULL),
  mk("ssclo",  getv(mg[mg$SSC_A<=166867,], am, "mfi"),
     expression(alpha*"-Syn MFI in SSC-A"^"low"*" microglia"), NULL))
f1 <- plot_grid(plotlist=P1, nrow=2, ncol=4, labels=letters[1:7], label_size=13) +
      theme(plot.background=element_rect(fill="white", colour=NA))
save_png(file.path(outd,"FIG1_aSyn_microglia.png"), f1, width=11.8, height=6.8, dpi=1200, bg="white")

nmCD45neg <- nm[nm$CD45<=c45,]; nmLy <- nm[nm$CD45>c45 & nm$CD11b<=c11,]
P2 <- list(
  mk("nUG",  getv(nm, an, "mfi"), expression(alpha*"-Syn MFI in non-microglia"), NULL),
  mk("neg",  getv(nmCD45neg, an, "mfi"),
     expression(alpha*"-Syn MFI in CD45"^"-"*" cells"), NULL),
  mk("ly",   getv(nmLy, an, "mfi"),
     expression(alpha*"-Syn MFI in CD45"^"+"*"CD11b"^"-"*" cells"), NULL))
f2 <- plot_grid(plotlist=P2, nrow=1, labels=letters[1:3], label_size=13) +
      theme(plot.background=element_rect(fill="white", colour=NA))
save_png(file.path(outd,"FIG2_aSyn_nonmicroglia.png"), f2, width=8.9, height=3.5, dpi=1200, bg="white")

# ---------- FIG 3: CD206 (CHI MFI) ----------
f3 <- mk("cd206", getv(mg, NA, "cd206mfi"), expression("CD206 MFI in microglia"), NULL) +
      theme(plot.background=element_rect(fill="white", colour=NA))
save_png(file.path(outd,"FIG3_CD206.png"), f3, width=3.1, height=3.5, dpi=1200, bg="white")

S <- do.call(rbind, ST); write.csv(S, file.path(outd,"bar_stats.csv"), row.names=FALSE)
print(format(S, digits=3), row.names=TRUE)
