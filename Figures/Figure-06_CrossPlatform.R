
suppressPackageStartupMessages({
  library(tidyverse); library(ggdist); library(ggtext); library(ggrepel); library(patchwork)
})
OUT <- Sys.getenv("CROSSPLATFORM_OUT", unset = file.path("CrossPlatform_aSyn", "output"))  # CSVs of 06_CrossPlatform_Integration/
FIG <- Sys.getenv("FIGURE_OUT",        unset = file.path("Figures", "output"))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
source(file.path("06_CrossPlatform_Integration","_Labels.R"))
ANCHOR <- "WB_aSyn_25kDa"; platform_of <- function(x) sub("_.*","",x)
m <- read_csv(file.path(OUT,"master_matched.csv"), show_col_types=FALSE) %>%
  mutate(timepoint=factor(timepoint,c("4h","8h","24h")), treatment=factor(treatment))

# ---- generic half-eye with INLINE stats (text on inner side of each ridge) ---
half_eye <- function(drw, sumd, pal, fill_name, xtitle, ylabs) {
  pos <- sumd %>% filter(rho >= 0); neg <- sumd %>% filter(rho < 0)
  ggplot(drw, aes(draw, row, fill=fillv)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey60", linewidth=0.4) +
    stat_halfeye(point_interval=median_hdi, .width=0.95, slab_alpha=0.60,
                 height=2.0, normalize="xy", point_size=1.2) +
    geom_richtext(data=pos, aes(x=CI_low-0.11, label=stat), hjust=1, size=4, #2.55
                  fill=NA, label.color=NA, label.padding=unit(0,"pt"), color="grey12") +
    geom_richtext(data=neg, aes(x=CI_high+0.11, label=stat), hjust=0, size=4, #2.55
                  fill=NA, label.color=NA, label.padding=unit(0,"pt"), color="grey12") +
    scale_fill_manual(values=pal, name=fill_name, drop=FALSE) +
    scale_y_discrete(labels=ylabs) +
    scale_x_continuous(limits=c(-1.95,1.95), breaks=seq(-1,1,0.5)) +
    labs(x=xtitle, y=NULL) +
    theme_bw(base_size=11) +
    theme(panel.grid.minor=element_blank(), legend.position="top",
          axis.text.y=element_text(size=15, face = "bold"), axis.title.x=element_text(size=15, face="bold"),
          plot.margin=margin(4,6,4,5))
}

# ============================ (A) TIER A =======================================
res   <- read_csv(file.path(OUT,"tierA_correlations.csv"), show_col_types=FALSE)
draws <- read_csv(file.path(OUT,"tierA_rescor_draws.csv"), show_col_types=FALSE)
ordA  <- res %>% filter(model=="within") %>% arrange(rho) %>% pull(marker)
ylabA <- setNames(disp_marker(ordA), ordA)
sumA  <- res %>% filter(model=="within") %>%
  mutate(row=factor(marker,levels=ordA),
         fillv=factor(disp_platform(platform_of(marker)),levels=PLAT_LEVELS),
         stat=stat_inline_md(rho,CI_low,CI_high,pd))
drwA  <- draws %>% mutate(row=factor(marker,levels=ordA),
                          fillv=factor(disp_platform(platform_of(as.character(marker))),levels=PLAT_LEVELS))
pA_eye <- half_eye(drwA, sumA, PLAT_COL, "Platform", "Within-Condition (partial)", ylabA)

hd <- res %>% filter(model%in%c("total","np_adj")) %>%
  mutate(model=factor(model,c("total","np_adj"),c("Bivariate","NP-adjusted")),
         row=factor(marker,levels=ordA), star=pd_stars(pd))
pA_heat <- ggplot(hd, aes(model, row, fill=rho)) + geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f%s",rho,star)), size=4) + #2.5
  scale_fill_gradient2(low="#5B8AB5",mid="white",high="#C56B6B",midpoint=0,limits=c(-1,1),guide="none") +
  labs(x=NULL,y=NULL) +
  theme_minimal(base_size=9) +
  theme(axis.text.y=element_blank(), panel.grid=element_blank(),
        axis.text.x=element_text(face="bold",size=15), plot.margin=margin(4,4,4,2)) #size 8.3
pA <- (pA_eye | pA_heat) + plot_layout(widths=c(4,0.9)) +
  plot_annotation(caption="Posterior Estimates of Correlation Coefficient (ρ) with 25-kDa α-Synuclein",
                  theme=theme(plot.margin = margin(b = 45), plot.caption=element_text(hjust=0.7,size=15,face="bold",color="black",margin=margin(t=4)))) #11.5
ggsave(file.path(FIG,"tierA_main.pdf"), pA, width=13, height=8.6, device=cairo_pdf)
ggsave(file.path(FIG,"tierA_main.png"), pA, width=13, height=8.6, dpi=1000)

# ============================ (B) WB-INTERNAL ==================================
resW <- read_csv(file.path(OUT,"wb_internal_corr.csv"), show_col_types=FALSE)
drwW <- read_csv(file.path(OUT,"wb_internal_rescor_draws.csv"), show_col_types=FALSE)
axis_of <- function(l) case_when(
  grepl("PSMB|polyub|proteas",l,ignore.case=TRUE)~"Proteostasis", grepl("IκB|p65|NF-κB",l)~"NF-κB",
  grepl("LC3|Cathepsin",l)~"Autophagy–lysosome", grepl("Synuclein",l)~"α-Synuclein species",
  grepl("ASC",l)~"Inflammasome", TRUE~"Other")
lab_lk <- resW %>% filter(model=="within") %>% distinct(target,label) %>% deframe()
ordW   <- resW %>% filter(model=="within") %>% arrange(rho) %>% pull(target)
sumW   <- resW %>% filter(model=="within") %>%
  mutate(row=factor(target,levels=ordW), fillv=factor(axis_of(label)),
         stat=stat_inline_md(rho,CI_low,CI_high,pd))
drwW   <- drwW %>% mutate(row=factor(target,levels=ordW), fillv=factor(axis_of(lab_lk[as.character(target)])))
pB <- half_eye(drwW, sumW, AXIS_COL, "Pathway axis",
               "Bayesian Intra-Lysate Partial Correlation (ρ) with 25-kDa α-Synuclein",
               setNames(lab_lk[ordW], ordW))
ggsave(file.path(FIG,"wb_internal_forest.pdf"), pB, width=12, height=6.6, device=cairo_pdf)
ggsave(file.path(FIG,"wb_internal_forest.png"), pB, width=12, height=6.6, dpi=1000)

# ============================ (C) TIER C BIPLOT ================================
# read Tier-C PCA outputs from 04 (single source — no recompute) ; WB α-Syn = purple
sc  <- read_csv(file.path(OUT,"tierC_pca_scores_early_qPCR.csv"), show_col_types=FALSE) %>%
  mutate(treatment=factor(treatment,levels=TRT_LEVELS), timepoint=factor(timepoint))
ld0 <- read_csv(file.path(OUT,"tierC_pca_loadings_early_qPCR.csv"), show_col_types=FALSE) %>%
  mutate(is_wb=grepl("^WB_",marker), is_asyn=marker==ANCHOR)
pmv <- read_csv(file.path(OUT,"tierC_permanova_early_qPCR.csv"), show_col_types=FALSE)
pcv <- read_csv(file.path(OUT,"tierC_var_early_qPCR.csv"), show_col_types=FALSE)$pct
getR <- function(t) pmv$R2[pmv$term==t]; getP <- function(t) pmv$`Pr(>F)`[pmv$term==t]
amp <- 0.9*max(abs(c(sc$PC1,sc$PC2))) / max(abs(c(ld0$PC1,ld0$PC2)))
ld  <- ld0 %>% mutate(x=PC1*amp, y=PC2*amp)
pC <- ggplot(sc, aes(PC1,PC2)) +
  stat_ellipse(aes(color=treatment), level=0.68, linewidth=0.5) +
  geom_point(aes(color=treatment, shape=timepoint), size=2.6) +
  geom_segment(data=ld, aes(x=0,y=0,xend=x,yend=y),
               arrow=arrow(length=unit(0.15,"cm")),
               color=ifelse(ld$is_wb,"#8E44AD","grey60"),
               alpha=ifelse(ld$is_wb,0.95,0.30),
               linewidth=ifelse(ld$is_asyn,1.1, ifelse(ld$is_wb,0.8,0.4))) +
  geom_text_repel(data=ld, aes(x,y,label=label, fontface=ifelse(is_asyn,"bold","plain")),
                  color=ifelse(ld$is_wb,"#8E44AD","grey25"), size=3, max.overlaps=30, #size 2.6
                  bg.color="white", bg.r=0.14, segment.color=NA) +
  scale_color_manual(values=TRT_COL, name="Treatment") +
  scale_shape_discrete(name="Time point") +
  labs(x=sprintf("PC1 (%.1f%%)",pcv[1]), y=sprintf("PC2 (%.1f%%)",pcv[2]),
       subtitle=sprintf("PERMANOVA Treatment R²=%.2f (p=%.4f), Time point R²=%.2f (p=%.4f)",
                        getR("treatment"),getP("treatment"),getR("timepoint"),getP("timepoint"))) +
  theme_bw(base_size=11) + theme(plot.subtitle=element_text(face="bold",size=9,color="grey20"),  axis.title = element_text(face="bold", size=15))
ggsave(file.path(FIG,"tierC_biplot_early_qPCR.pdf"), pC, width=8, height=6.5, device=cairo_pdf)
ggsave(file.path(FIG,"tierC_biplot_early_qPCR.png"), pC, width=8, height=6.5, dpi=1000)

# ============================ SUPP: scatter + trajectory ======================
resid_rhs <- function(d,v,rhs) residuals(lm(as.formula(paste0("`",v,"` ~ ",rhs)),data=d))
top3 <- res %>% filter(model=="within") %>% arrange(desc(abs(rho))) %>% slice(1:3) %>% pull(marker)
scat <- map_dfr(top3, function(mk){ d<-m%>%select(all_of(c(ANCHOR,mk)),timepoint,treatment)%>%drop_na()
  d$ra<-resid_rhs(d,ANCHOR,"timepoint*treatment"); d$rm<-resid_rhs(d,mk,"timepoint*treatment")
  transmute(d, marker=disp_marker(mk), treatment=factor(disp_treat(as.character(treatment)),levels=TRT_LEVELS), timepoint, ra, rm)}) %>%
  mutate(marker=factor(marker,levels=disp_marker(top3)))
pScat <- ggplot(scat,aes(ra,rm)) + geom_smooth(method="lm",se=TRUE,color="grey30",fill="grey80",linewidth=0.7) +
  geom_point(aes(color=treatment, shape=timepoint), size=2.6) + facet_wrap(~marker,scales="free",nrow=1) +
  scale_color_manual(values=TRT_COL, name="Treatment") +
  labs(x="25-kDa α-Synuclein (within-condition residual)", y="Marker (within-condition residual)") +
  theme_bw(base_size=10)+theme(strip.text=element_text(face="bold"), strip.background=element_rect(fill="grey88"))

cond <- m %>% group_by(timepoint, treatment) %>%
  summarize(across(c(starts_with("WB_"), starts_with("ELISA_"), starts_with("qPCR_")), ~mean(.x,na.rm=TRUE)), .groups="drop")
sel <- c(ANCHOR,"WB_aSyn_multimer","qPCR_MX1","qPCR_RSAD2","qPCR_IRF7","ELISA_IL_1b","ELISA_TNFa")
traj <- cond %>% select(timepoint, treatment, any_of(sel)) %>%
  pivot_longer(-c(timepoint,treatment), names_to="f0", values_to="val") %>%
  group_by(f0) %>% mutate(z=as.numeric(scale(val))) %>% ungroup() %>%
  mutate(is_asyn=f0==ANCHOR, is_wb=grepl("^WB_",f0),
         feat=case_when(f0==ANCHOR ~ "25-kDa α-Synuclein (Western blot)",
                        f0=="WB_aSyn_multimer" ~ "α-Synuclein multimer (Western blot)",
                        TRUE ~ disp_marker(f0)),
         treatment=factor(disp_treat(as.character(treatment)),levels=TRT_LEVELS)) %>%
  filter(treatment!="Control")
lab_pts <- traj %>% filter(!is.na(z)) %>% group_by(treatment,feat) %>% slice_tail(n=1)
pTraj <- ggplot(traj, aes(timepoint,z,group=feat)) +
  geom_line(data=filter(traj,!is_wb),color="grey65",linewidth=0.6) +
  geom_point(data=filter(traj,!is_wb),color="grey65",size=1.4) +
  geom_line(data=filter(traj,is_wb&!is_asyn),color="#8E44AD",linewidth=1.1,linetype="22") +
  geom_point(data=filter(traj,is_wb&!is_asyn),color="#8E44AD",size=2) +
  geom_line(data=filter(traj,is_asyn),color="#8E44AD",linewidth=1.8) +
  geom_point(data=filter(traj,is_asyn),color="#8E44AD",size=3) +
  geom_text_repel(data=lab_pts, aes(label=feat,color=is_wb), size=2.6, direction="y", hjust=0,
                  nudge_x=0.15, segment.size=0.2, show.legend=FALSE) +
  facet_wrap(~treatment,nrow=1) +
  scale_color_manual(values=c("FALSE"="grey45","TRUE"="#8E44AD"),guide="none") +
  scale_x_discrete(expand=expansion(mult=c(0.05,0.4))) +
  labs(x=NULL,y="z-scored level (across time course)") +
  theme_bw(base_size=11) + theme(strip.background=element_rect(fill="grey88"),strip.text=element_text(face="bold"))
ggsave(file.path(FIG,"tierA_scatter_top.pdf"), pScat, width=11, height=4, device=cairo_pdf)
ggsave(file.path(FIG,"tierA_scatter_top.png"), pScat, width=11, height=4, dpi=1000)
ggsave(file.path(FIG,"tierB_trajectory.pdf"), pTraj, width=12, height=4.3, device=cairo_pdf)
ggsave(file.path(FIG,"tierB_trajectory.png"), pTraj, width=12, height=4.3, dpi=1000)

# ============================ MERGED COMPOSITES ===============================
bc   <- (wrap_elements(pB) | wrap_elements(pC)) + plot_layout(widths=c(1.35,1))
main <- wrap_elements(pA) / bc +                                         
      plot_annotation(tag_levels="a") + plot_layout(heights=c(1.5,1.0)) &    
      theme(plot.tag = element_text(size = 30, face = "bold"))     
ggsave(file.path(FIG,"MAIN_crossplatform.pdf"), main, width=17, height=15, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(FIG,"MAIN_crossplatform.png"), main, width=17, height=15, dpi=500, limitsize=FALSE)

supp <- wrap_elements(pScat) / wrap_elements(pTraj) +                    
      plot_annotation(tag_levels="a") + plot_layout(heights=c(1,1)) &        
      theme(plot.tag = element_text(size = 30, face = "bold"))   
ggsave(file.path(FIG,"SUPP_crossplatform.pdf"), supp, width=12, height=8.5, device=cairo_pdf)
ggsave(file.path(FIG,"SUPP_crossplatform.png"), supp, width=12, height=8.5, dpi=600)
message("[05_figures] individual + MAIN + SUPP composites done")
