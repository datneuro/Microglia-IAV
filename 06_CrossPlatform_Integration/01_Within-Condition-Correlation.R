# ==============================================================================
# 01_Within-Condition-Correlation.R  —  Within-condition Bayesian PARTIAL correlation
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(brms); library(bayestestR); library(correlation)
  library(ggdist); library(ggtext); library(patchwork)
})
PROJ <- "."
DIR  <- file.path(PROJ,"CrossPlatform_aSyn"); OUT <- file.path(DIR,"output"); FIG <- file.path(DIR,"figures")
CODE <- "06_CrossPlatform_Integration"  
source(file.path(CODE,"_Labels.R"))
source(file.path(CODE,"_Bayes-Common.R"))   
ANCHOR <- "WB_aSyn_25kDa"; SEED <- CP_SEED
`%||%` <- function(a,b) if(is.null(a)||(is.data.frame(a)&&nrow(a)==0)) b else a

m <- read_csv(file.path(OUT,"master_matched.csv"), show_col_types=FALSE) %>%
  mutate(timepoint=factor(timepoint,c("4h","8h","24h")), treatment=factor(treatment))
# include the other WB α-Syn species (multimer, phospho) as markers; exclude the anchor
markers <- setdiff(grep("^(ELISA_|FACS_|qPCR_|WB_)", names(m), value=TRUE), ANCHOR)
platform_of <- function(x) sub("_.*","",x)

zdat <- function(mk) m %>% transmute(timepoint, treatment,
                                     a=as.numeric(scale(.data[[ANCHOR]])), b=as.numeric(scale(.data[[mk]]))) %>% drop_na()
message("[rescor] compiling base multivariate model …")
base_fit <- cp_brm(bf(mvbind(a,b) ~ timepoint*treatment) + set_rescor(TRUE),
                   data=zdat(markers[1]))   
rescor_fit <- function(mk){ d<-zdat(mk); if(nrow(d)<12||sd(d$b)==0) return(NULL)
  fit<-tryCatch(update(base_fit,newdata=d,refresh=0,silent=2),error=function(e) NULL); if(is.null(fit)) return(NULL)
  dr<-as.numeric(as_draws_df(fit)[["rescor__a__b"]])
  summ<-cp_describe_rho(dr); summ$n<-nrow(d)                    # median, HDI, pd, ROPE%, BF
  list(summary=summ, draws=dr, qc=cp_mcmc_qc(fit, paste0("tierA:", mk))) }

bcor <- function(x,y){ d<-tibble(x,y)%>%drop_na(); if(nrow(d)<10||sd(d$x)==0||sd(d$y)==0) return(NULL)
  cc<-as.data.frame(correlation(d,method="pearson",bayesian=TRUE))
  tibble(rho=cc$rho,CI_low=cc$CI_low,CI_high=cc$CI_high,pd=cc$pd,BF=cc$BF,n=cc$n_Obs) }
resid_rhs <- function(d,v,rhs) residuals(lm(as.formula(paste0("`",v,"` ~ ",rhs)),data=d))
pcor <- function(df,mk,rhs){ vars<-c(ANCHOR,mk,all.vars(as.formula(paste0("~",rhs))))
  d<-df%>%select(any_of(vars))%>%drop_na(); if(nrow(d)<12||n_distinct(d[[mk]])<3) return(NULL)
  bcor(resid_rhs(d,ANCHOR,rhs),resid_rhs(d,mk,rhs)) }

# ---- run ---------------------------------------------------------------------
res_list <- list(); draws_list <- list(); qc_list <- list()
for (mk in markers) {
  is_q <- platform_of(mk)=="qPCR"
  rf <- rescor_fit(mk)
  if (!is.null(rf)) { draws_list[[mk]] <- tibble(marker=mk, draw=rf$draws); qc_list[[mk]] <- rf$qc }
  rows <- bind_rows(
    mutate((if(!is.null(rf)) rf$summary else tibble()), model="within"),
    mutate(pcor(m,mk,"timepoint*treatment")          %||% tibble(), model="within_2step"),
    mutate(bcor(m[[ANCHOR]],m[[mk]])                  %||% tibble(), model="total"),
    if(!is_q) mutate(pcor(filter(m,timepoint=="24h"),mk,"treatment") %||% tibble(), model="t24"),
    if(mk!="FACS_NP") mutate(pcor(m,mk,"timepoint*treatment + FACS_NP") %||% tibble(), model="np_adj")
  ) %>% mutate(marker=mk, platform=platform_of(mk), label=disp_marker(mk), .before=1)
  res_list[[mk]] <- rows
}
res <- bind_rows(res_list); draws <- bind_rows(draws_list); qc <- bind_rows(qc_list)
write_csv(res, file.path(OUT,"tierA_correlations.csv"))
write_csv(draws, file.path(OUT,"tierA_rescor_draws.csv"))
write_csv(qc, file.path(OUT,"tierA_qc.csv"))
cat("\n== Tier A MCMC QC (rescor estimand) ==\n")
qc %>% transmute(marker=sub("tierA:","",model), divergent, max_rhat=round(max_rhat,3),
                 rescor_ess_bulk=round(rescor_ess_bulk), rescor_ess_tail=round(rescor_ess_tail),
                 min_ess_all=round(min_ess_bulk), pass) %>% as.data.frame() %>% print(row.names=FALSE)
cat(sprintf("[QC] divergences total=%d | all rescor ESS>10k: %s | max Rhat=%.3f\n",
            sum(qc$divergent), all(qc$rescor_ess_bulk>=CP_ESS_MIN & qc$rescor_ess_tail>=CP_ESS_MIN),
            max(qc$max_rhat)))
cat("\n== Tier A within (rescor) — |rho| ranked ==\n")
res %>% filter(model=="within") %>% arrange(desc(abs(rho))) %>%
  transmute(label,platform,rho=round(rho,2),CI=sprintf("[%.2f,%.2f]",CI_low,CI_high),pd=round(pd,3)) %>%
  as.data.frame() %>% print(row.names=FALSE)

# ============================== MAIN FIGURE ===================================
ord  <- res %>% filter(model=="within") %>% arrange(rho) %>% pull(marker)
ylab <- setNames(disp_marker(ord), ord)
sumA <- res %>% filter(model=="within") %>%
  left_join(filter(res,model=="within_2step") %>% select(marker, BF2=BF), by="marker") %>%
  mutate(marker=factor(marker,levels=ord),
         stat=stat_md(rho,CI_low,CI_high,pd,BF2))
drw <- draws %>% mutate(marker=factor(marker,levels=ord),
                        platform=factor(disp_platform(platform_of(as.character(marker))),levels=PLAT_LEVELS))

p_eye <- ggplot(drw, aes(draw, marker, fill=platform)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey55") +
  stat_halfeye(point_interval=median_hdi, .width=0.95, slab_alpha=0.75,
               height=1.7, normalize="xy", point_size=1.3) +
  scale_fill_manual(values=PLAT_COL, name="Platform", drop=FALSE) +
  scale_y_discrete(labels=ylab) +
  scale_x_continuous(limits=c(-1,1), breaks=seq(-1,1,0.5)) +
  labs(x="Posterior Estimates of Partial Correlation (ρ) with 25-kDa α-Synuclein species", y=NULL) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        axis.text.y=element_text(size=9), plot.margin=margin(5,2,5,5))

hd <- res %>% filter(model%in%c("total","np_adj")) %>%
  mutate(model=factor(model,c("total","np_adj"),c("Total","NP-adj")), marker=factor(marker,levels=ord))
p_heat <- ggplot(hd, aes(model, marker, fill=rho)) + geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f",rho)), size=2.5) +
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-1,1),guide="none") +
  labs(x=NULL,y=NULL,title="Robustness") +
  theme_minimal(base_size=9) +
  theme(axis.text.y=element_blank(), panel.grid=element_blank(),
        axis.text.x=element_text(face="bold",size=8), plot.title=element_text(size=9,hjust=0.5),
        plot.margin=margin(5,2,5,2))

p_stat <- ggplot(sumA, aes(0, marker)) +
  geom_richtext(aes(label=stat), hjust=0, size=2.7, fill=NA, label.color=NA, label.padding=unit(0,"pt")) +
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_void() + theme(plot.margin=margin(5,5,5,2))

fig <- p_eye + p_heat + p_stat + plot_layout(widths=c(3.0, 0.7, 2.9))
ggsave(file.path(FIG,"tierA_main.pdf"), fig, width=13.5, height=8.5, device=cairo_pdf)
ggsave(file.path(FIG,"tierA_main.png"), fig, width=13.5, height=8.5, dpi=300)

# ---- SUPP scatter (top-3 within, residual space) — academic labels ----------
top3 <- res %>% filter(model=="within") %>% arrange(desc(abs(rho))) %>% slice(1:3) %>% pull(marker)
sc <- map_dfr(top3, function(mk){ d<-m%>%select(all_of(c(ANCHOR,mk)),timepoint,treatment)%>%drop_na()
  d$ra<-resid_rhs(d,ANCHOR,"timepoint*treatment"); d$rm<-resid_rhs(d,mk,"timepoint*treatment")
  transmute(d, marker=disp_marker(mk), treatment=disp_treat(as.character(treatment)), timepoint, ra, rm) }) %>%
  mutate(treatment=factor(treatment,levels=TRT_LEVELS), marker=factor(marker,levels=disp_marker(top3)))
fig_sc <- ggplot(sc,aes(ra,rm)) + geom_smooth(method="lm",se=TRUE,color="grey30",fill="grey80",linewidth=0.7) +
  geom_point(aes(color=treatment,shape=timepoint),size=2.6) + facet_wrap(~marker,scales="free",nrow=1) +
  scale_color_manual(values=TRT_COL, name="Treatment") +
  labs(x="25-kDa α-Synuclein (within-condition residual)", y="Marker (within-condition residual)",
       title="Top within-condition associations (individual dishes, residual space)") +
  theme_bw(base_size=10)+theme(plot.title=element_text(face="bold",size=11), strip.text=element_text(face="bold"))
ggsave(file.path(FIG,"tierA_scatter_top.pdf"),fig_sc,width=11,height=4,device=cairo_pdf)
ggsave(file.path(FIG,"tierA_scatter_top.png"),fig_sc,width=11,height=4,dpi=300)
message("\n[01_Within-Condition-Correlation] half-eye main figure + scatter done")
