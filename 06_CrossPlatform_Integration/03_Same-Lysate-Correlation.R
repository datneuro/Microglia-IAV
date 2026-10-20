# ==============================================================================
# 03_Same-Lysate-Correlation.R 
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(brms); library(bayestestR); library(correlation)
  library(ggdist); library(ggtext); library(patchwork)
})
PROJ <- "."
DIR  <- file.path(PROJ,"CrossPlatform_aSyn"); OUT <- file.path(DIR,"output"); FIG <- file.path(DIR,"figures")
CODE <- "06_CrossPlatform_Integration"   # module scripts (repo-relative)
source(file.path("04_Western_Blot","00_Config-and-Ingestion.R")); source(file.path(CODE,"_Labels.R"))
source(file.path(CODE,"_Bayes-Common.R"))   # MCMC config, ROPE, QC gate
`%||%` <- function(a,b) if(is.null(a)||(is.data.frame(a)&&nrow(a)==0)) b else a
ANCHOR <- "syn1.25"; SEED <- CP_SEED

p1 <- ingest_phase1()
wide <- p1 %>% select(id,group,target,y) %>% distinct(id,target,.keep_all=TRUE) %>%
  pivot_wider(names_from=target,values_from=y) %>% left_join(distinct(p1,id,group),by=c("id","group"))
drop_raw <- c("pikb","ikb","pP65","P65","lc3.1","lc3.2")
markers <- setdiff(names(wide), c("id","group",ANCHOR,drop_raw))

zd <- function(mk) wide %>% transmute(group, a=as.numeric(scale(.data[[ANCHOR]])), b=as.numeric(scale(.data[[mk]]))) %>% drop_na()
message("[rescor] compiling base WB-internal model …")
# n=18 small sample: iter = 8000 (4 x 7000 = 28k draws) so that the rescor ESS
# clears the >10k bulk gate for all targets (Tier A uses the default 6000).
base_fit <- cp_brm(bf(mvbind(a,b) ~ group) + set_rescor(TRUE), data=zd(markers[1]), iter=8000)
rescor_fit <- function(mk){ d<-zd(mk); if(nrow(d)<10||sd(d$b)==0) return(NULL)
  fit<-tryCatch(update(base_fit,newdata=d,refresh=0,silent=2),error=function(e) NULL); if(is.null(fit)) return(NULL)
  dr<-as.numeric(as_draws_df(fit)[["rescor__a__b"]])
  summ<-cp_describe_rho(dr); summ$n<-nrow(d)
  list(summary=summ, draws=dr, qc=cp_mcmc_qc(fit, paste0("wb_internal:", mk))) }
bcor <- function(x,y){ d<-tibble(x,y)%>%drop_na(); if(nrow(d)<8||sd(d$x)==0||sd(d$y)==0) return(NULL)
  cc<-as.data.frame(correlation(d,method="pearson",bayesian=TRUE))
  tibble(rho=cc$rho,CI_low=cc$CI_low,CI_high=cc$CI_high,pd=cc$pd,BF=cc$BF,n=cc$n_Obs) }
resid_grp <- function(d,v) residuals(lm(reformulate("group",v),data=d))

res_list<-list(); draws_list<-list(); qc_list<-list()
for (mk in markers) {
  d <- wide %>% select(all_of(c(ANCHOR,mk)),group) %>% drop_na(); if(nrow(d)<10) next
  rf <- rescor_fit(mk); if(!is.null(rf)) { draws_list[[mk]] <- tibble(target=mk, draw=rf$draws); qc_list[[mk]] <- rf$qc }
  d2<-d; names(d2)[1:2]<-c("A","M")
  rows <- bind_rows(mutate((if(!is.null(rf)) rf$summary else tibble()), model="within"),
            mutate(bcor(resid_grp(d2,"A"),resid_grp(d2,"M")) %||% tibble(), model="within_2step"),
            mutate(bcor(d[[ANCHOR]],d[[mk]]) %||% tibble(), model="total")) %>%
    mutate(target=mk, label=lab_target(mk), .before=1)
  res_list[[mk]] <- rows
}
res <- bind_rows(res_list); draws <- bind_rows(draws_list); qc <- bind_rows(qc_list)
write_csv(res, file.path(OUT,"wb_internal_corr.csv"))
write_csv(draws, file.path(OUT,"wb_internal_rescor_draws.csv"))
write_csv(qc, file.path(OUT,"wb_internal_qc.csv"))
cat("\n== WB-internal MCMC QC (rescor estimand) ==\n")
qc %>% transmute(target=sub("wb_internal:","",model), divergent, max_rhat=round(max_rhat,3),
                 rescor_ess_bulk=round(rescor_ess_bulk), rescor_ess_tail=round(rescor_ess_tail),
                 min_ess_all=round(min_ess_bulk), pass) %>% as.data.frame() %>% print(row.names=FALSE)
cat(sprintf("[QC] divergences total=%d | all rescor ESS>10k: %s | max Rhat=%.3f\n",
            sum(qc$divergent), all(qc$rescor_ess_bulk>=CP_ESS_MIN & qc$rescor_ess_tail>=CP_ESS_MIN),
            max(qc$max_rhat)))
cat("== WB same-lysate within (rescor) — |rho| ranked ==\n")
res %>% filter(model=="within") %>% arrange(desc(abs(rho))) %>%
  transmute(label,rho=round(rho,2),CI=sprintf("[%.2f,%.2f]",CI_low,CI_high),pd=round(pd,3)) %>%
  as.data.frame() %>% print(row.names=FALSE)

# ============================== FIGURE (half-eye) =============================
axis_of <- function(l) case_when(
  grepl("PSMB|polyub|proteas",l,ignore.case=TRUE)~"Proteostasis", grepl("IκB|p65|NF-κB",l)~"NF-κB",
  grepl("LC3|Cathepsin",l)~"Autophagy–lysosome", grepl("Synuclein",l)~"α-Synuclein species",
  grepl("ASC",l)~"Inflammasome", TRUE~"Other")
AX <- c("Proteostasis"="#16A085","NF-κB"="#C0392B","Autophagy–lysosome"="#2980B9",
        "α-Synuclein species"="#8E44AD","Inflammasome"="#E67E22","Other"="grey50")
ord  <- res %>% filter(model=="within") %>% arrange(rho) %>% pull(target)
ylab <- setNames(lab_target(ord), ord)
sumW <- res %>% filter(model=="within") %>%
  left_join(filter(res,model=="within_2step")%>%select(target,BF2=BF),by="target") %>%
  mutate(target=factor(target,levels=ord), axis=axis_of(label), stat=stat_md(rho,CI_low,CI_high,pd,BF2))
drw <- draws %>% mutate(target=factor(target,levels=ord),
                        axis=axis_of(lab_target(as.character(target))))

p_eye <- ggplot(drw, aes(draw, target, fill=axis)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey55") +
  stat_halfeye(point_interval=median_hdi, .width=0.95, slab_alpha=0.75,
               height=1.7, normalize="xy", point_size=1.3) +
  scale_fill_manual(values=AX, name="Pathway axis") +
  scale_y_discrete(labels=ylab) + scale_x_continuous(limits=c(-1,1), breaks=seq(-1,1,0.5)) +
  labs(x="Posterior Estimates of Partial Correlation (ρ) with 25-kDa α-Synuclein species\n(same lysate, n = 18)", y=NULL) +
  theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="top",
        axis.text.y=element_text(size=9), plot.margin=margin(5,2,5,5))
p_stat <- ggplot(sumW, aes(0,target)) +
  geom_richtext(aes(label=stat), hjust=0, size=2.7, fill=NA, label.color=NA, label.padding=unit(0,"pt")) +
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) + theme_void() + theme(plot.margin=margin(5,5,5,2))
fig <- p_eye + p_stat + plot_layout(widths=c(3.2, 2.7))
ggsave(file.path(FIG,"wb_internal_forest.pdf"),fig,width=12,height=6.5,device=cairo_pdf)
ggsave(file.path(FIG,"wb_internal_forest.png"),fig,width=12,height=6.5,dpi=300)
message("\n[03_Same-Lysate-Correlation] half-eye figure done")
