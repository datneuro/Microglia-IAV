# ==============================================================================
# 05_Sensitivity-Analysis.R  —
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(brms); library(bayestestR) })
PROJ <- "."
DIR  <- file.path(PROJ,"CrossPlatform_aSyn"); OUT <- file.path(DIR,"output")
CODE <- "06_CrossPlatform_Integration"  
source(file.path("04_Western_Blot","00_Config-and-Ingestion.R")) 
source(file.path(CODE,"_Labels.R")); source(file.path(CODE,"_Bayes-Common.R"))

rope_pct <- function(dr,w) as.numeric(bayestestR::rope(dr,range=c(-w,w),ci=1)$ROPE_Percentage)
classify <- function(rp,lo,hi,w){ rp<-ifelse(!is.na(rp)&rp>1,rp/100,rp)
  dplyr::case_when((lo> w)|(hi< -w) ~ "reject", (!is.na(rp)&rp>=0.95) ~ "accept", TRUE ~ "undecided") }

describe_sens <- function(fit){
  dr <- as.numeric(as_draws_df(fit)[["rescor__a__b"]]); h <- bayestestR::hdi(dr,ci=.95)
  tibble(median=median(dr), HDI_low=as.numeric(h$CI_low), HDI_high=as.numeric(h$CI_high),
         pd=as.numeric(bayestestR::p_direction(dr)),
         ROPE_pct_05=rope_pct(dr,0.05), ROPE_pct_10=rope_pct(dr,0.10), ROPE_pct_15=rope_pct(dr,0.15)) }

run_tier <- function(tier, markers, zfun, fml){
  out <- list()
  for (pn in names(CP_PRIOR_SETS)) {
    message("  [", tier, "] prior=", pn, " — compiling base …")
    base <- cp_brm(fml, data=zfun(markers[1]), prior=CP_PRIOR_SETS[[pn]], iter=CP_ITER_SENS)
    for (mk in markers) {
      d <- zfun(mk); if (nrow(d)<10 || sd(d$b)==0) next
      fit <- tryCatch(update(base, newdata=d, refresh=0, silent=2), error=function(e) NULL)
      if (is.null(fit)) next
      r <- describe_sens(fit); r$tier<-tier; r$marker<-mk; r$prior_set<-pn
      out[[paste(tier,pn,mk,sep="|")]] <- r
    }
  }
  bind_rows(out)
}

ANCHOR_A <- "WB_aSyn_25kDa"
m <- read_csv(file.path(OUT,"master_matched.csv"), show_col_types=FALSE) %>%
  mutate(timepoint=factor(timepoint,c("4h","8h","24h")), treatment=factor(treatment))
markers_A <- setdiff(grep("^(ELISA_|FACS_|qPCR_|WB_)",names(m),value=TRUE), ANCHOR_A)
zA <- function(mk) m %>% transmute(timepoint,treatment,
  a=as.numeric(scale(.data[[ANCHOR_A]])), b=as.numeric(scale(.data[[mk]]))) %>% drop_na()

ANCHOR_W <- "syn1.25"
p1 <- ingest_phase1()
wide <- p1 %>% select(id,group,target,y) %>% distinct(id,target,.keep_all=TRUE) %>%
  pivot_wider(names_from=target,values_from=y) %>% left_join(distinct(p1,id,group),by=c("id","group"))
drop_raw <- c("pikb","ikb","pP65","P65","lc3.1","lc3.2")
markers_W <- setdiff(names(wide), c("id","group",ANCHOR_W,drop_raw))
zW <- function(mk) wide %>% transmute(group, a=as.numeric(scale(.data[[ANCHOR_W]])),
  b=as.numeric(scale(.data[[mk]]))) %>% drop_na()

# ---- run ---------------------------------------------------------------------
sens <- bind_rows(
  run_tier("tierA",       markers_A, zA, bf(mvbind(a,b) ~ timepoint*treatment) + set_rescor(TRUE)),
  run_tier("wb_internal", markers_W, zW, bf(mvbind(a,b) ~ group)               + set_rescor(TRUE))
)

sens <- sens %>% mutate(
  class_05 = if_else(prior_set=="default", classify(ROPE_pct_05,HDI_low,HDI_high,0.05), NA_character_),
  class_10 = classify(ROPE_pct_10,HDI_low,HDI_high,0.10),
  class_15 = if_else(prior_set=="default", classify(ROPE_pct_15,HDI_low,HDI_high,0.15), NA_character_))

rope_stab <- sens %>% filter(prior_set=="default") %>%
  transmute(tier,marker, rope_stable=(class_05==class_10)&(class_10==class_15))
prior_stab <- sens %>% group_by(tier,marker) %>%
  summarize(prior_stable=(n_distinct(sign(median))==1)&(n_distinct(class_10)==1), .groups="drop")

sens_final <- sens %>%
  left_join(rope_stab, by=c("tier","marker")) %>%
  left_join(prior_stab, by=c("tier","marker")) %>%
  mutate(label = ifelse(tier=="tierA", disp_marker(marker), lab_target(marker))) %>%
  arrange(tier, marker, factor(prior_set, levels=c("default","broader","narrower"))) %>%
  select(tier,marker,label,prior_set,median,HDI_low,HDI_high,pd,
         ROPE_pct_05,ROPE_pct_10,ROPE_pct_15,class_05,class_10,class_15,rope_stable,prior_stable)

write_csv(sens_final, file.path(OUT,"crossplatform_sensitivity.csv"))
cat("\n== Cross-platform sensitivity (default prior, ROPE |0.1|) ==\n")
sens_final %>% filter(prior_set=="default") %>%
  transmute(tier,label,rho=round(median,2),CI=sprintf("[%.2f,%.2f]",HDI_low,HDI_high),
            class_10, rope_stable, prior_stable) %>% as.data.frame() %>% print(row.names=FALSE)
n_stable <- sens_final %>% filter(prior_set=="default") %>% summarize(s=sum(rope_stable&prior_stable),n=n())
cat(sprintf("\n[sensitivity] %d/%d markers stable on BOTH prior and ROPE-width. Saved crossplatform_sensitivity.csv\n",
            n_stable$s, n_stable$n))
message("\n[05_Sensitivity-Analysis] DONE.")
