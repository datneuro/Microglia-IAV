# ==============================================================================
# 02_Temporal-Trajectories.R  —  Temporal precedence 
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(correlation); library(ggrepel) })
PROJ <- "."
DIR  <- file.path(PROJ,"CrossPlatform_aSyn"); OUT <- file.path(DIR,"output"); FIG <- file.path(DIR,"figures")
CODE <- "06_CrossPlatform_Integration"   # module scripts (repo-relative)
source(file.path(CODE,"_Labels.R"))
ANCHOR <- "WB_aSyn_25kDa"
m <- read_csv(file.path(OUT,"master_matched.csv"), show_col_types=FALSE) %>%
  mutate(timepoint=factor(timepoint,c("4h","8h","24h")), treatment=factor(treatment,c("Control","V1","V2","PolyIC")))
markers <- grep("^(ELISA_|FACS_|qPCR_)", names(m), value=TRUE)
pretty_mk <- function(x) x %>% sub("^ELISA_","",.) %>% sub("^FACS_","",.) %>% sub("^qPCR_","",.) %>%
  str_replace("IL_","IL-") %>% str_replace("_"," ")

# condition means (per timepoint × treatment)
cond <- m %>% group_by(timepoint, treatment) %>%
  summarize(across(all_of(c(ANCHOR, markers)), ~mean(.x, na.rm=TRUE)), .groups="drop")

# ---- early (4h/8h) marker pattern vs LATE αSyn-25kDa (24h), across treatments -
late <- cond %>% filter(timepoint=="24h") %>% transmute(treatment, late = .data[[ANCHOR]])
earlyLate <- map_dfr(markers, function(mk) map_dfr(c("4h","8h"), function(tp) {
  e <- cond %>% filter(timepoint==tp) %>% transmute(treatment, ev = .data[[mk]])
  j <- inner_join(e, late, by="treatment") %>% drop_na()
  if (nrow(j) < 4 || sd(j$ev)==0) return(NULL)
  cc <- as.data.frame(correlation(select(j, ev, late), method="pearson", bayesian=TRUE))
  tibble(marker=pretty_mk(mk), platform=sub("_.*","",mk), early_tp=tp,
         rho=cc$rho, CI_low=cc$CI_low, CI_high=cc$CI_high, pd=cc$pd, n=cc$n_Obs)
}))
write_csv(earlyLate, file.path(OUT,"tierB_earlyLate.csv"))
cat("== Tier B: early(4h/8h) → late αSyn-25kDa(24h), condition-level n=4 (EXPLORATORY) ==\n")
earlyLate %>% arrange(desc(abs(rho))) %>% slice(1:12) %>%
  transmute(marker, platform, early_tp, rho=round(rho,2),
            CI=sprintf("[%.2f,%.2f]",CI_low,CI_high), pd=round(pd,3)) %>%
  as.data.frame() %>% print(row.names=FALSE)

# ============================== FIGURES =======================================
TRT_COL <- c(Control="#7F8C8D", V1="#4285F4", V2="#C0392B", PolyIC="#27AE60")

# ---- (1) TRAJECTORY OVERLAY (z-scored): αSyn vs early immune/transcription ----
# pick αSyn + canonical early markers + top-2 Tier-A protein hits
sel <- c(ANCHOR, "qPCR_MX1","qPCR_RSAD2","qPCR_IRF7","ELISA_IL_1b","ELISA_TNFa")
traj <- cond %>% select(timepoint, treatment, any_of(sel)) %>%
  pivot_longer(-c(timepoint,treatment), names_to="feat0", values_to="val") %>%
  group_by(feat0) %>% mutate(z = as.numeric(scale(val))) %>% ungroup() %>%
  mutate(is_asyn = feat0 == ANCHOR,
         feat = ifelse(is_asyn, "25-kDa α-Synuclein (Western blot)", disp_marker(feat0)),
         treatment = factor(disp_treat(as.character(treatment)), levels=TRT_LEVELS)) %>%
  filter(treatment != "Control")    # show virus/PolyIC responses
lab_pts <- traj %>% filter(!is.na(z)) %>% group_by(treatment, feat) %>% slice_tail(n=1)

fig_traj <- ggplot(traj, aes(timepoint, z, group=feat)) +
  geom_line(data=filter(traj,!is_asyn), color="grey65", linewidth=0.6) +
  geom_point(data=filter(traj,!is_asyn), color="grey65", size=1.4) +
  geom_line(data=filter(traj, is_asyn), color="#8E44AD", linewidth=1.7) +
  geom_point(data=filter(traj, is_asyn), color="#8E44AD", size=3) +
  ggrepel::geom_text_repel(data=lab_pts, aes(label=feat, color=is_asyn), size=2.6,
                           direction="y", hjust=0, nudge_x=0.15, segment.size=0.2,
                           show.legend=FALSE) +
  facet_wrap(~treatment, nrow=1) +
  scale_color_manual(values=c("FALSE"="grey45","TRUE"="#8E44AD"), guide="none") +
  scale_x_discrete(expand=expansion(mult=c(0.05,0.4))) +
  labs(x=NULL, y="z-scored level (across timecourse)",
       title="Temporal precedence: early immune / transcriptional signals vs late 25-kDa α-Synuclein",
       subtitle="Descriptive (between-subjects timepoints) — early markers peak at 4–8 h, 25-kDa α-Synuclein peaks at 24 h. NOT causal mediation.") +
  theme_bw(base_size=11) +
  theme(strip.background=element_rect(fill="grey85"), strip.text=element_text(face="bold"),
        plot.title=element_text(face="bold", size=11), plot.subtitle=element_text(size=8, color="grey35"))
ggsave(file.path(FIG,"tierB_trajectory.pdf"), fig_traj, width=12, height=4.5, device=cairo_pdf)
ggsave(file.path(FIG,"tierB_trajectory.png"), fig_traj, width=12, height=4.5, dpi=300)


message("\n[02_Temporal-Trajectories] SUPPLEMENT (descriptive only): wrote tierB_earlyLate.csv + trajectory")
