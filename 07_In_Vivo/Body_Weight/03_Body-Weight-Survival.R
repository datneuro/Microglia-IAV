#!/usr/bin/env Rscript
# =============================================================================
# 03_Body-Weight-Survival.R
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(stringr); library(purrr); library(tibble); library(scales)
})
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)))
have <- function(p) requireNamespace(p, quietly = TRUE)

# ------------------------------ CONFIGURATION --------------------------------
xlsx_name   <- "20260826 Virus infection 5 days for HMC3.xlsx"
sheet_name  <- 1
out_fig     <- file.path("Body weight", "figures")   # Body weight/ = data folder (repository root)
out_tab     <- file.path("Body weight", "tables")

group_levels <- c("Naïve", "PR8", "S88")     # display order
ref_group    <- "Naïve"                       # reference group
group_cols   <- c("Naïve" = "#9E9E9E", "PR8" = "#4C72B0", "S88" = "#C44E52")
group_rename <- c("PR8" = "A/PR/8/34",      # display label of each group
                  "S33" = "A/NWS/33",
                  "S88" = "A/NWS/33")

err_type     <- "SEM"      
show_indiv   <- TRUE       

# --- SURVIVAL ----------------------------------------------------------------
death_events <- tribble(~Group, ~Mouse, ~time, ~status)

wl_endpoint  <- 20       
wl_sensitivity <- c(20, 10, 5)  
km_offset    <- 0.9     
last_day     <- NA         

dpi <- 400
set.seed(2026)


## ============================ 1. READ DATA ====================================
find_file <- function(f) {
  cand <- c(f, file.path("Body weight", f), file.path("..", f))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("File not found: ", f, call. = FALSE)
  hit[1]
}
read_excel_safe <- function(path, ...) {
  out <- try(read_excel(path, ...), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  tmp <- file.path(tempdir(), "bw_tmp.xlsx"); file.copy(path, tmp, overwrite = TRUE)
  read_excel(tmp, ...)
}

xlsx_path <- find_file(xlsx_name)
message("Reading: ", xlsx_path)
raw <- read_excel_safe(xlsx_path, sheet = sheet_name, col_names = FALSE,
                       .name_repair = "minimal")
raw <- as.data.frame(raw, stringsAsFactors = FALSE)

# locate the header rows containing "Day 0"
is_hdr <- apply(raw, 1, function(r) any(str_detect(str_squish(as.character(r)), "^Day ?0$"), na.rm = TRUE))
hdr_rows <- which(is_hdr)
if (!length(hdr_rows)) stop("No header row containing 'Day 0' found.", call. = FALSE)

parse_block <- function(h) {
  hdr  <- str_squish(as.character(unlist(raw[h, ])))
  dcol <- which(str_detect(hdr, "^Day ?[0-9]+$"))
  days <- as.integer(str_extract(hdr[dcol], "[0-9]+"))
  gcol <- 1; mcol <- 2
  r <- h + 1; keep <- integer(0)
  while (r <= nrow(raw)) {
    vals <- suppressWarnings(as.numeric(unlist(raw[r, dcol])))
    if (all(is.na(vals))) break
    keep <- c(keep, r); r <- r + 1
  }
  blk <- raw[keep, , drop = FALSE]
  tibble(Group = str_squish(as.character(blk[[gcol]])),
         Mouse = suppressWarnings(as.integer(blk[[mcol]]))) %>%
    bind_cols(as_tibble(lapply(setNames(dcol, days),
                               function(j) suppressWarnings(as.numeric(blk[[j]]))))) %>%
    mutate(Group = ifelse(Group == "" | is.na(Group), NA, Group)) %>%
    fill(Group, .direction = "down") %>%
    pivot_longer(-c(Group, Mouse), names_to = "Day", values_to = "value") %>%
    mutate(Day = as.integer(Day))
}

bw_raw <- parse_block(hdr_rows[1]) %>% rename(Weight = value)

# the % block (where present) is used only as a cross-check
bw_pct_file <- if (length(hdr_rows) > 1) parse_block(hdr_rows[2]) %>% rename(Pct_file = value) else NULL

dat <- bw_raw %>%
  filter(!is.na(Weight), !is.na(Mouse)) %>%
  group_by(Group, Mouse) %>%
  arrange(Day, .by_group = TRUE) %>%
  mutate(Baseline = Weight[Day == min(Day)],
         Pct      = 100 * Weight / Baseline,
         Loss     = 100 - Pct) %>%
  ungroup()

if (length(group_rename)) {
  dat$Group <- recode(dat$Group, !!!as.list(group_rename))
  group_levels <- recode(group_levels, !!!as.list(group_rename))
  ref_group    <- unname(ifelse(ref_group %in% names(group_rename), group_rename[[ref_group]], ref_group))
  names(group_cols) <- recode(names(group_cols), !!!as.list(group_rename))
}

dat <- dat %>%
  mutate(Group = factor(Group, levels = intersect(group_levels, unique(Group))),
         ID    = factor(paste0(Group, "_", Mouse)))
if (is.na(last_day)) last_day <- max(dat$Day)

# cross-check against the % block in the file
if (!is.null(bw_pct_file)) {
  chk <- dat %>% left_join(bw_pct_file %>% rename(Pct_file = Pct_file), by = c("Group", "Mouse", "Day"))
  d <- max(abs(chk$Pct - chk$Pct_file), na.rm = TRUE)
  message(sprintf("  Cross-check against the %% block in the file: max difference = %.4g %%", d))
}
message(sprintf("  %d mice | %d groups (%s) | Day %d-%d",
                n_distinct(dat$ID), nlevels(dat$Group),
                paste(levels(dat$Group), collapse = ", "), min(dat$Day), max(dat$Day)))

dir.create(out_fig, showWarnings = FALSE, recursive = TRUE); dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)


## ============================ 2. SHARED FUNCTIONS ============================
games_howell <- function(y, g) {                       # no equal-variance assumption
  g <- droplevels(factor(g)); k <- nlevels(g)
  m <- tapply(y, g, mean); v <- tapply(y, g, var); nn <- tapply(y, g, length)
  cmb <- combn(levels(g), 2)
  map_dfr(seq_len(ncol(cmb)), function(i) {
    a <- cmb[1, i]; b <- cmb[2, i]
    va <- v[[a]] / nn[[a]]; vb <- v[[b]] / nn[[b]]
    est <- m[[b]] - m[[a]]; se <- sqrt(va + vb)
    df <- (va + vb)^2 / (va^2 / (nn[[a]] - 1) + vb^2 / (nn[[b]] - 1))
    tibble(comparison = paste0(b, "-", a), estimate = est, df = df,
           p = ptukey(abs(est) / se * sqrt(2), nmeans = k, df = df, lower.tail = FALSE))
  })
}
perm_pairwise <- function(y, g) {                    
  g <- droplevels(factor(g)); cmb <- combn(levels(g), 2)
  map_dfr(seq_len(ncol(cmb)), function(i) {
    a <- cmb[1, i]; b <- cmb[2, i]
    ya <- y[g == a]; yb <- y[g == b]; v <- c(ya, yb); na <- length(ya)
    d0 <- mean(yb) - mean(ya)
    cb <- combn(na + length(yb), na)
    dp <- apply(cb, 2, function(ii) mean(v[-ii]) - mean(v[ii]))
    tibble(comparison = paste0(b, "-", a), estimate = d0,
           p_raw = mean(abs(dp) >= abs(d0) - 1e-12), n_perm = ncol(cb))
  }) %>% mutate(p = p.adjust(p_raw, "BH"))
}
safe_name <- function(x) x %>% str_replace_all("ï", "i") %>% str_replace_all("[^A-Za-z0-9._-]+", "_")
pick_pdf_device <- function() {
  probe <- function(dev) {
    f <- tempfile(fileext = ".pdf"); n0 <- dev.cur(); ok <- FALSE
    try({ suppressWarnings(dev(f, width = 2, height = 2))
          if (dev.cur() != n0) { plot.new(); dev.off(); ok <- file.exists(f) && file.size(f) > 0 } },
        silent = TRUE)
    while (dev.cur() != n0 && dev.cur() > 1) dev.off(); unlink(f); ok
  }
  if (probe(function(f, width, height) grDevices::cairo_pdf(f, width, height)))
    return(function(filename, width, height, ...) grDevices::cairo_pdf(filename, width, height))
  if (isTRUE(capabilities("aqua")))
    return(function(filename, width, height, ...)
      grDevices::quartz(file = filename, type = "pdf", width = width, height = height))
  function(filename, width, height, ...) grDevices::pdf(filename, width, height, useDingbats = FALSE)
}
dev_pdf <- pick_pdf_device()
save_plot <- function(p, file, w, h) {
  ggsave(file.path(out_fig, paste0(file, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  ok <- tryCatch({ ggsave(file.path(out_fig, paste0(file, ".pdf")), p, width = w, height = h,
                          device = dev_pdf); TRUE }, error = function(e) FALSE)
  message("  + ", out_fig, "/", file, ".png", if (ok) " / .pdf" else "")
}
write_csv_utf8 <- function(x, path) {
  fmt <- function(col) {
    if (is.numeric(col)) ifelse(is.na(col), "", format(col, trim = TRUE, digits = 6))
    else paste0('"', gsub('"', '""', ifelse(is.na(col), "", as.character(col)), fixed = TRUE), '"')
  }
  txt <- paste0(paste(sprintf('"%s"', names(x)), collapse = ","), "\n",
                paste(do.call(paste, c(lapply(x, fmt), sep = ",")), collapse = "\n"), "\n")
  con <- file(path, open = "wb"); on.exit(close(con))
  writeBin(charToRaw(enc2utf8(paste0("﻿", txt))), con)
}
put <- function(x, f) { write_csv_utf8(x, file.path(out_tab, f)); message("  + ", out_tab, "/", f) }

theme_bw2 <- function(base = 12) {
  theme_classic(base_size = base) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(), legend.position = "top",
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = rel(.9)),
          plot.caption = element_text(colour = "grey45", size = rel(.75), hjust = 0))
}


## ============================ 3. BODY WEIGHT: DESCRIPTION ====================
summ <- dat %>%
  group_by(Group, Day) %>%
  summarise(n = n(), mean_g = mean(Weight), sd_g = sd(Weight),
            mean_pct = mean(Pct), sd_pct = sd(Pct), .groups = "drop") %>%
  mutate(sem_pct = sd_pct / sqrt(n), sem_g = sd_g / sqrt(n),
         err_pct = if (err_type == "SEM") sem_pct else sd_pct,
         err_g   = if (err_type == "SEM") sem_g   else sd_g)
put(summ, "BW01_summary_by_group_day.csv")

pos <- position_dodge(width = 0.18)
cap_bw <- paste0("Mean +/- ", err_type, ", n = 4/group",
                 if (show_indiv) "; faint lines = individual mice" else "")

p_pct <- ggplot(summ, aes(Day, mean_pct, colour = Group, fill = Group)) +
  { if (show_indiv) geom_line(data = dat, aes(y = Pct, group = ID), alpha = .28, linewidth = .4) } +
  geom_hline(yintercept = 100, linetype = 2, colour = "grey65") +
  geom_errorbar(aes(ymin = mean_pct - err_pct, ymax = mean_pct + err_pct),
                width = .12, linewidth = .5, position = pos) +
  geom_line(linewidth = .9, position = pos) +
  geom_point(size = 2.4, shape = 21, colour = "black", stroke = .4, position = pos) +
  scale_colour_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  scale_x_continuous(breaks = sort(unique(dat$Day))) +
  labs(title = "Body weight change after infection",
       subtitle = "Percent of Day 0", x = "Day post infection",
       y = "Body weight (% of Day 0)", caption = cap_bw) +
  theme_bw2()
save_plot(p_pct, "BW_percent_of_day0", 6.6, 4.6)

p_g <- ggplot(summ, aes(Day, mean_g, colour = Group, fill = Group)) +
  { if (show_indiv) geom_line(data = dat, aes(y = Weight, group = ID), alpha = .28, linewidth = .4) } +
  geom_errorbar(aes(ymin = mean_g - err_g, ymax = mean_g + err_g), width = .12,
                linewidth = .5, position = pos) +
  geom_line(linewidth = .9, position = pos) +
  geom_point(size = 2.4, shape = 21, colour = "black", stroke = .4, position = pos) +
  scale_colour_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  scale_x_continuous(breaks = sort(unique(dat$Day))) +
  labs(title = "Body weight (gram)", x = "Day post infection", y = "Body weight (g)",
       caption = cap_bw) +
  theme_bw2()
save_plot(p_g, "BW_raw_grams", 6.6, 4.6)

p_ind <- ggplot(dat, aes(Day, Pct, colour = Group)) +
  geom_hline(yintercept = 100, linetype = 2, colour = "grey65") +
  geom_line(aes(group = ID), linewidth = .6) + geom_point(size = 1.5) +
  facet_wrap(~ Group) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_x_continuous(breaks = sort(unique(dat$Day))) +
  labs(title = "Individual body-weight trajectories", x = "Day post infection",
       y = "Body weight (% of Day 0)") +
  theme_bw2()
save_plot(p_ind, "BW_individual_trajectories", 8.4, 3.4)


## ============================ 4. BODY WEIGHT: STATISTICS =====================
per_mouse <- dat %>%
  group_by(Group, ID) %>%
  summarise(nadir_pct  = min(Pct),
            day_nadir  = Day[which.min(Pct)],
            pct_lastday = Pct[Day == max(Day)],
            AUC_loss   = { d <- Day[order(Day)]; l <- pmax(Loss[order(Day)], 0)
                           sum(diff(d) * (head(l, -1) + tail(l, -1)) / 2) },
            .groups = "drop")
put(per_mouse, "BW02_per_mouse_summary.csv")

mouse_tests <- map_dfr(c("nadir_pct", "pct_lastday", "AUC_loss"), function(v) {
  y <- per_mouse[[v]]; g <- droplevels(per_mouse$Group)
  bind_rows(
    tibble(metric = v, method = "Welch ANOVA", comparison = "omnibus",
           estimate = NA_real_, p = oneway.test(y ~ g, var.equal = FALSE)$p.value),
    tibble(metric = v, method = "Kruskal-Wallis", comparison = "omnibus",
           estimate = NA_real_, p = kruskal.test(y, g)$p.value),
    games_howell(y, g)  %>% transmute(metric = v, method = "Games-Howell", comparison, estimate, p),
    perm_pairwise(y, g) %>% transmute(metric = v, method = "Permutation exact (BH)", comparison, estimate, p)
  )
})
put(mouse_tests, "BW03_tests_nadir_lastday_AUC.csv")

perday <- map_dfr(sort(unique(dat$Day))[-1], function(dd) {
  d <- dat %>% filter(Day == dd)
  gh <- games_howell(d$Pct, droplevels(d$Group)) %>% mutate(method = "Games-Howell")
  pp <- perm_pairwise(d$Pct, droplevels(d$Group)) %>% transmute(comparison, estimate, p,
                                                                method = "Permutation exact (BH)")
  bind_rows(gh %>% select(comparison, estimate, p, method), pp) %>% mutate(Day = dd, .before = 1)
}) %>% group_by(method) %>% mutate(q_BH_allday = p.adjust(p, "BH")) %>% ungroup()
put(perday, "BW04_per-day_vs_groups.csv")

lmm_out <- tibble(); lmm_slope <- tibble()
if (have("nlme")) {
  suppressPackageStartupMessages(library(nlme))
  ctl <- nlme::lmeControl(opt = "optim", msMaxIter = 300, returnObject = TRUE)
  dm  <- dat %>% mutate(Dayf = factor(Day), Group = relevel(Group, ref = ref_group))

  # (i) Day as a factor -> effect at each timepoint
  m_fac <- try(nlme::lme(Pct ~ Group * Dayf, random = ~ 1 | ID,
                         correlation = nlme::corAR1(form = ~ Day | ID),
                         data = dm, control = ctl, method = "REML"), silent = TRUE)
  if (inherits(m_fac, "try-error"))
    m_fac <- nlme::lme(Pct ~ Group * Dayf, random = ~ 1 | ID, data = dm, control = ctl)
  an <- as.data.frame(anova(m_fac, type = "marginal")) %>% rownames_to_column("term") %>%
    filter(term != "(Intercept)")
  lmm_out <- tibble(model = "Pct ~ Group*Day(factor), 1|ID, AR1", term = an$term,
                    numDF = an$numDF, denDF = an$denDF, F = an$`F-value`, p = an$`p-value`)

  m_sl <- try(nlme::lme(Pct ~ Group * Day, random = ~ Day | ID, data = dm, control = ctl), silent = TRUE)
  if (inherits(m_sl, "try-error"))
    m_sl <- nlme::lme(Pct ~ Group * Day, random = ~ 1 | ID, data = dm, control = ctl)
  an2 <- as.data.frame(anova(m_sl, type = "marginal")) %>% rownames_to_column("term") %>%
    filter(term != "(Intercept)")
  lmm_out <- bind_rows(lmm_out,
    tibble(model = "Pct ~ Group*Day(numeric), Day|ID", term = an2$term, numDF = an2$numDF,
           denDF = an2$denDF, F = an2$`F-value`, p = an2$`p-value`))

  lmm_slope <- map_dfr(levels(dm$Group), function(rf) {
    d <- dm %>% mutate(Group = relevel(Group, ref = rf))
    m <- try(nlme::lme(Pct ~ Group * Day, random = ~ Day | ID, data = d, control = ctl), silent = TRUE)
    if (inherits(m, "try-error")) m <- nlme::lme(Pct ~ Group * Day, random = ~ 1 | ID, data = d, control = ctl)
    tt <- as.data.frame(summary(m)$tTable) %>% rownames_to_column("term")
    slope_ref <- tt %>% filter(term == "Day")
    tibble(reference = rf,
           term = c(paste0("slope[", rf, "]"), tt$term[str_detect(tt$term, "^Group.*:Day$")]),
           estimate = c(slope_ref$Value, tt$Value[str_detect(tt$term, "^Group.*:Day$")]),
           SE = c(slope_ref$Std.Error, tt$Std.Error[str_detect(tt$term, "^Group.*:Day$")]),
           p  = c(slope_ref$`p-value`, tt$`p-value`[str_detect(tt$term, "^Group.*:Day$")]))
  })
  put(lmm_out,   "BW05_LMM_anova.csv")
  put(lmm_slope, "BW06_LMM_slopes.csv")
}


sum_long <- per_mouse %>%
  transmute(Group, ID,
            `Nadir (% of Day 0)` = nadir_pct,
            `Day 5 (% of Day 0)` = pct_lastday,
            `AUC weight loss (%*day)` = AUC_loss) %>%
  pivot_longer(-c(Group, ID)) %>%
  mutate(name = factor(name, levels = c("Nadir (% of Day 0)", "Day 5 (% of Day 0)",
                                        "AUC weight loss (%*day)")))
href <- tibble(name = factor(c("Nadir (% of Day 0)", "Day 5 (% of Day 0)"),
                             levels = levels(sum_long$name)), yint = 100)

p_sum <- ggplot(sum_long, aes(Group, value, colour = Group)) +
  geom_hline(data = href, aes(yintercept = yint), linetype = 2, colour = "grey65") +
  stat_summary(fun = mean, geom = "crossbar", width = .5, linewidth = .4, fatten = 1.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .18, linewidth = .5) +
  geom_jitter(width = .11, height = 0, shape = 21, fill = "white", size = 2.1, stroke = .6) +
  facet_wrap(~ name, scales = "free_y") +
  scale_colour_manual(values = group_cols, guide = "none") +
  labs(title = "Per-mouse body-weight summaries", x = NULL, y = NULL,
       caption = "Bar = mean, whiskers = SEM, dots = individual mice. AUC = area under the weight-loss curve (Day 0-5).") +
  theme_bw2()
save_plot(p_sum, "BW_summary_metrics", 8.4, 3.8)


## ============================ 5. SURVIVAL ====================================
message("\nSurvival ...")
mice <- dat %>% distinct(Group, Mouse, ID)

surv_real <- mice %>%
  left_join(death_events %>% mutate(Group = as.character(Group)), by = c("Group", "Mouse")) %>%
  mutate(time = ifelse(is.na(time), last_day, time),
         status = ifelse(is.na(status), 0L, as.integer(status)),
         source = "observed events")

surv_wl <- function(thr) {
  dat %>% group_by(Group, ID, Mouse) %>%
    summarise(hit = any(Loss >= thr),
              time = if (any(Loss >= thr)) min(Day[Loss >= thr]) else last_day,
              status = as.integer(any(Loss >= thr)), .groups = "drop") %>%
    mutate(threshold = thr, source = paste0("weight loss >= ", thr, "%")) %>%
    select(Group, Mouse, ID, time, status, threshold, source)
}
surv_endpoint <- map_dfr(wl_sensitivity, surv_wl)

surv_all <- bind_rows(surv_real %>% mutate(threshold = NA_real_), surv_endpoint)
put(surv_all, "SURV01_survival_table.csv")

surv_summary <- surv_all %>% group_by(source, Group) %>%
  summarise(n = n(), events = sum(status), survival_pct = 100 * (1 - sum(status) / n()),
            .groups = "drop")
put(surv_summary, "SURV02_survival_summary.csv")
cat("\n===== SURVIVAL SUMMARY =====\n"); print(as.data.frame(surv_summary), row.names = FALSE)

km_curve <- function(sdf) {                       # KM computed directly (no survminer)
  sdf %>% group_by(Group) %>% group_modify(~ {
    d <- .x %>% arrange(time)
    n0 <- nrow(d); tt <- sort(unique(d$time))
    s <- 1; out <- tibble(time = 0, surv = 1, n_risk = n0, n_event = 0L, n_censor = 0L)
    at_risk <- n0
    for (t in tt) {
      ev <- sum(d$time == t & d$status == 1); cs <- sum(d$time == t & d$status == 0)
      if (at_risk > 0 && ev > 0) s <- s * (1 - ev / at_risk)
      out <- bind_rows(out, tibble(time = t, surv = s, n_risk = at_risk,
                                   n_event = as.integer(ev), n_censor = as.integer(cs)))
      at_risk <- at_risk - ev - cs
    }
    out
  }) %>% ungroup()
}

logrank <- function(sdf) {
  if (!have("survival") || sum(sdf$status) == 0)
    return(tibble(test = "log-rank", chisq = NA_real_, df = NA_integer_, p = NA_real_,
                  note = "no events -> not testable"))
  sd <- survival::survdiff(survival::Surv(time, status) ~ Group, data = sdf)
  tibble(test = "log-rank", chisq = unname(sd$chisq), df = length(sd$n) - 1L,
         p = pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE), note = "")
}

lr_tbl <- surv_all %>% group_by(source) %>% group_modify(~ logrank(.x)) %>% ungroup()
put(lr_tbl, "SURV03_logrank.csv")

add_offset <- function(km) {
  ends <- km %>% group_by(Group) %>% slice_max(time, n = 1, with_ties = FALSE) %>% pull(surv)
  overlap <- n_distinct(round(ends, 8)) < length(ends)
  km %>% mutate(y = 100 * surv +
                  if (overlap) (as.integer(Group) - (nlevels(Group) + 1) / 2) * km_offset else 0,
                .overlap = overlap)
}

plot_km <- function(sdf, ttl, sub) {
  km <- km_curve(sdf) %>% mutate(Group = factor(Group, levels = levels(dat$Group))) %>% add_offset()
  ev <- sum(sdf$status)
  cens <- km %>% filter(n_censor > 0)
  note <- paste0(if (ev == 0)
                   "No events -> 100% survival in all three groups; + = censored at the end of follow-up."
                 else "+ = censored.",
                 if (any(km$.overlap)) "\nOverlapping curves are nudged vertically for visibility (values unchanged)." else "")
  ggplot(km, aes(time, y, colour = Group, linetype = Group)) +
    geom_step(linewidth = 1, direction = "hv") +
    { if (nrow(cens)) geom_point(data = cens, shape = 3, size = 2.6, stroke = 1) } +
    scale_colour_manual(values = group_cols) +
    scale_linetype_manual(values = c("solid", "22", "42")[seq_len(nlevels(km$Group))]) +
    scale_x_continuous(breaks = sort(unique(dat$Day)), limits = c(0, last_day)) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
    labs(title = ttl, subtitle = sub, x = "Day post infection", y = "Survival (%)",
         caption = note) +
    theme_bw2()
}

km_plots <- surv_all %>% distinct(source) %>% pull(source)
for (src in km_plots) {
  sdf <- surv_all %>% filter(source == src)
  lr  <- lr_tbl %>% filter(source == src)
  sub <- if (is.na(lr$p)) "No events - log-rank not applicable"
         else sprintf("Log-rank: chi2 = %.2f, df = %d, p = %.3g", lr$chisq, lr$df, lr$p)
  p <- plot_km(sdf, paste0("Kaplan-Meier - ", src), sub)
  save_plot(p, paste0("SURV_KM_", safe_name(src)), 6.2, 4.4)
}

p_multi <- surv_endpoint %>% group_by(threshold) %>% group_modify(~ km_curve(.x)) %>%
  mutate(Group = factor(Group, levels = levels(dat$Group))) %>%
  group_modify(~ add_offset(.x)) %>% ungroup() %>%
  mutate(threshold = factor(paste0("Endpoint: weight loss >= ", threshold, "%"),
                            levels = paste0("Endpoint: weight loss >= ", wl_sensitivity, "%"))) %>%
  ggplot(aes(time, y, colour = Group, linetype = Group)) +
  geom_step(linewidth = 1, direction = "hv") +
  facet_wrap(~ threshold) +
  scale_colour_manual(values = group_cols) +
  scale_linetype_manual(values = c("solid", "22", "42")) +
  scale_x_continuous(breaks = sort(unique(dat$Day))) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity analysis over humane-endpoint thresholds",
       subtitle = "Proportion of mice NOT yet past the given weight-loss threshold",
       x = "Day post infection", y = "% mice before endpoint",
       caption = paste0("Maximum weight loss observed in the experiment: ",
                        sprintf("%.1f%%", max(dat$Loss)),
                        " -> no mouse reached the 20% threshold (standard humane endpoint for influenza models).",
                        "\nOverlapping curves are nudged vertically for visibility.")) +
  theme_bw2()
save_plot(p_multi, "SURV_endpoint_sensitivity", 9.2, 3.8)


## ============================ 6. CONSOLE SUMMARY ============================
cat("\n===== BODY WEIGHT: % of Day 0 (mean) =====\n")
summ %>% transmute(Group, Day, n, mean_pct = round(mean_pct, 2), sd_pct = round(sd_pct, 2)) %>%
  pivot_wider(id_cols = Group, names_from = Day, names_prefix = "D",
              values_from = mean_pct) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n===== Maximum weight loss per mouse =====\n")
per_mouse %>% mutate(max_loss_pct = round(100 - nadir_pct, 2)) %>%
  select(Group, ID, max_loss_pct, day_nadir, AUC_loss) %>%
  mutate(AUC_loss = round(AUC_loss, 2)) %>%
  arrange(desc(max_loss_pct)) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n===== Tests on the per-mouse summaries (p < 0.1) =====\n")
mouse_tests %>% filter(p < 0.1) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n===== Per-day comparisons (p < 0.1) =====\n")
perday %>% filter(p < 0.1) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (nrow(lmm_out)) {
  cat("\n===== LMM =====\n")
  lmm_out %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame() %>% print(row.names = FALSE)
  cat("\n----- Rate of weight change (%/day) -----\n")
  lmm_slope %>% filter(str_detect(term, "^slope|:Day$")) %>%
    mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame() %>% print(row.names = FALSE)
}

cat("\n===== LOG-RANK =====\n"); print(as.data.frame(lr_tbl), row.names = FALSE)
cat("\nDone.\n")
