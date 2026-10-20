#!/usr/bin/env Rscript
# =============================================================================
# 01_Body-Weight-Curve.R — Body-weight curve
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(stringr)
})
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)))

# ------------------------------ SETTINGS -------------------------------------
xlsx_name    <- "20260826 Virus infection 5 days for HMC3.xlsx"
out_fig      <- file.path("Body weight", "figures")   # Body weight/ = data folder (repository root)
out_name     <- "BW_curve_clean"

group_levels <- c("Naïve", "PR8", "S88")
group_cols   <- c("Naïve" = "#9E9E9E", "PR8" = "#4C72B0", "S88" = "#C44E52")
group_rename <- c("PR8" = "A/PR/8/34",      # display label of each group
                  "S33" = "A/NWS/33",
                  "S88" = "A/NWS/33")      # "S88" = group code of A/NWS/33 in the source sheet

y_lab        <- "Body weight (% of Day 0)"
x_lab        <- "Day post infection"
y_limits     <- NULL         
bold_titles  <- TRUE        

err_type     <- "SEM"        
show_indiv   <- TRUE         
show_hline   <- TRUE        
label_mode   <- "direct"     
label_size   <- 4.2
label_nudge  <- 0.12         
dodge_w      <- 0.18

fig_w <- 6.4; fig_h <- 4.4; dpi <- 400
# -----------------------------------------------------------------------------


## ---- read data --------------------------------------------------------------
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

raw <- as.data.frame(read_excel_safe(find_file(xlsx_name), sheet = 1,
                                     col_names = FALSE, .name_repair = "minimal"))
h <- which(apply(raw, 1, function(r)
  any(str_detect(str_squish(as.character(r)), "^Day ?0$"), na.rm = TRUE)))[1]
hdr  <- str_squish(as.character(unlist(raw[h, ])))
dcol <- which(str_detect(hdr, "^Day ?[0-9]+$"))
days <- as.integer(str_extract(hdr[dcol], "[0-9]+"))

r <- h + 1; keep <- integer(0)
while (r <= nrow(raw) && !all(is.na(suppressWarnings(as.numeric(unlist(raw[r, dcol])))))) {
  keep <- c(keep, r); r <- r + 1
}
blk <- raw[keep, , drop = FALSE]

dat <- tibble(Group = str_squish(as.character(blk[[1]])),
              Mouse = suppressWarnings(as.integer(blk[[2]]))) %>%
  bind_cols(as_tibble(lapply(setNames(dcol, days),
                             function(j) suppressWarnings(as.numeric(blk[[j]]))))) %>%
  mutate(Group = ifelse(Group == "" | is.na(Group), NA, Group)) %>%
  fill(Group, .direction = "down") %>%
  pivot_longer(-c(Group, Mouse), names_to = "Day", values_to = "Weight") %>%
  mutate(Day = as.integer(Day)) %>%
  filter(!is.na(Weight), !is.na(Mouse))

if (length(group_rename)) {
  dat$Group    <- recode(dat$Group, !!!as.list(group_rename))
  group_levels <- recode(group_levels, !!!as.list(group_rename))
  names(group_cols) <- recode(names(group_cols), !!!as.list(group_rename))
}

dat <- dat %>%
  group_by(Group, Mouse) %>% arrange(Day, .by_group = TRUE) %>%
  mutate(Pct = 100 * Weight / Weight[Day == min(Day)]) %>% ungroup() %>%
  mutate(Group = factor(Group, levels = intersect(group_levels, unique(Group))),
         ID = factor(paste0(Group, "_", Mouse)))

summ <- dat %>% group_by(Group, Day) %>%
  summarise(n = n(), mean_pct = mean(Pct), sd = sd(Pct), .groups = "drop") %>%
  mutate(err = if (err_type == "SEM") sd / sqrt(n) else sd)

lab_df <- summ %>% filter(Day == max(Day)) %>% mutate(x = Day + label_nudge, y_lab_pos = mean_pct)

if (!is.null(y_limits)) {
  pad <- diff(y_limits) * 0.03
  lab_df <- lab_df %>% mutate(y_lab_pos = pmin(pmax(mean_pct, y_limits[1] + pad), y_limits[2] - pad),
                              clamped   = abs(y_lab_pos - mean_pct) > 1e-9)
  hidden <- dat %>% filter(Pct > y_limits[2] | Pct < y_limits[1])
  if (nrow(hidden))
    message(sprintf("  WARNING: y_limits = [%g, %g] hides %d/%d points (observed range %.2f - %.2f).",
                    y_limits[1], y_limits[2], nrow(hidden), nrow(dat), min(dat$Pct), max(dat$Pct)))
  if (any(lab_df$clamped))
    message("  WARNING: labels of group(s) ", paste(lab_df$Group[lab_df$clamped], collapse = ", "),
            " fell outside the axis and were moved inside.")
  gap <- diff(y_limits) * 0.045
  lab_df <- lab_df %>% arrange(desc(y_lab_pos))
  if (nrow(lab_df) > 1)
    for (i in 2:nrow(lab_df))
      lab_df$y_lab_pos[i] <- min(lab_df$y_lab_pos[i], lab_df$y_lab_pos[i - 1] - gap)
}


## ---- plot -------------------------------------------------------------------
pos <- position_dodge(width = dodge_w)

p <- ggplot(summ, aes(Day, mean_pct, colour = Group, fill = Group)) +
  { if (show_hline) geom_hline(yintercept = 100, linetype = 2, colour = "grey65") } +
  { if (show_indiv) geom_line(data = dat, aes(y = Pct, group = ID),
                              alpha = .28, linewidth = .4) } +
  geom_errorbar(aes(ymin = mean_pct - err, ymax = mean_pct + err),
                width = .12, linewidth = .5, position = pos) +
  geom_line(linewidth = .9, position = pos) +
  geom_point(size = 2.4, shape = 21, colour = "black", stroke = .4, position = pos) +
  scale_colour_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  scale_x_continuous(breaks = sort(unique(dat$Day)),
                     expand = expansion(mult = c(0.03, if (label_mode == "direct") 0.21 else 0.05))) +
  labs(x = x_lab, y = y_lab) +           # axis titles only
  theme_classic(base_size = 12) +
  theme(axis.text  = element_text(colour = "black"),
        axis.title = element_text(colour = "black",
                                  face = if (bold_titles) "bold" else "plain"),
        plot.title = element_blank(), plot.subtitle = element_blank(),
        plot.caption = element_blank())

if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)

if (label_mode == "direct") {
  p <- p + geom_text(data = lab_df, aes(x = x, y = y_lab_pos, label = Group, colour = Group),
                     inherit.aes = FALSE, hjust = 0, vjust = 0.5,
                     size = label_size, fontface = "bold", show.legend = FALSE) +
    theme(legend.position = "none")
} else {
  p <- p + theme(legend.position = "top", legend.title = element_blank())
}


## ---- save -------------------------------------------------------------------
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_fig, paste0(out_name, ".png")), p,
       width = fig_w, height = fig_h, dpi = dpi, bg = "white")

pdf_dev <- local({
  probe <- function(dev) {
    f <- tempfile(fileext = ".pdf"); n0 <- dev.cur(); ok <- FALSE
    try({ suppressWarnings(dev(f, width = 2, height = 2))
          if (dev.cur() != n0) { plot.new(); dev.off(); ok <- file.exists(f) && file.size(f) > 0 } },
        silent = TRUE)
    while (dev.cur() != n0 && dev.cur() > 1) dev.off(); unlink(f); ok
  }
  if (probe(function(f, width, height) grDevices::cairo_pdf(f, width, height)))
    function(filename, width, height, ...) grDevices::cairo_pdf(filename, width, height)
  else if (isTRUE(capabilities("aqua")))
    function(filename, width, height, ...)
      grDevices::quartz(file = filename, type = "pdf", width = width, height = height)
  else function(filename, width, height, ...) grDevices::pdf(filename, width, height)
})
ok <- tryCatch({ ggsave(file.path(out_fig, paste0(out_name, ".pdf")), p,
                        width = fig_w, height = fig_h, device = pdf_dev); TRUE },
               error = function(e) FALSE)

message("  + figures/", out_name, ".png", if (ok) " / .pdf" else "")
message("Done.")
