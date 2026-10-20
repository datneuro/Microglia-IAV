#!/usr/bin/env Rscript
# =============================================================================
# 02_Body-Weight-AUC.R — Area under the weight-loss curve
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(purrr)
})
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)))

# ------------------------------ SETTINGS -------------------------------------
xlsx_name    <- "20260826 Virus infection 5 days for HMC3.xlsx"
out_fig      <- file.path("Body weight", "figures")   # Body weight/ = data folder (repository root)
out_name     <- "BW_AUC_weight_loss"

group_levels <- c("Naïve", "PR8", "S88")
group_cols   <- c("Naïve" = "#9E9E9E", "PR8" = "#4C72B0", "S88" = "#C44E52")
group_rename <- c("PR8" = "A/PR/8/34",      
                  "S33" = "A/NWS/33",
                  "S88" = "A/NWS/33")     

y_lab        <- "AUC weight loss (% × day)"
plot_type    <- "bar"           
err_type     <- "SEM"           
clip_gain    <- TRUE            
annotate_p   <- FALSE           
p_method     <- "Games-Howell"  
show_ns      <- FALSE

fig_w <- 4.0; fig_h <- 3.8; dpi <- 400
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
  mutate(Pct = 100 * Weight / Weight[Day == min(Day)],
         Loss = 100 - Pct) %>% ungroup() %>%
  mutate(Group = factor(Group, levels = intersect(group_levels, unique(Group))),
         ID = factor(paste0(Group, "_", Mouse)))


## ---- AUC (trapezoidal rule) -------------------------------------------------
auc <- dat %>%
  group_by(Group, ID) %>%
  summarise(AUC = { d <- Day[order(Day)]
                    l <- Loss[order(Day)]
                    if (clip_gain) l <- pmax(l, 0)
                    sum(diff(d) * (head(l, -1) + tail(l, -1)) / 2) },
            .groups = "drop")

summ <- auc %>% group_by(Group) %>%
  summarise(n = n(), mean = mean(AUC), sd = sd(AUC), .groups = "drop") %>%
  mutate(sem = sd / sqrt(n), err = if (err_type == "SEM") sem else sd)

cat("== AUC weight loss by group ==\n")
summ %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame() %>% print(row.names = FALSE)


## ---- statistics -------------------------------------------------------------
gh <- function(y, g) {
  g <- droplevels(factor(g)); k <- nlevels(g)
  m <- tapply(y, g, mean); v <- tapply(y, g, var); nn <- tapply(y, g, length)
  cmb <- combn(levels(g), 2)
  map_dfr(seq_len(ncol(cmb)), function(i) {
    a <- cmb[1, i]; b <- cmb[2, i]
    va <- v[[a]] / nn[[a]]; vb <- v[[b]] / nn[[b]]
    est <- m[[b]] - m[[a]]; se <- sqrt(va + vb)
    df <- (va + vb)^2 / (va^2 / (nn[[a]] - 1) + vb^2 / (nn[[b]] - 1))
    tibble(comparison = paste0(b, "-", a), g1 = a, g2 = b, estimate = est,
           p = ptukey(abs(est) / se * sqrt(2), nmeans = k, df = df, lower.tail = FALSE))
  })
}
pm <- function(y, g) {
  g <- droplevels(factor(g)); cmb <- combn(levels(g), 2)
  map_dfr(seq_len(ncol(cmb)), function(i) {
    a <- cmb[1, i]; b <- cmb[2, i]
    ya <- y[g == a]; yb <- y[g == b]; v <- c(ya, yb); na <- length(ya)
    d0 <- mean(yb) - mean(ya); cb <- combn(na + length(yb), na)
    dp <- apply(cb, 2, function(ii) mean(v[-ii]) - mean(v[ii]))
    tibble(comparison = paste0(b, "-", a), g1 = a, g2 = b, estimate = d0,
           p_raw = mean(abs(dp) >= abs(d0) - 1e-12))
  }) %>% mutate(p = p.adjust(p_raw, "BH"))
}

ph <- if (p_method == "Permutation") pm(auc$AUC, auc$Group) else gh(auc$AUC, auc$Group)
cat("\n== Omnibus ==\n")
cat(sprintf("  Welch ANOVA    p = %.3f\n", oneway.test(AUC ~ Group, data = auc, var.equal = FALSE)$p.value))
cat(sprintf("  Kruskal-Wallis p = %.3f\n", kruskal.test(AUC ~ Group, data = auc)$p.value))
cat(sprintf("\n== Post-hoc (%s) ==\n", p_method))
ph %>% transmute(comparison, estimate = round(estimate, 2), p = signif(p, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                     ifelse(p < 0.05, "*", if (show_ns) "ns" else NA_character_)))


## ---- plot -------------------------------------------------------------------
top <- max(c(auc$AUC, summ$mean + summ$err))

p <- ggplot(summ, aes(Group, mean))

if (plot_type == "bar") {
  p <- p +
    geom_col(aes(fill = Group), width = .66, colour = "black", linewidth = .3) +
    geom_errorbar(aes(ymin = pmax(mean - err, 0), ymax = mean + err),
                  width = .2, linewidth = .45, colour = "black") +
    scale_fill_manual(values = group_cols)
} else {
  p <- p +
    geom_crossbar(aes(ymin = mean, ymax = mean, colour = Group),
                  width = .5, linewidth = .5, fatten = 0) +
    geom_errorbar(aes(ymin = mean - err, ymax = mean + err, colour = Group),
                  width = .18, linewidth = .5) +
    scale_colour_manual(values = group_cols)
}

p <- p +
  geom_jitter(data = auc, aes(Group, AUC), inherit.aes = FALSE,
              width = .12, height = 0, shape = 21, fill = "white", colour = "black",
              size = 2.1, stroke = .5) +
  scale_y_continuous(expand = expansion(mult = c(0, if (annotate_p) 0.22 else 0.08))) +
  labs(x = NULL, y = y_lab) +          # y title only; group names on the x axis
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text  = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        plot.title = element_blank(), plot.subtitle = element_blank(),
        plot.caption = element_blank())

if (annotate_p) {
  sig <- ph %>% mutate(lab = stars(p)) %>% filter(!is.na(lab)) %>%
    mutate(x1 = match(g1, levels(auc$Group)), x2 = match(g2, levels(auc$Group)),
           y = top * (1 + 0.07 * seq_len(n())))
  if (nrow(sig)) {
    p <- p +
      geom_segment(data = sig, aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = .4) +
      geom_text(data = sig, aes(x = (x1 + x2) / 2, y = y, label = lab),
                inherit.aes = FALSE, vjust = -0.2, size = 4.4)
  } else message("  (no pair reaches p < 0.05 -> no stars drawn)")
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
message("\n  + figures/", out_name, ".png", if (ok) " / .pdf" else "")
