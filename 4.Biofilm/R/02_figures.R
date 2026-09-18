################################################################################
# 02_figures.R
#
# Figures for biofilm production:
#   fig_biofilm_relative_dist.pdf   posterior density relative to the ancestor
#   fig_biofilm_absolute_dist.pdf   posterior density on the measured scale
#   fig_biofilm_mean_sd.pdf         per-strain geometric mean, wells behind
#
# All three carry nested group brackets: mutation (sinR, slrC), cell type
# (Spore, Vegetative), and origin (Ancestor, Evolved). Strain names are off the
# figures by default -- see show_clones.
#
# Everything is drawn on a log10 axis. Biofilm OD spans two orders of magnitude
# across these strains, the model is fit on the log scale, and one strain's
# mean minus its SD is negative, so a linear axis would be both unreadable and
# unable to carry an error bar.
#
# Input : output/posteriors.rds (R/01_bayes_biofilm.R), data/biofil.csv
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the biofilm project root (or R/).")
source(.setup[1])

posterior_file <- file.path(OUT_DIR, "posteriors.rds")
if (!file.exists(posterior_file)) {
  stop("output/posteriors.rds not found -- run R/01_bayes_biofilm.R first.")
}
post  <- readRDS(posterior_file)$biofilm
assay <- read.csv(file.path(DATA_DIR, "biofil.csv"), stringsAsFactors = FALSE)

# Tick labels for a log10 axis, written in the original OD units.
od_labels <- function(x) formatC(10^x, format = "g", digits = 2)

## ---- ridgeline plots -------------------------------------------------------

to_long <- function(draws, order) {
  as.data.frame(draws[, order, drop = FALSE]) |>
    pivot_longer(everything(), names_to = "clone", values_to = "value") |>
    mutate(clone = factor(clone, levels = order))
}

#' Ridgeline plot of posterior densities, one row per strain.
#'
#' Fill shades the two-tailed probability of each point, so the pale ends of a
#' ridge are its tails and the dark middle is its bulk. Group brackets sit
#' where the strain labels would otherwise be.
ridge_plot <- function(draws, order, xlab, xlim, scale = 2, ref = 0,
                       rows = NULL, labels_fn = waiver(), fig_width = 6.2,
                       show_clones = FALSE) {
  # Trim to the axis range up front. scale_*_continuous(limits=) would drop the
  # same rows, but it also clips the brackets, which sit outside the panel.
  d <- to_long(draws, order)
  d <- d[d$value >= xlim[1] & d$value <= xlim[2], ]

  span   <- diff(xlim)
  labels <- if (show_clones) short_label(order) else rep("", length(order))

  # Bracket geometry, in fractions of the x span, measured leftwards.
  b_start <- 0.05
  b_step  <- 0.105
  b_gap   <- 0.05

  p <- ggplot(d, aes(x = value, y = as.numeric(clone), group = clone,
                     fill = 0.5 - abs(0.5 - after_stat(ecdf))))

  if (!is.na(ref)) {
    p <- p + geom_vline(xintercept = ref, linetype = "dashed", colour = "grey25")
  }

  p <- p +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                        rel_min_height = 0.005, alpha = 0.5, scale = scale,
                        colour = "grey25") +
    scale_fill_gradient(low = "white", high = "grey50", name = "Tail prob.") +
    scale_x_continuous(breaks = scales::pretty_breaks(5), labels = labels_fn,
                       sec.axis = dup_axis(labels = labels_fn)) +
    scale_y_continuous(breaks = seq_along(labels), labels = labels,
                       sec.axis = sec_axis(~ ., breaks = seq_along(labels),
                                           labels = NULL)) +
    labs(x = xlab) +
    mytheme +
    theme(axis.title.y = element_blank(),
          # Strains are identified by the brackets, so the ticks that used to
          # point at each row have nothing left to mark.
          axis.ticks.y = element_blank(),
          axis.ticks.y.right = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.9, 0.82))

  left_pt <- 5
  if (!is.null(rows)) {
    p <- p + nested_brackets(rows, axis = "y",
                             start = xlim[1] - b_start * span,
                             step  = -b_step * span,
                             gap   = -b_gap * span,
                             size  = 4)
    frac    <- b_start + b_step * (length(rows) - 1) + b_gap
    left_pt <- bracket_margin(frac, fig_width, overhead_in = 0.45) + 6
  }

  p +
    coord_cartesian(xlim = xlim, clip = "off") +
    theme(plot.margin = margin(5, 10, 5, left_pt))
}

## ---- geometric mean +/- SD strip plot --------------------------------------

#' Per-strain geometric mean with a one-geometric-SD error bar over the eight
#' wells, individual wells in grey behind, strains left to right in
#' STRIP_ORDER, nested group brackets below.
#'
#' Mean and SD are taken on the log10 scale and the axis is labelled in OD
#' units, which is what the log-normal model assumes and what keeps the lower
#' error bar positive.
mean_sd_plot <- function(values, clones, ylab, fig_height = 5.4,
                         show_clones = FALSE) {
  d <- data.frame(clone = factor(clones, levels = STRIP_ORDER),
                  value = log10(values))
  d <- d[!is.na(d$clone), ]
  d$x <- as.numeric(d$clone)

  s <- d |>
    group_by(clone, x) |>
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")

  meta <- strain_meta(STRIP_ORDER)
  rows <- list(meta$mutation, meta$cell, meta$origin)

  lo  <- min(d$value, s$mean - s$sd)
  hi  <- max(d$value, s$mean + s$sd)
  pad <- hi - lo

  # Bracket geometry, in fractions of the y span, measured downwards.
  b_start <- if (show_clones) 0.22 else 0.10
  b_step  <- 0.11
  b_gap   <- 0.035
  frac    <- b_start + b_step * (length(rows) - 1) + b_gap

  p <- ggplot(d, aes(x = x, y = value)) +
    geom_jitter(width = 0.13, height = 0, colour = "grey72", size = 2.6,
                alpha = 0.75) +
    geom_errorbar(data = s, aes(x = x, y = mean, ymin = mean - sd,
                                ymax = mean + sd),
                  width = 0.16, linewidth = 0.6, colour = "black") +
    geom_point(data = s, aes(x = x, y = mean), size = 3.6, shape = 21,
               colour = "black", fill = "white", stroke = 0.9)

  if (show_clones) {
    p <- p + annotate("text", x = seq_along(STRIP_ORDER), y = lo - 0.10 * pad,
                      label = short_label(STRIP_ORDER), colour = "grey55",
                      size = 3.6, vjust = 1)
  }

  p +
    nested_brackets(rows, axis = "x", start = lo - b_start * pad,
                    step = -b_step * pad, gap = -b_gap * pad, size = 4) +
    scale_x_continuous(breaks = seq_along(STRIP_ORDER), sec.axis = dup_axis()) +
    scale_y_continuous(labels = od_labels,
                       sec.axis = dup_axis(labels = od_labels)) +
    labs(y = ylab) +
    coord_cartesian(ylim = c(lo - 0.02 * pad, hi + 0.02 * pad), clip = "off") +
    mytheme +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.x.top = element_blank(),
          plot.margin  = margin(
            5, 5,
            bracket_margin(frac, fig_height, overhead_in = 0.55) + 10,
            5))
}

## ---- draw ------------------------------------------------------------------

# Bracket rows, innermost first. The relative panel has no ancestor row -- the
# ancestor is the reference, so every strain shown is evolved.
meta_rel <- strain_meta(PLOT_ORDER)
rows_rel <- list(meta_rel$mutation, meta_rel$cell)

ABS_ORDER <- c("ancestor", PLOT_ORDER)
meta_abs  <- strain_meta(ABS_ORDER)
rows_abs  <- list(meta_abs$mutation, meta_abs$cell, meta_abs$origin)

rel <- ridge_plot(
  log10(post$relative), PLOT_ORDER,
  xlab = expression(paste("Relative biofilm, ", italic(log[10]))),
  xlim = c(-1.5, 1.5), scale = 2, ref = 0, rows = rows_rel)
ggsave(file.path(OUT_DIR, "fig_biofilm_relative_dist.pdf"), rel,
       width = 6.2, height = 5)

# The absolute panel gets its axis range from the credible intervals rather
# than the extreme quantiles: the vague Gamma prior on the precision leaves
# long tails on strains whose wells disagree.
abs_draws <- log10(post$absolute[, ABS_ORDER])
abs_xlim  <- range(log10(post$quantiles$absolute[ABS_ORDER, c("2.5%", "97.5%")]))
abs_xlim  <- abs_xlim + c(-1, 1) * 0.12 * diff(abs_xlim)

ab <- ridge_plot(abs_draws, ABS_ORDER, xlab = expression(Biofilm ~ (OD[550])),
                 xlim = abs_xlim, scale = 2, ref = NA, rows = rows_abs,
                 labels_fn = od_labels)
ggsave(file.path(OUT_DIR, "fig_biofilm_absolute_dist.pdf"), ab,
       width = 6.2, height = 5.4)

sd_plot <- mean_sd_plot(assay$OD550_Corrected, assay$clones,
                        ylab = expression(Biofilm ~ (OD[550])))
ggsave(file.path(OUT_DIR, "fig_biofilm_mean_sd.pdf"), sd_plot,
       width = 6.5, height = 5.4)

message("Wrote output/fig_biofilm_{relative_dist,absolute_dist,mean_sd}.pdf")
