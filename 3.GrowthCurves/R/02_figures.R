################################################################################
# 02_figures.R
#
# Figures for the growth parameters, three per parameter:
#   fig_<p>_relative_dist.pdf   posterior density relative to the ancestor
#   fig_<p>_absolute_dist.pdf   posterior density on the measured scale
#   fig_<p>_mean_sd.pdf         per-strain mean +/- SD over the six curve fits
#
# All three carry the manuscript's nested group brackets: mutation (sinR,
# ywcC), cell type (Spore, Total), and origin (Ancestor, Evolved).
#
# Input : output/posteriors.rds        (R/01_bayes_fitness.R)
#         data/comp_data_annotated.csv (R/91_build_comp_data.R)
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

posterior_file <- file.path(OUT_DIR, "posteriors.rds")
if (!file.exists(posterior_file)) {
  stop("output/posteriors.rds not found -- run R/01_bayes_fitness.R first.")
}
post <- readRDS(posterior_file)
fits <- read.csv(file.path(DATA_DIR, "comp_data_annotated.csv"),
                 stringsAsFactors = FALSE)

## ---- ridgeline plots -------------------------------------------------------

# Posterior draws (columns = strains) -> long data frame in plotting order.
to_long <- function(draws, order) {
  as.data.frame(draws[, order, drop = FALSE]) |>
    pivot_longer(everything(), names_to = "clone", values_to = "value") |>
    mutate(clone = factor(clone, levels = order))
}

#' Ridgeline plot of posterior densities, one row per strain.
#'
#' Fill shades the two-tailed probability of each point, so the pale ends of a
#' ridge are its tails and the dark middle is its bulk. Group brackets go on
#' the right, where the secondary axis has ticks but no labels.
#'
#' @param draws matrix of posterior draws, one column per strain
#' @param order strain order, bottom to top
#' @param rows  bracket rows (innermost first), or NULL for none
#' @param ref   x position of the dashed reference line, NA to omit
ridge_plot <- function(draws, order, xlab, xlim, scale = 2, ref = 1,
                       rows = NULL, fig_width = 6.2, show_clones = FALSE) {
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
    scale_x_continuous(breaks = scales::pretty_breaks(5),
                       sec.axis = dup_axis()) +
    scale_y_continuous(breaks = seq_along(labels), labels = labels,
                       sec.axis = sec_axis(~ ., breaks = seq_along(labels),
                                           labels = NULL)) +
    labs(x = xlab) +
    mytheme +
    theme(axis.title.y = element_blank(),
          # Strains are identified by the brackets now, so the ticks that used
          # to point at each row have nothing left to mark.
          axis.ticks.y = element_blank(),
          axis.ticks.y.right = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.9, 0.82))

  left_pt <- 5
  if (!is.null(rows)) {
    # Brackets sit where the strain labels used to: group membership is what
    # the panel is about, individual clone identity is in the data files.
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

## ---- mean +/- SD strip plot ------------------------------------------------

#' Per-strain mean with an SD error bar over the individual curve fits, in the
#' style of the manuscript's strain panel: raw curves in grey behind, strains
#' left to right in STRIP_ORDER, nested group brackets below.
#'
#' Strain labels are drawn as annotations rather than axis text so they stack
#' cleanly with the brackets underneath.
mean_sd_plot <- function(values, clones, ylab, fig_height = 5.4,
                         show_clones = FALSE) {
  d <- data.frame(clone = factor(clones, levels = STRIP_ORDER), value = values)
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

  # Bracket geometry, in fractions of the y span, measured downwards. The
  # brackets start lower when the strain labels are shown, to clear them.
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
    scale_x_continuous(breaks = seq_along(STRIP_ORDER),
                       sec.axis = dup_axis()) +
    scale_y_continuous(sec.axis = dup_axis()) +
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

## ---- what to draw ----------------------------------------------------------

PANELS <- list(
  umax = list(stem = "umax", column = "umax",
              relative_lab = expression("Relative " * mu * "max"),
              absolute_lab = expression(mu * "max"),
              strip_lab    = "Maximum growth rate",
              relative_xlim = c(0, 2),   relative_scale = 2),
  A    = list(stem = "yield", column = "A",
              relative_lab = "Relative yield",
              absolute_lab = "Maximum yield",
              strip_lab    = "Maximum yield",
              relative_xlim = c(0.2, 2), relative_scale = 2),
  L    = list(stem = "lag", column = "L",
              relative_lab = "Relative lag time",
              absolute_lab = "Lag time",
              strip_lab    = "Lag time",
              relative_xlim = c(0, 2),   relative_scale = 0.6)
)

# Bracket rows, innermost first. The relative panels have no ancestor row --
# the ancestor is the reference, so every strain shown is evolved.
meta_rel <- strain_meta(PLOT_ORDER)
rows_rel <- list(meta_rel$mutation, meta_rel$cell)

ABS_ORDER <- c("ancestor", PLOT_ORDER)
meta_abs  <- strain_meta(ABS_ORDER)
rows_abs  <- list(meta_abs$mutation, meta_abs$cell, meta_abs$origin)

for (nm in names(PANELS)) {
  cfg <- PANELS[[nm]]

  rel <- ridge_plot(post[[nm]]$relative, PLOT_ORDER, cfg$relative_lab,
                    xlim = cfg$relative_xlim, scale = cfg$relative_scale,
                    ref = 1, rows = rows_rel)
  ggsave(file.path(OUT_DIR, sprintf("fig_%s_relative_dist.pdf", cfg$stem)),
         rel, width = 6.2, height = 5)

  # Absolute panels get their axis range from the credible intervals rather
  # than the extreme quantiles: the vague Gamma prior on the precision leaves
  # long tails on strains whose six replicates disagree, and those tails would
  # otherwise stretch the axis far past anything worth reading.
  abs_draws <- post[[nm]]$absolute[, ABS_ORDER]
  abs_xlim  <- range(post[[nm]]$quantiles$absolute[ABS_ORDER, c("2.5%", "97.5%")])
  abs_xlim  <- abs_xlim + c(-1, 1) * 0.12 * diff(abs_xlim)

  ab <- ridge_plot(abs_draws, ABS_ORDER, cfg$absolute_lab,
                   xlim = abs_xlim, scale = 2, ref = NA, rows = rows_abs)
  ggsave(file.path(OUT_DIR, sprintf("fig_%s_absolute_dist.pdf", cfg$stem)),
         ab, width = 6.2, height = 5.4)

  sd_plot <- mean_sd_plot(fits[[cfg$column]], fits$clones, cfg$strip_lab)
  ggsave(file.path(OUT_DIR, sprintf("fig_%s_mean_sd.pdf", cfg$stem)),
         sd_plot, width = 6.5, height = 5.4)

  message("Wrote output/fig_", cfg$stem, "_{relative_dist,absolute_dist,mean_sd}.pdf")
}
