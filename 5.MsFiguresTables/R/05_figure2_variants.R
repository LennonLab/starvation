################################################################################
# 05_figure2_variants.R
#
# Three versions of Figure 2, to choose between.
#
# The figure is narrowed to the two parameters that carry the result: maximum
# growth rate and biofilm. Yield and lag move to the supplement -- they are
# reported, but nothing in the argument turns on them.
#
#   v1  posterior ridges        what the model actually estimates, densities and
#                               all. Truthful, and unfamiliar to many readers.
#   v2  mean +/- SD by clone    the traditional plot: every replicate shown,
#                               mean and spread. Familiar, but it throws away
#                               the inverse-variance weighting the estimates
#                               are built on and shows only replicate scatter.
#   v3  posterior median +/-    the middle: the model's estimate and its 95%
#       95% credible interval   credible interval, drawn as a point and a bar.
#                               Reads like a conventional figure; the number
#                               plotted is still the weighted one.
#
# Output: output/figure2_v1_posterior_ridges.pdf / .png
#         output/figure2_v2_mean_sd.pdf / .png
#         output/figure2_v3_posterior_interval.pdf / .png
#         output/figureS_growth_yield_lag.pdf / .png
################################################################################

.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the manuscript-figures project root (or R/).")
source(.setup[1])

g <- load_project("growth",  "02_figures.R")
b <- load_project("biofilm", "02_figures.R")
require_objects(g, "ridge_plot", "mean_sd_plot", "post", "ABS_ORDER", "rows_abs",
                "STRIP_ORDER", "PANELS", "fits", what = "3.GrowthCurves")
require_objects(b, "ridge_plot", "mean_sd_plot", "post", "ABS_ORDER", "rows_abs",
                "STRIP_ORDER", "od_labels", "assay", what = "4.Biofilm")

abs_xlim_for <- function(quant, order, pad = 0.12) {
  r <- range(quant[order, c("2.5%", "97.5%")])
  r + c(-1, 1) * pad * diff(r)
}
ancestor_ref <- function(draws) mean(draws[, "ancestor"])

## ---- v1: posterior ridges ---------------------------------------------------

# Manuscript form of the ridges, made to say exactly what v3 says:
#   * scale 0.9, so each density stays inside its own row. At scale 2 a ridge
#     rose two rows above its baseline and the brackets, which are placed at the
#     baselines, looked shifted onto the wrong clones.
#   * a reference line in A only. B has none, for the same reason v3 has none:
#     the biofilm reference well is B. subtilis 168 delta 6, not the ancestor.
#   * no tail-probability legend; the shading is explained in the caption.
#   * brackets on A; B shares A's row order (both ancestor, M4, M13, M79,
#     M17-M54, M23, M26), so they read straight across.
RIDGE_SCALE <- 0.9

v1_umax <- g$ridge_plot(
  g$post$umax$absolute[, g$ABS_ORDER], g$ABS_ORDER, g$PANELS$umax$absolute_lab,
  xlim = abs_xlim_for(g$post$umax$quantiles$absolute, g$ABS_ORDER),
  scale = RIDGE_SCALE, ref = ancestor_ref(g$post$umax$absolute), rows = g$rows_abs) +
  theme(legend.position = "none")

v1_bio <- b$ridge_plot(
  log10(b$post$absolute[, b$ABS_ORDER]), b$ABS_ORDER,
  xlab = expression(Biofilm ~ (OD[550])),
  xlim = abs_xlim_for(log10(b$post$quantiles$absolute), b$ABS_ORDER),
  scale = RIDGE_SCALE, ref = NA,
  rows = NULL, labels_fn = b$od_labels) +
  theme(legend.position = "none")

# The biofilm reference well (B. subtilis 168 delta 6, not the ancestor) stays
# in panel B's bottom row for now. Whether to keep it or drop it is with the
# coauthors; the caption says what it is.
## ---- v2: mean +/- SD over replicates, clones on x ---------------------------

v2_umax <- g$mean_sd_plot(g$fits$umax, g$fits$clones, g$PANELS$umax$strip_lab)
v2_bio  <- b$mean_sd_plot(b$assay$OD550_Corrected, b$assay$clones,
                          expression(Biofilm ~ (OD[550])))

## ---- v3: posterior median and 95% credible interval, clones on x ------------

#' The estimate the model produces, drawn the way a reader expects to see one.
#'
#' Same quantity as the ridges -- the inverse-variance weighted posterior --
#' reduced to a median and an interval. Raw replicates go behind it in grey so
#' the reader can still see the data the estimate came from.
interval_plot <- function(env, quant, raw_value, raw_clone, ylab,
                          log_y = FALSE, fig_height = 5.4,
                          ancestor_line = TRUE) {
  ord <- env$STRIP_ORDER

  # Biofilm spans two orders of magnitude, so it is drawn on log10 values with
  # a linear axis whose labels are back-transformed -- the same arrangement
  # mean_sd_plot() uses. A real log scale would break the brackets, whose
  # offsets below the panel are additive.
  tr <- if (log_y) log10 else identity
  q  <- tr(quant[ord, , drop = FALSE])

  est <- data.frame(x = seq_along(ord),
                    lo = q[, "2.5%"], mid = q[, "50%"], hi = q[, "97.5%"])

  raw <- data.frame(clone = factor(raw_clone, levels = ord),
                    value = tr(raw_value))
  raw <- raw[!is.na(raw$clone), ]
  raw$x <- as.numeric(raw$clone)

  meta <- env$strain_meta(ord)
  rows <- list(meta$mutation, meta$cell, meta$origin)

  lo  <- min(raw$value, est$lo)
  hi  <- max(raw$value, est$hi)
  pad <- hi - lo

  b_start <- 0.10; b_step <- 0.11; b_gap <- 0.035
  frac <- b_start + b_step * (length(rows) - 1) + b_gap

  y_scale <- if (log_y) {
    scale_y_continuous(labels = env$od_labels,
                       sec.axis = dup_axis(labels = env$od_labels))
  } else {
    scale_y_continuous(sec.axis = dup_axis())
  }

  ggplot() +
    geom_jitter(data = raw, aes(x, value), width = 0.13, height = 0,
                colour = "grey78", size = 2.4, alpha = 0.75) +
    (if (ancestor_line)
       geom_hline(yintercept = est$mid[match("ancestor", ord)],
                  linetype = "dashed", colour = "grey35") else NULL) +
    geom_linerange(data = est, aes(x = x, ymin = lo, ymax = hi),
                   linewidth = 0.7, colour = "black") +
    geom_point(data = est, aes(x = x, y = mid), size = 3.6, shape = 21,
               colour = "black", fill = "white", stroke = 0.9) +
    env$nested_brackets(rows, axis = "x", start = lo - b_start * pad,
                        step = -b_step * pad, gap = -b_gap * pad, size = 4) +
    scale_x_continuous(breaks = seq_along(ord), sec.axis = dup_axis()) +
    y_scale +
    labs(y = ylab) +
    coord_cartesian(ylim = c(lo - 0.02 * pad, hi + 0.02 * pad), clip = "off") +
    env$mytheme +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.x.top = element_blank(),
          plot.margin = margin(
            5, 5,
            env$bracket_margin(frac, fig_height, overhead_in = 0.55) + 10, 5))
}

v3_umax <- interval_plot(g, g$post$umax$quantiles$absolute,
                         g$fits$umax, g$fits$clones, g$PANELS$umax$strip_lab)
# No ancestor line on the biofilm panel. The reference well is recorded as
# B. subtilis 168 delta 6 and is probably not this experiment's ancestor, and
# the claim this panel supports -- sinR against ywcC -- does not use it. The
# growth panel keeps its line, where the ancestor is not in doubt.
v3_bio  <- interval_plot(b, b$post$quantiles$absolute,
                         b$assay$OD550_Corrected, b$assay$clones,
                         expression(Biofilm ~ (OD[550])), log_y = TRUE,
                         ancestor_line = FALSE)

## ---- write ------------------------------------------------------------------

pair <- function(p1, p2, collect = FALSE) {
  fig <- (p1 | p2) + plot_annotation(tag_levels = "A") & TAG_THEME
  if (collect) fig <- fig + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  fig
}

save_figure(pair(v1_umax, v1_bio), "figure2_v1_posterior_ridges", width = 13, height = 5.6)
save_figure(pair(v2_umax, v2_bio), "figure2_v2_mean_sd", width = 13, height = 5.4)
save_figure(pair(v3_umax, v3_bio), "figure2_v3_posterior_interval",
            width = 13, height = 5.4)

# Whichever variant FIGURE2_VARIANT names is also written under the name the
# manuscript includes, so main.tex's \includegraphics{figures/figure2_phenotypes.pdf}
# gets it and not the superseded four-panel version from R/02_figure2.R. All
# three are on disk either way, so switching is one word in R/00_setup.R.
chosen <- switch(FIGURE2_VARIANT,
  v1 = list(fig = pair(v1_umax, v1_bio), h = 5.6),
  v2 = list(fig = pair(v2_umax, v2_bio), h = 5.4),
  v3 = list(fig = pair(v3_umax, v3_bio), h = 5.4),
  stop("FIGURE2_VARIANT must be \"v1\", \"v2\" or \"v3\", not \"",
       FIGURE2_VARIANT, "\"", call. = FALSE))
save_figure(chosen$fig, "figure2_phenotypes", width = 13, height = chosen$h)
message(sprintf("Figure 2 in the manuscript is variant %s (%s).", FIGURE2_VARIANT,
                c(v1 = "posterior ridges", v2 = "mean +/- SD",
                  v3 = "posterior interval")[[FIGURE2_VARIANT]]))

## ---- yield and lag, for the supplement --------------------------------------

supp <- (g$ridge_plot(g$post$A$absolute[, g$ABS_ORDER], g$ABS_ORDER,
                      g$PANELS$A$absolute_lab,
                      xlim = abs_xlim_for(g$post$A$quantiles$absolute, g$ABS_ORDER),
                      scale = 2, ref = ancestor_ref(g$post$A$absolute),
                      rows = g$rows_abs) |
         (g$ridge_plot(g$post$L$absolute[, g$ABS_ORDER], g$ABS_ORDER,
                       g$PANELS$L$absolute_lab,
                       xlim = abs_xlim_for(g$post$L$quantiles$absolute, g$ABS_ORDER),
                       scale = 2, ref = ancestor_ref(g$post$L$absolute),
                       rows = NULL) + theme(legend.position = "none"))) +
  plot_annotation(tag_levels = "A") & TAG_THEME

save_figure(supp, "figureS_growth_yield_lag", width = 13, height = 5.2)

cat("Wrote figure2_v1/v2/v3 and figureS_growth_yield_lag.\n")
