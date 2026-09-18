################################################################################
# 02_figure2.R
#
# Figure 2. Phenotypes of the evolved clones, grouped by the mutation they
# carry -- the two genes Figure 1C nominates.
#
#   A  maximum growth rate   3.GrowthCurves/R/02_figures.R
#   B  maximum yield         3.GrowthCurves/R/02_figures.R
#   C  lag time              3.GrowthCurves/R/02_figures.R
#   D  biofilm               4.Biofilm/R/02_figures.R
#
# Output: output/figure2_phenotypes.pdf / .png
#         output/figure2_caption_notes.md
#
# The panels are the projects' own ridge_plot() output, built from the same
# posterior draws the project figures use -- see load_project() in R/00_setup.R.
# Every strain, ordering and bracket decision is inherited. What is decided
# here: which scale (absolute, see PHENOTYPE_SCALE in R/00_setup.R), that
# group brackets are drawn once per row rather than on all four panels, and
# the ancestor footnote described below.
################################################################################

.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the 5.MsFigures project root (or R/).")
source(.setup[1])

stopifnot(PHENOTYPE_SCALE %in% c("absolute", "relative"))

## ---- growth parameters ------------------------------------------------------

g <- load_project("growth", "02_figures.R")
require_objects(g, "ridge_plot", "post", "ABS_ORDER", "rows_abs", "rows_rel",
                "PLOT_ORDER", "PANELS", what = "3.GrowthCurves")

#' One absolute ridge panel, with the project's own axis-range rule.
#'
#' The range comes from the credible intervals rather than the extreme
#' quantiles, because the vague Gamma prior on the precision leaves long tails
#' on strains whose replicates disagree. This is the rule 3.GrowthCurves uses;
#' it is repeated here only because the project applies it inline.
abs_xlim_for <- function(quant, order, pad = 0.12) {
  r <- range(quant[order, c("2.5%", "97.5%")])
  r + c(-1, 1) * pad * diff(r)
}

# On the measured scale the panels have no natural reference, which makes a
# clone hard to place against the strain it descends from. The dashed line is
# the ancestor's posterior mean -- the same quantity the relative panels divide
# by, drawn rather than divided out. On the relative panels it stays at 1,
# which is where the ancestor already is.
ancestor_ref <- function(draws) mean(draws[, "ancestor"])

growth_panel <- function(nm, with_brackets) {
  cfg <- g$PANELS[[nm]]
  if (PHENOTYPE_SCALE == "absolute") {
    g$ridge_plot(g$post[[nm]]$absolute[, g$ABS_ORDER], g$ABS_ORDER,
                 cfg$absolute_lab,
                 xlim = abs_xlim_for(g$post[[nm]]$quantiles$absolute, g$ABS_ORDER),
                 scale = 2, ref = ancestor_ref(g$post[[nm]]$absolute),
                 rows = if (with_brackets) g$rows_abs else NULL)
  } else {
    g$ridge_plot(g$post[[nm]]$relative, g$PLOT_ORDER, cfg$relative_lab,
                 xlim = cfg$relative_xlim, scale = cfg$relative_scale, ref = 1,
                 rows = if (with_brackets) g$rows_rel else NULL)
  }
}

panel_a <- growth_panel("umax", with_brackets = TRUE)
panel_b <- growth_panel("A",    with_brackets = FALSE)
panel_c <- growth_panel("L",    with_brackets = TRUE)

## ---- biofilm ----------------------------------------------------------------

b <- load_project("biofilm", "02_figures.R")
require_objects(b, "ridge_plot", "post", "ABS_ORDER", "rows_abs", "rows_rel",
                "PLOT_ORDER", "od_labels", "ANCESTOR_IS_PROVISIONAL",
                what = "4.Biofilm")

# data/biofil.csv records the ancestor as B. subtilis 168 delta 6, which is
# probably not this experiment's ancestor. The absolute reading is a real
# measurement of a real strain, so the row stays; what cannot stand in a
# manuscript figure is calling it the ancestor without qualification.
#
# The mark goes on this panel's own caption rather than on the shared group
# brackets. The brackets are drawn once per row and are read by the growth
# panel beside it, whose ancestor is not in doubt -- marking them there would
# say the wrong thing about the wrong measurement. Resolve the strain, set
# ANCESTOR_IS_PROVISIONAL to FALSE in 4.Biofilm, and the mark disappears.
# The note names the row rather than carrying a marker into the panel: the
# group brackets are shared with the growth panel beside it, so there is
# nowhere inside this panel to hang a marker that would not also read as a
# statement about the lag-time ancestor.
biofilm_caption <- if (isTRUE(b$ANCESTOR_IS_PROVISIONAL) &&
                       isTRUE(SHOW_ANCESTOR_FOOTNOTE)) {
  paste0(FOOTNOTE_MARK, " biofilm ancestor (bottom row): strain unconfirmed")
} else NULL

panel_d <- if (PHENOTYPE_SCALE == "absolute") {
  # The biofilm panel is drawn on log10 draws, so the reference is the mean of
  # the log10 ancestor draws, not the log of the mean.
  b$ridge_plot(log10(b$post$absolute[, b$ABS_ORDER]), b$ABS_ORDER,
               xlab = expression(Biofilm ~ (OD[550])),
               xlim = abs_xlim_for(log10(b$post$quantiles$absolute), b$ABS_ORDER),
               scale = 2, ref = ancestor_ref(log10(b$post$absolute)),
               rows = NULL, labels_fn = b$od_labels)
} else {
  b$ridge_plot(log10(b$post$relative), b$PLOT_ORDER,
               xlab = expression(paste("Relative biofilm, ", italic(log[10]))),
               xlim = c(-1.5, 1.5), scale = 2, ref = 0, rows = NULL)
}

if (!is.null(biofilm_caption) && PHENOTYPE_SCALE == "absolute") {
  panel_d <- panel_d +
    labs(caption = biofilm_caption) +
    theme(plot.caption = element_text(size = 8.5, colour = "grey30", hjust = 0))
}

## ---- assemble ---------------------------------------------------------------

# ridge_plot() puts the tail-probability key inside the panel, where it works
# on a figure of one panel and collides with the ridges on a figure of four.
# The four keys are identical, so they are moved out and collected into one.
#' Prepare one panel for the 2x2 grid.
#'
#' `pad_left` separates the two columns. It is applied only to the right-hand
#' panels, which carry no brackets and so have ridge_plot()'s default left
#' margin of 5pt; the left-hand panels keep the margin ridge_plot() computed
#' for their brackets, which setting plot.margin here would discard.
for_grid <- function(p, pad_left = NULL) {
  p <- p + theme(legend.position   = "bottom",
                 legend.direction  = "horizontal",
                 legend.title      = element_text(size = 11.5, hjust = 0.5),
                 legend.text       = element_text(size = 9.5),
                 legend.title.position = "top",
                 legend.key.height = unit(0.3, "cm"),
                 legend.key.width  = unit(1.6, "cm"))
  if (!is.null(pad_left)) {
    p <- p + theme(plot.margin = margin(5, 10, 5, pad_left))
  }
  p
}

PAD_LEFT <- 30   # pt between the columns

figure2 <- ((for_grid(panel_a) | for_grid(panel_b, PAD_LEFT)) /
            (for_grid(panel_c) | for_grid(panel_d, PAD_LEFT))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom") &
  TAG_THEME

save_figure(figure2, "figure2_phenotypes", width = 13, height = 8.8)

## ---- caption notes ----------------------------------------------------------

notes <- c(
  "# Figure 2 — notes for the caption",
  "",
  sprintf("Scale: **%s**.", PHENOTYPE_SCALE),
  "",
  "Ridge heights are posterior densities; shading is the two-tailed tail",
  "probability, so the pale ends of a ridge are its tails and the dark middle",
  "its bulk. Brackets group the clones by the mutation they carry, by cell",
  "type, and by origin; they are drawn once per row and apply to both panels",
  "in that row.",
  "")

if (isTRUE(b$ANCESTOR_IS_PROVISIONAL)) {
  notes <- c(notes,
    if (!isTRUE(SHOW_ANCESTOR_FOOTNOTE))
      paste("The caveat below is **not printed on the figure**",
            "(`SHOW_ANCESTOR_FOOTNOTE` is FALSE in `R/00_setup.R`).",
            "It still applies.") else NULL,
    if (!isTRUE(SHOW_ANCESTOR_FOOTNOTE)) "" else NULL,
    # the mark is escaped: an unescaped "*" abutting "**" is read as emphasis
    sprintf("%s **The biofilm ancestor is provisional.**",
            gsub("*", "\\*", FOOTNOTE_MARK, fixed = TRUE)),
    "`4.Biofilm/data/biofil.csv`",
    "labels the reference well \"Ancestor\", but the original plate record",
    "(Behringer_Plate1.xlsx, column 11) names it *B. subtilis* 168 Δ6, which",
    "is probably not this experiment's ancestor. The absolute reading is a",
    "genuine measurement, so the row is shown; the mark notes that the",
    "strain's identity as the ancestor is unconfirmed. The evolved clones and",
    "every comparison among them are unaffected. See `4.Biofilm/README.md`.",
    "")
}

writeLines(notes, file.path(OUT_DIR, "figure2_caption_notes.md"))

cat("Wrote output/figure2_phenotypes.pdf / .png and figure2_caption_notes.md\n")
