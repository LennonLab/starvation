################################################################################
# 00_setup.R -- paths, packages, plotting theme, small helpers
#
# Sourced by every other script; not meant to be run on its own.
################################################################################

## ---- project root ----------------------------------------------------------
# Walk up from the working directory until we find the project marker, so the
# scripts run the same whether you are in the project root or in R/.
find_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "data", "comp_data.csv"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not locate the growthCurves project root.")
    path <- parent
  }
}

PROJ     <- find_root()
DATA_DIR <- file.path(PROJ, "data")
OUT_DIR  <- file.path(PROJ, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- packages --------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
})

## ---- clone metadata --------------------------------------------------------
# The ancestor is the reference; the ten evolved clones are grouped by the
# mutation they carry:
#   ywcC    : M4, M13, M79   (BSU_38220; RefSeq now calls it slrC)
#   sinR    : M17, M19, M21, M41, M54
#   spormut : M23, M26   (sporulation mutants, scored as spores)
#
# MODEL_ORDER is the order the original JAGS code fed to the sampler, kept so
# that column indices stay comparable with the old outputs. PLOT_ORDER is the
# bottom-to-top order of the published panel, which groups clones by mutation.
MODEL_ORDER <- c("ancestor", "M26", "M23", "M54", "M41",
                 "M21", "M19", "M17", "M79", "M13", "M4")
PLOT_ORDER  <- c("M4", "M13", "M79", "M17", "M19",
                 "M21", "M41", "M54", "M23", "M26")

N_REPS <- 6L   # replicate growth curves retained per clone (the lowest-RMSE six)

## ---- gene naming -----------------------------------------------------------
# BSU_38220 is displayed as ywcC, not as the current RefSeq symbol slrC.
# The rename is genuine, but the people who know this regulatory network know
# the gene as ywcC, and slrC sits one letter from slrA and slrR, which appear
# in the same sentences. The locus tag carries the identification; the symbol
# only has to be recognisable. Cite as "ywcC (BSU_38220)" at first mention.
# The data files already say ywcC, so nothing needs mapping.
GENE_SYNONYMS <- character(0)

gene_display <- function(x) {
  ifelse(is.na(x), x, ifelse(x %in% names(GENE_SYNONYMS), GENE_SYNONYMS[x], x))
}

# Left-to-right order of the manuscript's strain panel: ancestor, then evolved
# clones grouped by mutation.
STRIP_ORDER <- c("ancestor", "M23", "M26", "M17", "M19",
                 "M21", "M41", "M54", "M4", "M13", "M79")

#' Grouping variables for a set of strains, in the order given.
#'
#' `origin`, `cell`, and `mutation` are the three nested annotation rows of the
#' manuscript figure. NA means no bracket is drawn for that strain in that row:
#' the ancestor has no cell type or mutation to contrast, and the sporulation
#' mutants are not part of the sinR/ywcC comparison.
strain_meta <- function(clones) {
  t <- read.csv(file.path(DATA_DIR, "treatments_original_corrected.csv"),
                stringsAsFactors = FALSE)
  t <- t[match(clones, t$clones), ]
  if (anyNA(t$clones)) stop("Unknown strain: ",
                            paste(setdiff(clones, t$clones), collapse = ", "))

  data.frame(
    clone    = clones,
    label    = tolower(clones),
    origin   = ifelse(t$evo.type == "ancestor", "Ancestor", "Evolved"),
    cell     = ifelse(t$evo.type == "ancestor", NA,
                      ifelse(t$cell.type == "spore", "Spore", "Total")),
    # The sporulation mutants get no label here: the cell-type row already
    # marks them as Spore.
    mutation = gene_display(c(nomut = NA, sinR = "sinR", ywcC = "ywcC",
                              spormut = NA)[t$mutation]),
    # full mutation factor, used by the models rather than the figures; also
    # mapped to current names so group labels and contrasts agree with the
    # figures
    mutation_full = gene_display(t$mutation),
    stringsAsFactors = FALSE
  )
}

## ---- plotting theme (Lennon lab house style) -------------------------------
mytheme <- theme_bw() +
  theme(
    axis.ticks.length = unit(0.25, "cm"),
    axis.text         = element_text(size = 14),
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),
    legend.text       = element_text(size = 12),
    legend.title      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(fill = NA, colour = "black", linewidth = 1),
    strip.text.x      = element_text(size = 14),
    # secondary axes are ticks only -- no duplicated labels
    axis.text.x.top   = element_blank(), axis.title.x.top   = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_blank()
  )

## ---- nested group brackets -------------------------------------------------

# Axis labels: lower case, and "ancestor" shortened the way the manuscript
# panel writes it.
short_label <- function(clone) ifelse(clone == "ancestor", "anc", tolower(clone))

# Space to reserve outside the panel for `n` bracket rows, as a margin in
# points. `frac` is how far the brackets reach expressed as a fraction of the
# panel's own extent, and the margin eats into that extent, hence the solve.
bracket_margin <- function(frac, size_in, overhead_in) {
  72 * frac * (size_in - overhead_in) / (1 + frac)
}

# Consecutive runs of the same non-NA label, as start/end positions.
label_runs <- function(x) {
  keep <- !is.na(x) & nzchar(x)
  if (!any(keep)) return(NULL)
  i <- which(keep)
  brk <- c(0, which(diff(i) != 1 | x[i][-1] != x[i][-length(i)]), length(i))
  do.call(rbind, lapply(seq_len(length(brk) - 1), function(k) {
    idx <- i[(brk[k] + 1):brk[k + 1]]
    data.frame(label = x[idx[1]], start = min(idx), end = max(idx))
  }))
}

#' Nested group brackets drawn outside the panel, as in the manuscript figure.
#'
#' Returns a list of ggplot layers. The plot must use coord_cartesian(clip =
#' "off") and leave margin on the bracket side for them to be visible.
#'
#' @param rows   list of label vectors, innermost row first; one entry per
#'               strain position, NA where no bracket belongs
#' @param axis   "x" draws brackets below the panel, "y" to the left of it
#' @param start  coordinate (on the other axis) of the first bracket row
#' @param step   signed distance between successive rows, pointing away
#' @param gap    signed distance from a bracket to its label
#' @param italic labels to set in italic (gene names)
nested_brackets <- function(rows, axis = c("x", "y"), start, step, gap,
                            size = 4.2, colour = "grey30",
                            italic = c("sinR", "ywcC", "slrC")) {
  axis   <- match.arg(axis)
  layers <- list()

  for (r in seq_along(rows)) {
    runs <- label_runs(rows[[r]])
    if (is.null(runs)) next
    off <- start + (r - 1) * step
    lo  <- runs$start - 0.45
    hi  <- runs$end   + 0.45
    mid <- (runs$start + runs$end) / 2

    # Caps turn each span into a bracket, so it is obvious where one group
    # ends and the next begins. They point back towards the panel.
    cap <- -sign(step) * abs(step) * 0.22

    layers <- c(layers, list(
      if (axis == "x") {
        annotate("segment", x = lo, xend = hi, y = off, yend = off,
                 colour = colour, linewidth = 0.7)
      } else {
        annotate("segment", y = lo, yend = hi, x = off, xend = off,
                 colour = colour, linewidth = 0.7)
      },
      if (axis == "x") {
        annotate("segment", x = c(lo, hi), xend = c(lo, hi),
                 y = off, yend = off + cap, colour = colour, linewidth = 0.7)
      } else {
        annotate("segment", y = c(lo, hi), yend = c(lo, hi),
                 x = off, xend = off + cap, colour = colour, linewidth = 0.7)
      }
    ))

    for (face in c("plain", "italic")) {
      sel <- if (face == "italic") runs$label %in% italic else !runs$label %in% italic
      if (!any(sel)) next
      layers <- c(layers, list(
        if (axis == "x") {
          annotate("text", x = mid[sel], y = off + gap, label = runs$label[sel],
                   colour = colour, size = size, fontface = face, vjust = 1)
        } else {
          annotate("text", y = mid[sel], x = off + gap, label = runs$label[sel],
                   colour = colour, size = size, fontface = face, angle = 90)
        }
      ))
    }
  }
  layers
}

## ---- helpers ---------------------------------------------------------------
sem   <- function(x) sqrt(var(x) / length(x))
cv    <- function(x) 100 * (sd(x) / mean(x))
LL.95 <- function(x) t.test(x)$conf.int[1]
UL.95 <- function(x) t.test(x)$conf.int[2]
