################################################################################
# 00_setup.R -- paths, packages, plotting theme, strain metadata, helpers
#
# Sourced by every other script; not meant to be run on its own.
################################################################################

## ---- project root ----------------------------------------------------------
# Walk up from the working directory until we find the project marker, so the
# scripts run the same whether you are in the project root or in R/.
find_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "data", "biofil.csv"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not locate the biofilm project root.")
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
#   ywcC      : m4, m13, m79   (the data file calls this ywcC/slrR)
#   sinR      : m17, m19, m21, m41, m54
#   spore     : m23, m26   (sporulation mutants, assayed as spores)
#
# MODEL_ORDER is the row order of data/biofil.csv. PLOT_ORDER is the
# bottom-to-top order of the published ridgeline panel, and STRIP_ORDER the
# left-to-right order of the strain panel; both group clones by mutation.
MODEL_ORDER <- c("m4", "m13", "m17", "m19", "m21",
                 "m23", "m26", "m41", "m54", "m79", "ancestor")
PLOT_ORDER  <- c("m4", "m13", "m79", "m17", "m19",
                 "m21", "m41", "m54", "m23", "m26")
STRIP_ORDER <- c("ancestor", "m23", "m26", "m17", "m19",
                 "m21", "m41", "m54", "m4", "m13", "m79")

N_REPS <- 8L   # replicate wells per strain (A-H on the single plate)

## ---- gene naming -----------------------------------------------------------
# data/biofil.csv labels the m4/m13/m79 mutation "ywcC/slrR" and the assay file
# is left exactly as it is; the mapping below is for display only.
#
# BSU_38220 is displayed as ywcC, not as the current RefSeq symbol slrC. The
# rename is genuine -- it sits beside slrA and its product is the regulator of
# it -- but the people who know this regulatory network know it as ywcC, and
# slrC sits one letter from slrA and slrR, which appear in the same sentences.
# The locus tag carries the identification. Cite as "ywcC (BSU_38220)" at
# first mention.
#
# The "/slrR" half of that label is NOT the slrR gene. An early manuscript
# draft names the group "ywcC/epsA-slrR mutants" and describes the split as
# "clones containing a mutation in sinR did not also contain mutations in ywcC
# or upstream of slrR" -- so it means the intergenic epsA-slrR variant at
# 3,529,981, the marker that splits the sequenced clones into two clades.
# R/06_lineage_groups.R runs that comparison.
# M23 and M26 are clones 23 and 26 of the endpoint TOTAL-fraction matrix, and
# carry no mutation there. They were labelled "spore" because they had no
# mutation -- inferring cell type from mutation status, which then makes any
# test of mutation against cell type circular. Megan Behringer confirmed (Oct
# 2026) they came from the mixed fraction. Their fraction is what is known, so
# they are labelled Total. The S isolates are genuinely spore-fraction: their
# mutations (yetA, ylmG, kdgA, levB-aspP) are the endpoint-spore calls and
# appear nowhere in the total-fraction matrix.
TOTAL_FRACTION_NO_MUTATION <- c("M23", "M26")

# The comparison group {M4, M13, M79} is a LINEAGE group, not a gene group.
# M4 and M13 carry ywcC frameshifts; M79 carries none -- its acquired mutations
# are yutK and cotI, on the epsA-slrR standing variant. The original metadata
# called the group "ywcC.slrR" (growth) and "ywcC/slrR" (biofilm), meaning
# "ywcC or the slrR variant"; the code had shortened that to "ywcC", which made
# M79 look like a ywcC mutant. Group comparisons use the lineage name, and the
# per-strain tables give each clone its own genes (column `genes`).
LINEAGE_YWCC <- "ywcC/epsA-slrR"
LINEAGE_ONLY_GENES <- c(M79 = "yutK, cotI")


# biofil.csv puts "spore" in the mutation column for m23 and m26. That is the
# fraction they came from, not a mutation: sequencing finds none in clones 23
# or 26. They are labelled "none" for the models, and the ancestor likewise.
MUTATION_LABEL <- c(ancestor = "none", sinR = "sinR",
                    `ywcC/slrR` = LINEAGE_YWCC, spore = "none")

#' Grouping variables for a set of strains, in the order given.
#'
#' Read from the assay file itself, which carries the treatment and mutation of
#' every well. `origin`, `cell` and `mutation` are the three nested annotation
#' rows of the figures. NA means no bracket is drawn for that strain in that
#' row: the ancestor has no cell type or mutation to contrast, and the
#' sporulation mutants are already marked Spore by the cell-type row.
strain_meta <- function(clones) {
  d <- read.csv(file.path(DATA_DIR, "biofil.csv"), stringsAsFactors = FALSE)
  d <- d[match(clones, d$clones), ]
  if (anyNA(d$clones)) stop("Unknown strain: ",
                            paste(clones[is.na(d$clones)], collapse = ", "))

  data.frame(
    clone    = clones,
    label    = clones,
    origin   = ifelse(d$Treatment == "Ancestor", "Ancestor", "Evolved"),
    # biofil.csv's Treatment column says "Vegetative", which describes how the
    # strain was assayed. The growth-curve project labels the same ten clones
    # "Total", which is the sequencing fraction they were isolated from. The
    # partition is identical either way -- the two spore-fraction clones against
    # the other eight -- so the figures use one label, and the fraction is the
    # one that is a property of the clone rather than of the assay. The assay
    # condition belongs in the Methods, not on a bracket. The data file is
    # left alone.
    cell     = ifelse(d$Treatment == "Ancestor", NA,
                      ifelse(toupper(d$clones) %in% TOTAL_FRACTION_NO_MUTATION, "Total",
                             ifelse(d$Treatment == "Vegetative", "Total", d$Treatment))),
    mutation = c(ancestor = NA, sinR = "sinR", `ywcC/slrR` = LINEAGE_YWCC,
                 spore = NA)[d$mutation],
    # full mutation factor, used by the models rather than the figures
    mutation_full = MUTATION_LABEL[d$mutation],
    # TRUE for the eight sinR / ywcC clones; defines model 3
    mutated = d$mutation %in% c("sinR", "ywcC/slrR"),
    stringsAsFactors = FALSE
  ) -> m
  m$genes <- ifelse(toupper(m$clone) %in% names(LINEAGE_ONLY_GENES),
                    LINEAGE_ONLY_GENES[toupper(m$clone)],
                    ifelse(m$mutation_full == LINEAGE_YWCC, "ywcC", m$mutation_full))
  m
}

## ---- the ancestor is under question ---------------------------------------
# The plate note in OneDrive.../Biofilm/Behringer_Plate1.xlsx calls column 11
# "B. subtilis 168 delta 6", where data/biofil.csv calls it "Ancestor". Delta 6
# is a domesticated, biofilm-impaired 168 derivative, not the strain the
# starvation experiment began from, so everything expressed *relative to the
# ancestor* is provisional. What does not depend on it: the absolute posteriors,
# the lineage-group comparison in R/06_lineage_groups.R, and every model fit to
# the evolved clones alone.
ANCESTOR_IS_PROVISIONAL <- TRUE

warn_ancestor <- function() {
  if (!ANCESTOR_IS_PROVISIONAL) return(invisible(NULL))
  message(strrep("-", 72))
  message("NOTE: the 'ancestor' in data/biofil.csv is recorded on the original plate")
  message("      as B. subtilis 168 delta 6, which is probably not this experiment's")
  message("      ancestor. Everything relative to it is provisional. Absolute values,")
  message("      the lineage-group comparison and the evolved-only models are fine.")
  message(strrep("-", 72))
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

# Axis labels for strains, should a figure want them.
short_label <- function(clone) ifelse(clone == "ancestor", "anc", clone)

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

#' Nested group brackets drawn outside the panel.
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
                            italic = c("sinR", "ywcC", "slrC", LINEAGE_YWCC)) {
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
