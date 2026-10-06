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
#   none    : M23, M26   (spore-fraction clones with no detected mutation)
#
# M23 and M26 are labelled "spormut" in treatments_original_corrected.csv, and
# this code used to read that as "sporulation mutant". They are not mutants:
# sequencing finds no mutation in clones 23 or 26 (2.Mutations output,
# clones_without_mutations.txt), and the lab's own treatments.csv records both
# as mutant = no, mutation = none. The data file is left as it is and the label
# is corrected here, in strain_meta().
#
# MODEL_ORDER is the order the original JAGS code fed to the sampler, kept so
# that column indices stay comparable with the old outputs. PLOT_ORDER is the
# bottom-to-top order of the published panel, which groups clones by mutation.
MODEL_ORDER <- c("ancestor", "M26", "M23", "M54", "M41",
                 "M21", "M19", "M17", "M79", "M13", "M4")
PLOT_ORDER  <- c("M4", "M13", "M79", "M17", "M19",
                 "M21", "M41", "M54", "M23", "M26")

# The six endpoint-spore isolates. Phenotyped in the same five runs as the
# clones above, six curves each, but left out of the published analysis. They
# are the only strains that separate carrying a mutation from the fraction a
# clone was sequenced from: three carry one and three do not, all within the
# spore fraction. See SPORE_FRACTION_MUTANTS below and R/04_group_models.R.
SPORE_ISOLATES <- c("S1", "S6", "S11", "S22", "S51", "S95")

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
#' the ancestor has no cell type or mutation to contrast, and M23 and M26 carry
#' no mutation to label.
#'
#' `mutated` is TRUE for the eight clones carrying a sinR or ywcC mutation and
#' FALSE for the ancestor, M23 and M26. It defines model 3.
# The six endpoint-spore isolates, and which of them carry a mutation.
#
# treatments_original_corrected.csv records all six as "nomut", following
# treatments_original.csv (the plate map), which also files them under
# evo.type "anc". The per-replicate file treatments.csv disagrees: S1 and S6
# carry "yetA.ymlG" and S95 "levB_yveA", with S11, S22 and S51 carrying none.
#
# The 96-clone endpoint-spore variant matrix (Megan Behringer, Oct 2026;
# 2.Mutations/data/LT_Heat_Bacillus.compare.tab.xlsx) settles which isolate
# carries what. Four mutations survive QC, in three clones:
#   clone 1  (S1)   yetA 777,008 and ylmG 1,611,379
#   clone 6  (S6)   kdgA 2,323,251
#   clone 95 (S95)  levB-aspP 3,539,121
# treatments.csv records S6 as "yetA.ymlG", copied from S1; the matrix says
# kdgA. All three mutant clones in that fraction were phenotyped.
#
# Corrected here rather than in the CSV, as everywhere else in this project.
SPORE_FRACTION_MUTANTS <- c(S1 = "yetA, ylmG", S6 = "kdgA", S95 = "levB-aspP")

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
                      ifelse(clones %in% TOTAL_FRACTION_NO_MUTATION, "Total",
                             ifelse(t$cell.type == "spore", "Spore", "Total"))),
    # M23 and M26 carry no mutation, so no label; the cell-type row already
    # marks them as Spore.
    mutation = ifelse(clones %in% names(SPORE_FRACTION_MUTANTS),
                      SPORE_FRACTION_MUTANTS[clones],
                      gene_display(c(nomut = NA, sinR = "sinR", ywcC = LINEAGE_YWCC,
                                     spormut = NA)[t$mutation])),
    # Full mutation factor, used by the models. "spormut" is the data file's
    # label for M23/M26 and is shown as "none", which is what they carry.
    mutation_full = ifelse(clones %in% names(SPORE_FRACTION_MUTANTS),
                           SPORE_FRACTION_MUTANTS[clones],
                           gene_display(c(nomut = "none", sinR = "sinR",
                                          ywcC = LINEAGE_YWCC,
                                          spormut = "none")[t$mutation])),
    # "carries a mutation", not "carries a biofilm-regulator mutation". On the
    # eleven strains the analysis uses, the two are the same set, so this does
    # not change model 3; it matters only when the spore isolates are included.
    mutated = t$mutation %in% c("sinR", "ywcC") |
              clones %in% names(SPORE_FRACTION_MUTANTS),
    stringsAsFactors = FALSE
  ) -> m
  # each clone's own acquired mutations, for per-strain tables
  m$genes <- ifelse(m$clone %in% names(LINEAGE_ONLY_GENES), LINEAGE_ONLY_GENES[m$clone],
                    ifelse(m$mutation_full == LINEAGE_YWCC, "ywcC", m$mutation_full))
  m
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
