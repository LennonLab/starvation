################################################################################
# 00_setup.R -- paths, packages, genome constants, the allele-frequency palette
#
# Sourced by every other script; not meant to be run on its own.
################################################################################

find_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, "data")) &&
        file.exists(file.path(path, "R", "00_setup.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not locate the mutations project root.")
    path <- parent
  }
}

PROJ     <- find_root()
DATA_DIR <- file.path(PROJ, "data")
OUT_DIR  <- file.path(PROJ, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(viridisLite)
  library(tidyr)
  library(ggplot2)
})

## ---- genome ----------------------------------------------------------------
GENOME     <- "NC_000964.3"          # Bacillus subtilis subsp. subtilis 168
GENOME_BP  <- 4215606L

## ---- the three sequenced fractions -----------------------------------------
# n is the number of clones sequenced in each, which is the denominator for
# every frequency reported.
FRACTIONS <- list(
  S10      = list(label = "Early spore fraction",    n = 52L),
  S1000    = list(label = "Endpoint spore fraction", n = 96L),
  `S+NS1000` = list(label = "Endpoint total fraction", n = 84L)
)
FRACTION_ORDER <- names(FRACTIONS)

# A variant carried by roughly half the clones is not something that arose
# during the experiment; it is excluded from the figure and every count, as it
# was in the published version.
EXCLUDE_POS <- 3529981L

# The Circos highlight tracks (*.mutations.txt, *.circos.txt) do not store
# point positions: each row is a window of +/-10 kb around the mutation, so
# start = position - 10000 and end = position + 10000. Anything read from a
# track must be shifted by this. Enforced by check_circos_offset().
CIRCOS_HALF_WINDOW <- 10000L

## ---- allele-frequency palette ----------------------------------------------
# The published Circos track used a rainbow, which has no perceptual order --
# a reader cannot tell from the colours alone which of two tiles is the more
# frequent. Magma is ordered and stays legible in greyscale: singletons dark,
# the most widespread mutation bright.
#
# The original mapping is kept in FREQ_COLOURS_ORIGINAL, since it is what
# decodes the archived track files.
FREQ_COLOURS_ORIGINAL <- c(
  "1" = "#0000E5", "2" = "#00BFFF", "3" = "#ADFF2F",
  "4" = "#FFD700", "6" = "#FFA500", "8" = "#FF0000")

MAGMA_RANGE <- c(0.12, 0.86)   # skip the near-black and near-white extremes

freq_colour <- function(count, max_count = 8) {
  x <- pmin(pmax((count - 1) / (max_count - 1), 0), 1)
  viridisLite::magma(256, begin = MAGMA_RANGE[1], end = MAGMA_RANGE[2])[
    1 + round(x * 255)]
}

# Ring background tints, outermost first, so each ring can be keyed to its
# label in the centre.
RING_TINTS <- c("#E2E2E2", "#ECECEC", "#F6F6F6")

## ---- gene names ------------------------------------------------------------
# Names are resolved once, here, so every table and figure agrees. Mutations
# are mapped to the CDS that contains them rather than matched by name, and the
# symbol is taken from the RefSeq annotation of NC_000964.3 except where
# data/gene_name_overrides.csv says otherwise -- see that file for why.
#
# Convention: current name everywhere, with the legacy name carried in its own
# column so it can be given in parentheses at first mention.

load_cds <- function() {
  cds <- read.delim(file.path(DATA_DIR, "NC_000964.3_cds.tsv"),
                    stringsAsFactors = FALSE)
  ov <- read.csv(file.path(DATA_DIR, "gene_name_overrides.csv"),
                 comment.char = "#", stringsAsFactors = FALSE)
  i <- match(cds$locus_tag, ov$locus_tag)
  cds$refseq_gene <- cds$gene
  cds$gene   <- ifelse(is.na(i), cds$gene, ov$preferred[i])
  cds$legacy <- ifelse(is.na(i), NA_character_, ov$legacy[i])
  cds
}

#' Resolve each position to its containing CDS.
#'
#' @return data frame with gene (current name), legacy, locus_tag, cds_length;
#'   all NA for an intergenic position.
resolve_genes <- function(positions, cds = load_cds()) {
  idx <- vapply(positions, function(p) {
    i <- which(cds$start <= p & cds$end >= p)
    if (length(i)) i[1] else NA_integer_
  }, integer(1))

  data.frame(
    gene_current = cds$gene[idx],
    gene_legacy  = cds$legacy[idx],
    locus_tag    = cds$locus_tag[idx],
    cds_length   = cds$length_bp[idx],
    stringsAsFactors = FALSE)
}

## ---- plotting theme --------------------------------------------------------
mytheme <- theme_bw() +
  theme(
    axis.ticks.length = unit(0.2, "cm"),
    axis.text         = element_text(size = 11),
    axis.title        = element_text(size = 12),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.8),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 11, hjust = 0),
    legend.title      = element_blank()
  )
