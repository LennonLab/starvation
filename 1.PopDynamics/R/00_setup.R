################################################################################
# 00_setup.R -- paths, packages, plotting theme, helpers
#
# Sourced by every other script; not meant to be run on its own.
################################################################################

## ---- project root ----------------------------------------------------------
find_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "data", "spore.transition.txt"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the populationDynamics project root.")
    }
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
})

## ---- experiment constants --------------------------------------------------
REPS <- c("A", "B", "C", "D")

# Counts are colonies on a plate from a 100 uL spread of a 10^-dil dilution.
# CFU/mL = colonies * 10^dil * 10, the 10 converting 100 uL to 1 mL.
PLATE_VOLUME_FACTOR <- 10

# The latent-state model was fit on counts divided by 10^4, to keep the numbers
# small enough not to lean on the gamma priors. Multiply back to get CFU/mL.
MODEL_SCALE <- 1e4

# The state model splits the series into a transition phase and an equilibrium
# phase at 240 h (10 d).
PHASE_BREAK_H <- 240

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
    axis.text.x.top   = element_blank(), axis.title.x.top   = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_blank()
  )

## ---- helpers ---------------------------------------------------------------
sem <- function(x) sqrt(var(x) / length(x))

# The sigmoidal curve the abundances are fit with, on log10 axes:
#   log10(N) = b + (a - b) / (1 + exp((m - log10(t)) / w))
# a is the level it settles at, b the level it starts from, m the midpoint in
# log10 days and w the width of the transition. At t = 0, log10(t) is -Inf and
# the curve returns b, so day-zero counts pull on the starting asymptote.
#
# Written with plogis rather than the expression above, which is the same
# function: 1 / (1 + exp((m - x) / w)) is plogis((x - m) / w). The direct form
# evaluates exp(Inf) at the day-zero rows, and although the value it returns is
# right, the Inf makes the finite-difference gradient unreliable and L-BFGS-B
# ends the pooled total fit with "ABNORMAL_TERMINATION_IN_LNSRCH". plogis
# handles an infinite argument without overflowing, and the estimates are the
# same to four decimal places.
sigmoid_log10 <- function(log_t, a, b, m, w) {
  b + (a - b) * stats::plogis((log_t - m) / w)
}
