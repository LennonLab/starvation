################################################################################
# 00_setup.R
#
# Shared setup for the two manuscript figures.
#
# The panels are not redrawn here. Each one is built by the code that already
# owns it in 1.PopDynamics, 2.Mutations, 3.GrowthCurves or 4.Biofilm, loaded
# through load_project() below. Nothing about a panel's statistics, ordering
# or annotation is restated in this project, so the manuscript figure and the
# project figure cannot drift apart. What this project decides is layout,
# lettering, and the few presentation choices listed under "House choices".
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

for (p in c("cowplot", "gridGraphics")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. install.packages(\"%s\")", p, p))
  }
}

## ---- where things are ------------------------------------------------------

MS_ROOT <- normalizePath(
  Filter(dir.exists, c(".", ".."))[
    which(vapply(Filter(dir.exists, c(".", "..")),
                 function(d) dir.exists(file.path(d, "R")) &&
                             file.exists(file.path(d, "R", "00_setup.R")),
                 logical(1)))[1]])

STARVATION <- normalizePath(file.path(MS_ROOT, ".."))
OUT_DIR    <- file.path(MS_ROOT, "output")
dir.create(OUT_DIR, showWarnings = FALSE)

PROJECTS <- c(pop     = "1.PopDynamics",
              mut     = "2.Mutations",
              growth  = "3.GrowthCurves",
              biofilm = "4.Biofilm")

project_dir <- function(key) {
  d <- file.path(STARVATION, PROJECTS[[key]])
  if (!dir.exists(d)) stop("Project not found: ", d)
  d
}

## ---- running a project's figure script without letting it write ------------

#' Source one of the projects' scripts into an isolated environment.
#'
#' The scripts build their plot objects and then save them. Here the saving
#' calls are stubbed out, so what comes back is the environment holding the
#' objects -- the panel functions, the posterior draws, the bracket rows --
#' with no files written and nothing added to the global environment.
#'
#' `source` is stubbed too: the scripts locate and source their own
#' `00_setup.R` with the default `local = FALSE`, which would evaluate it in
#' the global environment and let two projects overwrite each other's
#' PLOT_ORDER and DATA_DIR. Redirecting it into `env` keeps each project's
#' constants with that project.
#'
#' @param key   name in PROJECTS
#' @param script file name inside the project's R/ directory
#' @return the environment the script was evaluated in
load_project <- function(key, script) {
  dir <- project_dir(key)
  path <- file.path(dir, "R", script)
  if (!file.exists(path)) stop("No such script: ", path)

  env <- new.env(parent = globalenv())
  env$source  <- function(file, ...) sys.source(file, envir = env)
  env$ggsave  <- function(...) invisible(NULL)
  env$pdf     <- function(...) invisible(NULL)
  env$png     <- function(...) invisible(NULL)
  env$dev.off <- function(...) invisible(NULL)
  env$message <- function(...) invisible(NULL)

  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)

  # With pdf()/dev.off() stubbed, a script that draws with base graphics --
  # 2.Mutations/R/02_circos.R does -- has no open device and R opens the
  # default one, leaving a stray Rplots.pdf in the working directory. A null
  # device absorbs the drawing instead. The panel itself is replayed later by
  # cowplot::as_grob, which opens its own device.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Some project scripts print tables as they go. That is right when you run
  # the project and noise when you are only borrowing a figure, so stdout is
  # captured. Errors and warnings are untouched and still propagate.
  utils::capture.output(sys.source(path, envir = env))
  env
}

#' Stop early and loudly if a project script did not produce what we expect.
require_objects <- function(env, ..., what) {
  need <- c(...)
  miss <- need[!vapply(need, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss)) {
    stop(sprintf("%s: expected object(s) not found: %s",
                 what, paste(miss, collapse = ", ")))
  }
  invisible(TRUE)
}

## ---- House choices ---------------------------------------------------------
# The few things this project decides rather than inherits.

# Panel letters, in the journal's style. "topleft" puts the letter in the
# plot's own margin; a fixed position would land it on the axis of whichever
# panel has the narrowest margin.
TAG_THEME <- theme(plot.tag = element_text(size = 13, face = "bold"),
                   plot.tag.position = "topleft")

# Footnote marker. The default pdf device's font has no dagger, so a glyph
# that survives the round trip is used instead.
FOOTNOTE_MARK <- "*" 

# Population dynamics: spore and non-spore are told apart by grey level alone,
# and the panel is smaller in the merged figure than it is on its own, so the
# two greys are pushed further apart than the standalone version uses.
CELL_GREYS  <- c(Spore = "grey20", `Non-spore` = "grey68")
CELL_SIZE   <- 2.4

# The total (S + V) is counted directly, while spore and non-spore come from
# the latent-state model. Drawing it over them puts the modelled split against
# the quantity that was actually measured; a colour rather than a third grey
# keeps it from reading as a third cell type.
TOTAL_COLOUR <- "#8C1515"

# Figure 2 works on the measured scale rather than relative to the ancestor.
# The biofilm ancestor is provisional -- data/biofil.csv records it as
# B. subtilis 168 delta 6, which is probably not this experiment's ancestor --
# so every relative biofilm value is provisional with it, while the absolute
# posteriors are unaffected. Growth and biofilm are therefore both shown
# absolute, which keeps the four panels on one footing and keeps the figure
# clear of that problem. Set to "relative" to switch back; see 4.Biofilm's
# README for what that would inherit.
PHENOTYPE_SCALE <- "absolute"

# Whether Figure 2 prints the ancestor caveat on the panel itself. The caveat
# is real either way and is always written to output/figure2_caption_notes.md;
# this only controls whether it is set on the figure, which is a decision for
# the authors rather than for the code.
SHOW_ANCESTOR_FOOTNOTE <- FALSE

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 300)
  invisible(file.path(OUT_DIR, stem))
}
