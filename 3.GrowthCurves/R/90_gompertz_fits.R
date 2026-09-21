################################################################################
# 90_gompertz_fits.R
#
# Fit the modified Gompertz model to the raw plate-reader curves.
#
# Input : data/raw/*.csv        -- plate-reader exports, one file per run
#                                (the exports do not name the instrument; the
#                                 group's methods document says Epoch2)
#         data/plate_runs.csv   -- which strain came from which run
# Output: output/gompertz_refit/<name>.fit.parms.txt  -- fitted b0, A, umax, L
#         output/gompertz_refit/<name>.fit.parms.pdf  -- per-curve diagnostics
#
# This is the slow step (~2 min) and it is NOT part of run_all.R. The archived
# fits in data/gompertz_fits/ are what the analysis uses, and the curves that
# were rejected by eye are recorded in R/91_build_comp_data.R as row positions
# in those files. Refitting therefore writes to output/gompertz_refit/ rather
# than over the archive; compare the two before adopting a new fit.
#
# Requires: bbmle, MuMIn, nlme, gtools
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

missing <- Filter(function(p) !requireNamespace(p, quietly = TRUE),
                  c("bbmle", "MuMIn", "nlme", "gtools"))
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
}

REFIT_DIR <- file.path(OUT_DIR, "gompertz_refit")
dir.create(REFIT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- raw plate runs --------------------------------------------------------
# Time is exported as "HH:MM:SS"; the original analysis reads it as HH.MM, so
# lag and growth rate are in those units rather than decimal hours. Kept as is
# so the fits match the archived ones.
read_plate <- function(file) {
  d <- read.csv(file.path(DATA_DIR, "raw", file))
  d$Time <- as.numeric(sub("^(\\d+):(\\d+).*", "\\1.\\2", d$Time))
  d
}

runs   <- read.csv(file.path(DATA_DIR, "plate_runs.csv"), comment.char = "#")
plates <- lapply(setNames(nm = unique(runs$raw_file)), read_plate)

## ---- fit -------------------------------------------------------------------
# growth.modGomp() resolves its helpers as "../bin/" and writes to a relative
# output.dir, so run it from inside R/.
old_wd <- setwd(file.path(PROJ, "R"))
on.exit(setwd(old_wd), add = TRUE)

source(file.path(PROJ, "bin", "modified_Gomp_diagnostic3.R"))

for (i in seq_len(nrow(runs))) {
  r     <- runs[i, ]
  wells <- plates[[r$raw_file]] |> dplyr::select(Time, starts_with(r$prefix))
  message(sprintf("[%2d/%d] %s %s -> %s (%d curves)",
                  i, nrow(runs), r$run, r$prefix, r$fit_name, ncol(wells) - 1))

  growth.modGomp(input = wells, output.name = paste0(r$fit_name, ".fit.parms"),
                 output.dir = "../output/gompertz_refit/",
                 synergy = FALSE, temp = FALSE, smooth = TRUE, trim = TRUE)
}

## ---- compare against the archived fits -------------------------------------
setwd(old_wd)

compare_one <- function(name) {
  f_new <- file.path(REFIT_DIR, paste0(name, ".fit.parms.txt"))
  f_old <- file.path(DATA_DIR, "gompertz_fits", paste0(name, ".fit.parms.txt"))
  if (!file.exists(f_new) || !file.exists(f_old)) return(NULL)
  a <- read.table(f_old, sep = ",", header = TRUE)
  b <- read.table(f_new, sep = ",", header = TRUE)
  data.frame(fit = name,
             archived = nrow(a), refit = nrow(b),
             identical = isTRUE(all.equal(a, b)))
}

cmp <- do.call(rbind, lapply(runs$fit_name, compare_one))
cat("\nRefit vs. archived fits:\n")
print(cmp, row.names = FALSE)
cat(sprintf("\n%d of %d fit files reproduce exactly.\n",
            sum(cmp$identical), nrow(cmp)))
message("\nRefits written to output/gompertz_refit/ (the archive is untouched).")
