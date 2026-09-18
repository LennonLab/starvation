################################################################################
# 91_build_comp_data.R
#
# Rebuild the analysis input from the archived modified-Gompertz fits, and label
# every curve with the plate run it came from.
#
# Input : data/gompertz_fits/*.fit.parms.txt   (one file per strain per run)
#         data/plate_runs.csv                  (which run each fit came from)
#         data/treatments_original_corrected.csv
# Output: data/comp_data_annotated.csv  -- the analysis input, with `run`,
#                                          `well` and `curve_id` added and the
#                                          S6 mix-up corrected
#         output/comp_data_rebuilt.csv  -- byte-for-byte rebuild of the
#                                          archived data/comp_data.csv
#
# Why the extra columns: the `Curve` names are raw well labels (`M21_1` =
# strain M21, well 1), and each strain was run on two plates whose wells are
# both numbered from 1, so a name like `M21_1` appears twice in the archived
# file with no way to tell the runs apart. `run` and `curve_id` fix that.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

FITS_DIR <- file.path(DATA_DIR, "gompertz_fits")
runs_tbl <- read.csv(file.path(DATA_DIR, "plate_runs.csv"), comment.char = "#")

# One archived fit file, tagged with the run it came from.
fit <- function(name) {
  d <- read.table(file.path(FITS_DIR, paste0(name, ".fit.parms.txt")),
                  sep = ",", header = TRUE)
  d$run <- runs_tbl$run[match(name, runs_tbl$fit_name)]
  d$fit_name <- name
  d
}

## ---- 1. combine runs, dropping curves rejected during fit diagnostics ------
# The negative indices are the curves that failed visual QC in the original
# analysis (see the diagnostic PDFs in data/gompertz_fits/). They are row
# positions within each archived file, so they only make sense against these
# exact files -- refitting from raw plates renumbers them.
#
# fix_s6: the original read S1_new where it meant S6_new, so the S6 rows in the
# archived comp_data.csv are a second copy of S1. FALSE reproduces the archive;
# TRUE reads the S6 fits that were there all along.
build_strains <- function(fix_s6) list(
  ancestor = fit("anc_new"),
  M4       = rbind(fit("M4_new"),               fit("M4_repeat.new")),
  M13      = rbind(fit("M13_new")[-c(1, 6), ],  fit("M13_repeat.new")),
  M17      = rbind(fit("M17_new")[-1, ],        fit("M17_repeat.new")[-c(2, 5, 6), ]),
  M19      = rbind(fit("M19_new"),              fit("M19_repeat.new")[-c(2, 10), ]),
  M21      = rbind(fit("M21_new")[-c(2, 6), ],  fit("M21_repeat.new")),
  M23      = rbind(fit("M23_repeat.new"),       fit("M23_repeat2.new")),
  M26      = rbind(fit("M26_new"),              fit("M26_repeat.new")[-1, ]),
  M41      = fit("M41_new"),
  M54      = rbind(fit("M54_new")[-3, ],        fit("M54_repeat.new")[-c(2, 4), ]),
  M79      = fit("M79_new"),
  S1       = rbind(fit("S1_new")[-5, ],         fit("S1_repeat.new")),
  S6       = if (fix_s6) fit("S6_new") else fit("S1_new"),
  S11      = rbind(fit("S11_new"),              fit("S11_repeat.new")),
  S22      = rbind(fit("S22_new")[-c(3, 4), ],  fit("S22_repeat.new")),
  S51      = rbind(fit("S51_new"),              fit("S51_repeat.new")),
  S95      = rbind(fit("S95_new"),              fit("S95_repeat.new"))
)

## ---- 2. keep the six best-fitting curves per strain ------------------------
best_six <- function(x) x[x$RSME %in% sort(x$RSME, partial = 1:6)[1:6], ]

build <- function(fix_s6) {
  strains <- build_strains(fix_s6)
  d <- bind_rows(lapply(strains, best_six), .id = "strain")

  # `column_label` is the strain's position in the list, kept for compatibility
  # with the archived file.
  d$column_label <- match(d$strain, names(strains))
  d$clones <- sub("_.*", "", d$Curve)
  d$clones[d$strain == "ancestor"] <- "ancestor"
  d$strain <- NULL

  treats <- read.csv(file.path(DATA_DIR, "treatments_original_corrected.csv"))
  left_join(d, treats, by = "clones")
}

## ---- 3. reproduce the archived file (provenance check) ---------------------
ARCHIVED_COLS <- names(read.csv(file.path(DATA_DIR, "comp_data.csv"), nrows = 1))

as_archived <- function(d) {
  d <- d[, ARCHIVED_COLS]
  rownames(d) <- NULL
  d
}

rebuilt <- as_archived(build(fix_s6 = FALSE))
write.csv(rebuilt, file.path(OUT_DIR, "comp_data_rebuilt.csv"), row.names = FALSE)

shipped <- read.csv(file.path(DATA_DIR, "comp_data.csv"), stringsAsFactors = FALSE)
check   <- all.equal(shipped, read.csv(file.path(OUT_DIR, "comp_data_rebuilt.csv"),
                                       stringsAsFactors = FALSE))
if (isTRUE(check)) {
  message("comp_data.csv reproduces exactly from the archived Gompertz fits.")
} else {
  message("Rebuilt comp_data.csv differs from the archived file:")
  print(check)
}

## ---- 4. write the annotated analysis input ---------------------------------
annotated <- build(fix_s6 = TRUE)

# `well` is the number the plate reader gave the well. R renames a duplicated
# column header by appending ".1", so a well like S95_4.1 is the second column
# labelled S95_4 in data/raw/20230824_...csv -- a plate-map typo, most likely
# meant to be S95_3.
annotated$well <- ifelse(annotated$clones == "ancestor",
                         sub("^ancestor", "", annotated$Curve),
                         sub("^[^_]*_", "", annotated$Curve))
annotated$curve_id <- paste(annotated$clones, annotated$run, annotated$well,
                            sep = "-")

stopifnot(!anyDuplicated(annotated$curve_id))

annotated <- annotated[, c("curve_id", "clones", "run", "well", "Curve",
                           "fit_name", "b0", "b0.se", "A", "A.se",
                           "umax", "umax.se", "L", "L.se", "RSME", "CV",
                           "outlier", "order", "evo.type", "cell.type",
                           "mutation")]

write.csv(annotated, file.path(DATA_DIR, "comp_data_annotated.csv"),
          row.names = FALSE)

changed <- sum(!(annotated$Curve[annotated$clones == "S6"] %in%
                   shipped$Curve[shipped$column_label == 13]))
message("Wrote data/comp_data_annotated.csv (", nrow(annotated), " curves; ",
        "S6 corrected, ", changed, " rows now differ from the archive).")
