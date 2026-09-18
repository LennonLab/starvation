################################################################################
# 03_validate.R
#
# Check this recreation against the original JAGS run.
#
# The archived reference values are the umax posterior quantiles relative to the
# ancestor, recovered from the SporeMut repo (see the header of
# data/reference_jags_umax_relative.csv). Agreement to within Monte Carlo noise
# means the rewritten sampler is estimating the same posterior as the original.
#
# Input : output/posteriors.rds, data/reference_jags_umax_relative.csv
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

posterior_file <- file.path(OUT_DIR, "posteriors.rds")
if (!file.exists(posterior_file)) {
  stop("output/posteriors.rds not found -- run R/01_bayes_fitness.R first.")
}

# Two independent MCMC runs of the same model will not agree exactly. Medians
# are well determined; the 2.5% and 97.5% tails are noisier.
TOL_MEDIAN <- 0.01   # 1%
TOL_TAIL   <- 0.03   # 3%

post <- readRDS(posterior_file)
here <- post$umax$quantiles$relative

ref <- read.csv(file.path(DATA_DIR, "reference_jags_umax_relative.csv"),
                comment.char = "#")
stopifnot(setequal(ref$clone, rownames(here)))
ref <- ref[match(rownames(here), ref$clone), ]

cmp <- data.frame(
  clone = ref$clone,
  ref_median  = ref$q50,   this_median  = here[, "50%"],
  ref_lower   = ref$q2.5,  this_lower   = here[, "2.5%"],
  ref_upper   = ref$q97.5, this_upper   = here[, "97.5%"],
  row.names = NULL
)
cmp$median_pct_diff <- 100 * abs(cmp$this_median - cmp$ref_median) / cmp$ref_median
cmp$lower_pct_diff  <- 100 * abs(cmp$this_lower  - cmp$ref_lower)  / cmp$ref_lower
cmp$upper_pct_diff  <- 100 * abs(cmp$this_upper  - cmp$ref_upper)  / cmp$ref_upper

cat("umax relative to ancestor -- this run vs. the original JAGS run\n\n")
print(cmp[, c("clone", "ref_median", "this_median", "median_pct_diff",
              "lower_pct_diff", "upper_pct_diff")],
      digits = 3, row.names = FALSE)

ok_median <- max(cmp$median_pct_diff) < TOL_MEDIAN * 100
ok_tails  <- max(cmp$lower_pct_diff, cmp$upper_pct_diff) < TOL_TAIL * 100

cat(sprintf("\nlargest difference: medians %.2f%% (tolerance %.0f%%), tails %.2f%% (tolerance %.0f%%)\n",
            max(cmp$median_pct_diff), TOL_MEDIAN * 100,
            max(cmp$lower_pct_diff, cmp$upper_pct_diff), TOL_TAIL * 100))

if (ok_median && ok_tails) {
  cat("PASS -- reproduces the original posterior within Monte Carlo noise.\n")
} else {
  stop("FAIL -- posterior differs from the archived JAGS run by more than tolerance.")
}
