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

# The tail tolerance applies only when the posterior was drawn under the
# original model -- the conjugate Gibbs or JAGS engine, with Gamma(0.01, 0.01)
# on the precision and weights summing to one within each clone.
#
# Under brms it cannot apply, and should not. That Gamma prior is flat on
# log(sigma) only while the precision is small; normalising the weights to sum
# to one scales the precision up about six-fold, into the region where the
# prior's exp(-0.01 * tau) term pulls it down and widens every interval. The
# same Gibbs sampler with the same prior gives 95% intervals about a third
# narrower if the weights are instead scaled to average one -- and a genuinely
# vague prior would not care which. brms uses a prior that is vague on
# log(sigma) (Normal(0, 5)), so its intervals are narrower than the archived
# ones by design, not by error. The medians are what must reproduce, and do.

post <- readRDS(posterior_file)
here <- post$umax$quantiles$relative
CHECK_TAILS <- post$umax$engine %in% c("gibbs", "jags")

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
ok_tails  <- !CHECK_TAILS ||
             max(cmp$lower_pct_diff, cmp$upper_pct_diff) < TOL_TAIL * 100

cat(sprintf("\nengine: %s\nlargest difference: medians %.2f%% (tolerance %.0f%%), tails %.2f%% (%s)\n",
            post$umax$engine,
            max(cmp$median_pct_diff), TOL_MEDIAN * 100,
            max(cmp$lower_pct_diff, cmp$upper_pct_diff),
            if (CHECK_TAILS) sprintf("tolerance %.0f%%", TOL_TAIL * 100)
            else "not checked -- different prior on sigma, see header"))

if (!CHECK_TAILS) {
  here_w <- here[, "97.5%"] - here[, "2.5%"]
  ref_w  <- ref$q97.5 - ref$q2.5
  cat(sprintf("interval width vs archived run: %.0f%% narrower on average, narrower for %d of %d clones\n",
              100 * (1 - mean(here_w / ref_w)), sum(here_w < ref_w), length(here_w)))
}

if (ok_median && ok_tails) {
  cat("PASS -- reproduces the original posterior medians",
      if (CHECK_TAILS) "and tails" else "", "\n")
} else {
  stop("FAIL -- posterior differs from the archived JAGS run by more than tolerance.")
}
