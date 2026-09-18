################################################################################
# 03_validate.R
#
# Check the Gibbs sampler against an exact calculation.
#
# For one strain the model has only two parameters, and tau can be integrated
# out in closed form:
#
#   p(mu, tau | y) ~ tau^(n/2 + a0 - 1) exp(-tau (b0 + S(mu)/2)) exp(-p0 mu^2/2)
#   integral over tau  =>  p(mu | y) ~ exp(-p0 mu^2 / 2) (b0 + S(mu)/2)^-(n/2+a0)
#
# with S(mu) = sum((log(x) - mu)^2), p0 = 0.001, a0 = b0 = 0.01. That leaves a
# one-dimensional density that can be normalised on a grid and read off to any
# accuracy, so the sampler is checked against the exact posterior rather than
# against another sampler.
#
# Input : output/posteriors.rds, data/biofil.csv
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the biofilm project root (or R/).")
source(.setup[1])

posterior_file <- file.path(OUT_DIR, "posteriors.rds")
if (!file.exists(posterior_file)) {
  stop("output/posteriors.rds not found -- run R/01_bayes_biofilm.R first.")
}
post  <- readRDS(posterior_file)$biofilm
assay <- read.csv(file.path(DATA_DIR, "biofil.csv"), stringsAsFactors = FALSE)

TOL_MEDIAN <- 0.005   # 0.5% on exp(mu)
TOL_TAIL   <- 0.02    # 2% on the interval bounds

# Exact marginal posterior of mu for one strain, evaluated on a grid.
exact_quantiles <- function(x, probs = c(.025, .5, .975),
                            p0 = 0.001, a0 = 0.01, b0 = 0.01, n_grid = 200001) {
  y <- log(x)
  n <- length(y)

  # A generous window around the sample mean; the density is negligible well
  # inside it, which the tail mass check below confirms.
  centre <- mean(y)
  half   <- 25 * sd(y)
  mu     <- seq(centre - half, centre + half, length.out = n_grid)

  S <- vapply(mu, function(m) sum((y - m)^2), numeric(1))
  log_dens <- -p0 * mu^2 / 2 - (n / 2 + a0) * log(b0 + S / 2)
  dens <- exp(log_dens - max(log_dens))

  cdf <- cumsum(dens) / sum(dens)
  if (cdf[10] > 1e-8 || cdf[n_grid - 10] < 1 - 1e-8) {
    stop("Grid too narrow for the exact posterior.")
  }
  stats::approx(cdf, mu, xout = probs)$y
}

strains <- colnames(post$absolute)
rows <- lapply(strains, function(s) {
  x <- assay$OD550_Corrected[assay$clones == s]
  ex <- exp(exact_quantiles(x))
  gi <- post$quantiles$absolute[s, ]
  data.frame(
    strain = s,
    exact_median = ex[2], gibbs_median = gi[["50%"]],
    median_pct = 100 * abs(gi[["50%"]] - ex[2]) / ex[2],
    lower_pct  = 100 * abs(gi[["2.5%"]]  - ex[1]) / ex[1],
    upper_pct  = 100 * abs(gi[["97.5%"]] - ex[3]) / ex[3]
  )
})
cmp <- do.call(rbind, rows)

cat("Gibbs sampler vs. exact marginal posterior of exp(mu)\n\n")
print(cmp, digits = 3, row.names = FALSE)

ok_median <- max(cmp$median_pct) < TOL_MEDIAN * 100
ok_tails  <- max(cmp$lower_pct, cmp$upper_pct) < TOL_TAIL * 100

cat(sprintf(
  "\nlargest difference: medians %.3f%% (tolerance %.1f%%), tails %.3f%% (tolerance %.1f%%)\n",
  max(cmp$median_pct), TOL_MEDIAN * 100,
  max(cmp$lower_pct, cmp$upper_pct), TOL_TAIL * 100))

if (ok_median && ok_tails) {
  cat("PASS -- the sampler reproduces the exact posterior within Monte Carlo noise.\n")
} else {
  stop("FAIL -- the sampler differs from the exact posterior by more than tolerance.")
}
