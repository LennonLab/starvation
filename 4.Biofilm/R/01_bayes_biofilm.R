################################################################################
# 01_bayes_biofilm.R
#
# Estimate a posterior distribution for each strain's biofilm production, then
# express it relative to the ancestor.
#
# Input : data/biofil.csv  -- one row per well: raw OD550, the plate blank, and
#                             the blank-corrected reading used here, plus the
#                             strain's treatment and mutation.
# Output: output/posteriors.rds, output/posterior_summary.csv
#
# Model, independently for each strain j:
#
#     x[i, j]  ~  LogNormal( mu[j], tau[j] )        i = 1..8 wells
#     mu[j]    ~  Normal(0, precision = 0.001)
#     tau[j]   ~  Gamma(0.01, 0.01)
#
# Unlike the growth-curve assay there is no per-well standard error to weight
# by -- each well is a single OD reading -- so the replicates enter unweighted
# and each strain gets its own variance. The log-normal keeps estimates
# positive during sampling and suits a response spanning two orders of
# magnitude. exp(mu[j]) is the strain's biofilm level and
# exp(mu[j]) / exp(mu[ancestor]) is its level relative to the ancestor.
#
# Two samplers give the same posterior:
#   engine = "jags"  -- the original rjags implementation (needs JAGS installed)
#   engine = "gibbs" -- the exact conjugate Gibbs sampler for the model above,
#                       written in base R, so the figures reproduce with no
#                       external dependencies
#   engine = "auto"  -- rjags if available, otherwise gibbs (default)
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the biofilm project root (or R/).")
source(.setup[1])

set.seed(20230826)

N_CHAINS <- 4L
N_BURN   <- 2000L    # 1000 adapt + 1000 update in the original
N_ITER   <- 10000L   # draws kept per chain

RESPONSE <- "OD550_Corrected"

## ---- 1. read the assay and shape it into an 8 x 11 matrix ------------------

read_assay <- function(path = file.path(DATA_DIR, "biofil.csv")) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$clones <- factor(d$clones, levels = MODEL_ORDER)
  if (anyNA(d$clones)) stop("Unexpected strain in ", basename(path))

  n <- table(d$clones)
  if (any(n != N_REPS)) {
    stop("Expected ", N_REPS, " wells per strain, got: ",
         paste(names(n), n, sep = "=", collapse = ", "))
  }
  d[order(d$clones), ]
}

as_matrix <- function(assay, column) {
  matrix(assay[[column]], nrow = N_REPS, ncol = length(MODEL_ORDER),
         dimnames = list(NULL, MODEL_ORDER))
}

## ---- 2. samplers -----------------------------------------------------------

# Conjugate Gibbs sampler for a single strain. Because the priors are
# independent across strains, the eleven groups factorise and can be sampled
# one at a time.
#
#   y = log(x), so   y[i] ~ Normal(mu, precision = tau * w[i])
#   mu  | tau, y  ~  Normal( tau * sum(w*y) / P,  precision P = 0.001 + tau*sum(w) )
#   tau | mu,  y  ~  Gamma( 0.01 + n/2,  rate = 0.01 + sum(w*(y-mu)^2)/2 )
#
# `w` is all ones for this assay; the argument is kept so the sampler matches
# the weighted version used for the growth curves.
gibbs_one <- function(x, w, n_iter, n_burn, mu_prior_prec = 0.001,
                      tau_prior_shape = 0.01, tau_prior_rate = 0.01) {
  y  <- log(x)
  n  <- length(y)
  sw <- sum(w)

  tau  <- 1 / var(y)
  keep <- numeric(n_iter)

  for (it in seq_len(n_iter + n_burn)) {
    prec <- mu_prior_prec + tau * sw
    mu   <- rnorm(1, mean = tau * sum(w * y) / prec, sd = 1 / sqrt(prec))
    tau  <- rgamma(1, shape = tau_prior_shape + n / 2,
                   rate  = tau_prior_rate + sum(w * (y - mu)^2) / 2)
    if (it > n_burn) keep[it - n_burn] <- mu
  }
  keep
}

sample_gibbs <- function(x_mat, w_mat, n_chains, n_iter, n_burn) {
  strains <- colnames(x_mat)
  draws <- lapply(seq_len(n_chains), function(chain) {
    vapply(strains,
           function(s) gibbs_one(x_mat[, s], w_mat[, s], n_iter, n_burn),
           numeric(n_iter))
  })
  do.call(rbind, draws)               # (n_chains * n_iter) x 11, on the log scale
}

sample_jags <- function(x_mat, w_mat, n_chains, n_iter, n_burn) {
  mod <- "model{
       for(i in 1:8){
         for(j in 1:10){
           M_mu[i,j] ~ dlnorm(mu_m[j], tau_m[j])
         }
         An_mu[i] ~ dlnorm(mu_a, tau_a)
       }
       for(i in 1:10){
         mu_m[i]  ~ dnorm(0, .001)
         tau_m[i] ~ dgamma(.01, .01)
       }
       mu_a  ~ dnorm(0, .001)
       tau_a ~ dgamma(.01, .01)
  }"

  # The ancestor is the last column of MODEL_ORDER for this assay.
  anc <- match("ancestor", colnames(x_mat))
  jdat <- list(M_mu = x_mat[, -anc, drop = FALSE], An_mu = x_mat[, anc])

  jmod <- rjags::jags.model(textConnection(mod), n.chains = n_chains,
                            data = jdat, n.adapt = n_burn / 2)
  update(jmod, n_burn / 2)
  ans <- rjags::coda.samples(jmod, c("mu_m", "mu_a"), n_iter)

  # coda orders the monitored nodes mu_a, mu_m[1] ... mu_m[10].
  vals <- do.call(rbind, lapply(ans, as.matrix))
  colnames(vals) <- c("ancestor", setdiff(MODEL_ORDER, "ancestor"))
  vals[, MODEL_ORDER]
}

## ---- 3. fit ----------------------------------------------------------------

fit_response <- function(assay, column, engine = c("auto", "gibbs", "jags")) {
  engine <- match.arg(engine)
  if (engine == "auto") {
    engine <- if (requireNamespace("rjags", quietly = TRUE)) "jags" else "gibbs"
  }

  x_mat <- as_matrix(assay, column)
  w_mat <- matrix(1, nrow(x_mat), ncol(x_mat), dimnames = dimnames(x_mat))

  message("  ", column, ": sampling with ", engine, " (",
          N_CHAINS, " chains x ", format(N_ITER, big.mark = ","), " draws)")

  log_draws <- switch(engine,
    gibbs = sample_gibbs(x_mat, w_mat, N_CHAINS, N_ITER, N_BURN),
    jags  = sample_jags(x_mat, w_mat, N_CHAINS, N_ITER, N_BURN))

  absolute <- exp(log_draws)
  # Provisional: see warn_ancestor() in 00_setup.R -- the denominator strain
  # may not be this experiment's ancestor.
  relative <- absolute[, setdiff(MODEL_ORDER, "ancestor")] / absolute[, "ancestor"]

  # With a vague prior the posterior median of an unweighted log-normal should
  # sit on the geometric mean of the replicates.
  geo_mean <- exp(colMeans(log(x_mat)))

  list(
    engine    = engine,
    absolute  = absolute,
    relative  = relative,
    geo_mean  = geo_mean,
    quantiles = list(
      absolute = t(apply(absolute, 2, quantile, c(.025, .5, .975))),
      relative = t(apply(relative, 2, quantile, c(.025, .5, .975)))
    )
  )
}

## ---- 4. run ----------------------------------------------------------------

local({
  assay <- read_assay()

  message("Fitting the biofilm model (", nrow(assay), " wells, ",
          length(MODEL_ORDER), " strains)")

  post <- list(biofilm = fit_response(assay, RESPONSE))

  med <- post$biofilm$quantiles$absolute[, "50%"]
  message(sprintf(
    "\nPosterior median vs. geometric mean of the wells: max difference %.2f%%",
    100 * max(abs(med - post$biofilm$geo_mean) / post$biofilm$geo_mean)))

  saveRDS(post, file.path(OUT_DIR, "posteriors.rds"))

  summary_df <- bind_rows(
    data.frame(parameter = "biofilm", scale = "absolute",
               clone = rownames(post$biofilm$quantiles$absolute),
               post$biofilm$quantiles$absolute,
               geo_mean = post$biofilm$geo_mean,
               check.names = FALSE, row.names = NULL),
    data.frame(parameter = "biofilm", scale = "relative",
               clone = rownames(post$biofilm$quantiles$relative),
               post$biofilm$quantiles$relative,
               geo_mean = NA_real_,
               check.names = FALSE, row.names = NULL))

  write.csv(summary_df, file.path(OUT_DIR, "posterior_summary.csv"),
            row.names = FALSE)

  message("Wrote output/posteriors.rds and output/posterior_summary.csv")
})

warn_ancestor()
