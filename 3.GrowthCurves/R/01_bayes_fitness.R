################################################################################
# 01_bayes_fitness.R
#
# Estimate a posterior distribution for each Gompertz growth parameter of each
# clone, then express it relative to the ancestor.
#
# Input : data/comp_data_annotated.csv -- one row per replicate growth curve,
#                                 holding the modified-Gompertz fit (A, umax, L),
#                                 the standard error of each fitted parameter,
#                                 and the plate run the curve came from.
#                                 Built by R/91_build_comp_data.R.
# Output: output/posteriors.rds   -- posterior draws + summaries for umax, A, L
#         output/posterior_summary.csv
#
# Model (independently for each parameter and each clone j):
#
#     x[i, j]  ~  LogNormal( mu[j],  tau[j] * w[i, j] )      i = 1..6 replicates
#     mu[j]    ~  Normal(0, precision = 0.001)
#     tau[j]   ~  Gamma(0.01, 0.01)
#
# Each replicate carries its own uncertainty from the curve fit, so replicates
# are weighted by inverse variance (w[i,j], normalised to sum to 1 within a
# clone) rather than averaged naively. The log-normal keeps the estimates
# positive during sampling. exp(mu[j]) is the clone's growth parameter, and
# exp(mu[j]) / exp(mu[ancestor]) is its relative fitness.
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
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

set.seed(20230826)   # the date of the last growth-curve run

N_CHAINS <- 4L
N_BURN   <- 2000L    # 1000 adapt + 1000 update in the original
N_ITER   <- 10000L   # draws kept per chain

## ---- 1. read the curve fits and shape them into 6 x 11 matrices ------------

read_fits <- function(path = file.path(DATA_DIR, "comp_data_annotated.csv")) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$clones %in% MODEL_ORDER, ]
  d$clones <- factor(d$clones, levels = MODEL_ORDER)

  n <- table(d$clones)
  if (any(n != N_REPS)) {
    stop("Expected ", N_REPS, " replicates per clone, got: ",
         paste(names(n), n, sep = "=", collapse = ", "))
  }
  d[order(d$clones), ]
}

# Pull one fitted parameter (and its SE) out as a replicate-by-clone matrix.
as_matrix <- function(fits, column) {
  m <- matrix(fits[[column]], nrow = N_REPS, ncol = length(MODEL_ORDER),
              dimnames = list(NULL, MODEL_ORDER))
  m
}

# Inverse-variance weights, normalised to sum to 1 within each clone.
inv_var_weights <- function(se_mat) {
  prec <- 1 / se_mat^2
  sweep(prec, MARGIN = 2, STATS = colSums(prec), FUN = "/")
}

## ---- 2. samplers -----------------------------------------------------------

# Conjugate Gibbs sampler for a single clone. Because the priors are
# independent across clones, the eleven groups factorise and can be sampled
# one at a time.
#
#   y = log(x), so   y[i] ~ Normal(mu, precision = tau * w[i])
#   mu  | tau, y  ~  Normal( tau * sum(w*y) / P,  precision P = 0.001 + tau*sum(w) )
#   tau | mu,  y  ~  Gamma( 0.01 + n/2,  rate = 0.01 + sum(w*(y-mu)^2)/2 )
gibbs_one <- function(x, w, n_iter, n_burn, mu_prior_prec = 0.001,
                      tau_prior_shape = 0.01, tau_prior_rate = 0.01) {
  y  <- log(x)
  n  <- length(y)
  sw <- sum(w)

  tau  <- 1 / var(y)                      # start at the unweighted precision
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
  clones <- colnames(x_mat)
  draws <- lapply(seq_len(n_chains), function(chain) {
    vapply(clones,
           function(cl) gibbs_one(x_mat[, cl], w_mat[, cl], n_iter, n_burn),
           numeric(n_iter))
  })
  do.call(rbind, draws)               # (n_chains * n_iter) x 11, on the log scale
}

sample_jags <- function(x_mat, w_mat, n_chains, n_iter, n_burn) {
  mod <- "model{
       for(i in 1:6){
         for(j in 1:10){
           M_mu[i,j] ~ dlnorm(mu_m[j], tau_m[j] * wts_m[i,j])
         }
         An_mu[i] ~ dlnorm(mu_a, tau_a * wts_a[i])
       }
       for(i in 1:10){
         mu_m[i]  ~ dnorm(0, .001)
         tau_m[i] ~ dgamma(.01, .01)
       }
       mu_a  ~ dnorm(0, .001)
       tau_a ~ dgamma(.01, .01)
  }"

  jdat <- list(M_mu  = x_mat[, -1, drop = FALSE], An_mu = x_mat[, 1],
               wts_m = w_mat[, -1, drop = FALSE], wts_a = w_mat[, 1])

  jmod <- rjags::jags.model(textConnection(mod), n.chains = n_chains,
                            data = jdat, n.adapt = n_burn / 2)
  update(jmod, n_burn / 2)
  ans <- rjags::coda.samples(jmod, c("mu_m", "mu_a"), n_iter)

  # coda orders the monitored nodes mu_a, mu_m[1] ... mu_m[10], which is
  # MODEL_ORDER.
  vals <- do.call(rbind, lapply(ans, as.matrix))
  colnames(vals) <- MODEL_ORDER
  vals
}

## ---- 3. fit one parameter --------------------------------------------------

#' @param fits     data frame from read_fits()
#' @param value    column holding the fitted parameter (e.g. "umax")
#' @param se       column holding its standard error (e.g. "umax.se")
#' @param relative TRUE to also return draws divided by the ancestor's
fit_parameter <- function(fits, value, se, engine = c("auto", "gibbs", "jags")) {
  engine <- match.arg(engine)
  if (engine == "auto") {
    engine <- if (requireNamespace("rjags", quietly = TRUE)) "jags" else "gibbs"
  }

  x_mat <- as_matrix(fits, value)
  w_mat <- inv_var_weights(as_matrix(fits, se))

  message("  ", value, ": sampling with ", engine, " (",
          N_CHAINS, " chains x ", format(N_ITER, big.mark = ","), " draws)")

  log_draws <- switch(engine,
    gibbs = sample_gibbs(x_mat, w_mat, N_CHAINS, N_ITER, N_BURN),
    jags  = sample_jags(x_mat, w_mat, N_CHAINS, N_ITER, N_BURN))

  absolute <- exp(log_draws)                       # clone-level parameter
  relative <- absolute[, -1] / absolute[, 1]       # relative to the ancestor

  # Reference point the original code printed alongside the posterior.
  weighted_mean <- colSums(x_mat * w_mat)

  list(
    engine        = engine,
    absolute      = absolute,
    relative      = relative,
    weighted_mean = weighted_mean,
    quantiles     = list(
      absolute = t(apply(absolute, 2, quantile, c(.025, .5, .975))),
      relative = t(apply(relative, 2, quantile, c(.025, .5, .975)))
    )
  )
}

## ---- 4. run ---------------------------------------------------------------

local({

  fits <- read_fits()

  message("Fitting Bayesian models (", nrow(fits), " curves, ",
          length(MODEL_ORDER), " strains)")

  PARAMS <- list(
    umax = c(value = "umax", se = "umax.se"),   # maximum growth rate
    A    = c(value = "A",    se = "A.se"),      # maximum yield (carrying capacity)
    L    = c(value = "L",    se = "L.se")       # lag time
  )

  post <- lapply(PARAMS, function(p) fit_parameter(fits, p[["value"]], p[["se"]]))

  # Diagnostic the original code also ran. The posterior median of a log-normal
  # is not the arithmetic weighted mean, so these agree closely only when the
  # replicates are tight (umax, A) and drift apart when they are not (L).
  message("\nPosterior median vs. inverse-variance weighted mean:")
  for (nm in names(post)) {
    med <- post[[nm]]$quantiles$absolute[, "50%"]
    message(sprintf("  %-4s max relative difference: %.3f%%", nm,
                    100 * max(abs(med - post[[nm]]$weighted_mean) /
                                post[[nm]]$weighted_mean)))
  }

  saveRDS(post, file.path(OUT_DIR, "posteriors.rds"))

  summary_df <- bind_rows(lapply(names(post), function(nm) {
    q_abs <- post[[nm]]$quantiles$absolute
    q_rel <- post[[nm]]$quantiles$relative
    bind_rows(
      data.frame(parameter = nm, scale = "absolute", clone = rownames(q_abs),
                 q_abs, weighted_mean = post[[nm]]$weighted_mean,
                 check.names = FALSE, row.names = NULL),
      data.frame(parameter = nm, scale = "relative", clone = rownames(q_rel),
                 q_rel, weighted_mean = NA_real_,
                 check.names = FALSE, row.names = NULL)
    )
  }))

  write.csv(summary_df, file.path(OUT_DIR, "posterior_summary.csv"),
            row.names = FALSE)

  message("\nWrote output/posteriors.rds and output/posterior_summary.csv")
})
