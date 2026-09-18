################################################################################
# 90_latent_states_jags.R
#
# Re-run the latent-state model that separates spores from vegetative cells.
#
# Input : data/spore.transition.txt
# Output: output/latent_states_refit/transition.Bayes.rep{A,B,C,D}.txt
#
# NOT part of run_all.R: it needs JAGS installed, takes a few minutes, and the
# archived fits in data/latent_states/ are what the analysis reads. Refits are
# written alongside the archive rather than over it, and compared at the end.
#
# The model. Total and spore abundances are counted; vegetative cells are not.
# Subtracting gives a negative answer at 174 of 1844 replicate-timepoints, so
# instead both are treated as unobserved states that the counts speak to:
#
#   S[i] ~ Gamma(s_alpha, s_beta)          spore state
#   V[i] ~ Gamma(v_alpha, v_beta)          vegetative state
#   Sobs[i] ~ Normal(S[i],        tau_spore_obs)
#   Tobs[i] ~ Normal(V[i] + S[i], tau_total_obs)
#
# The counted total constrains the sum, the counted spores constrain one part,
# and the gamma priors keep both parts positive. Prior shape and rate come
# from moment-matching the observed abundances -- separately for the
# transition phase (before 240 h) and the equilibrium phase after it, since
# the two differ by orders of magnitude.
#
# Two tidy-ups against the original model text, neither of which changes the
# posterior: sig.v/sig.s and their precisions were declared but never used by
# the likelihood, and have been dropped; and the two observation precisions
# were named the other way round from the quantity they were applied to
# (tau.to on Sobs, tau.so on Tobs), which only mislabelled the output.
#
# Requires: rjags, and a JAGS installation (https://mcmc-jags.sourceforge.io).
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the populationDynamics project root (or R/).")
source(.setup[1])

if (!requireNamespace("rjags", quietly = TRUE)) {
  stop("Package 'rjags' is required, along with a JAGS installation.\n",
       "  macOS:  brew install jags && Rscript -e 'install.packages(\"rjags\")'\n",
       "  or download JAGS from https://mcmc-jags.sourceforge.io")
}

set.seed(20151112)

N_CHAINS <- 5L
N_BURN   <- 1000L
N_ITER   <- 5000L

REFIT_DIR <- file.path(OUT_DIR, "latent_states_refit")
dir.create(REFIT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- data ------------------------------------------------------------------

counts <- read.table(file.path(DATA_DIR, "spore.transition.txt"), header = TRUE)
names(counts)[1] <- "time.h"

# Counts go into the model divided by 10^4, which keeps them small enough not
# to lean on the gamma priors. R/01_abundances.R multiplies back.
counts <- counts %>%
  mutate(time  = time.h / 24,
         total = total * 10^total_dil * PLATE_VOLUME_FACTOR / MODEL_SCALE,
         spore = spore * 10^spor_dil  * PLATE_VOLUME_FACTOR / MODEL_SCALE)

# Gamma shape and rate from the mean and variance of a sample, by moments:
#   mean = a / b, var = a / b^2. Taking logs makes it a 2x2 linear solve.
gamma_moments <- function(x) {
  as.vector(exp(solve(matrix(c(1, 1, -1, -2), 2, 2),
                      matrix(log(c(mean(x), var(x))), 2, 1))))
}

MODEL <- "model{
  for(i in 1:n_transition){
    V[i] ~ dgamma(v_alpha1, v_beta1)
    S[i] ~ dgamma(s_alpha1, s_beta1)
    Sobs[i] ~ dnorm(S[i],        tau_spore_obs)
    Tobs[i] ~ dnorm(V[i] + S[i], tau_total_obs)
  }
  for(i in (n_transition + 1):n){
    V[i] ~ dgamma(v_alpha2, v_beta2)
    S[i] ~ dgamma(s_alpha2, s_beta2)
    Sobs[i] ~ dnorm(S[i],        tau_spore_obs)
    Tobs[i] ~ dnorm(V[i] + S[i], tau_total_obs)
  }
  sig_spore_obs ~ dunif(0, 1000)
  sig_total_obs ~ dunif(0, 1000)
  tau_spore_obs <- pow(sig_spore_obs, -2)
  tau_total_obs <- pow(sig_total_obs, -2)
}"

fit_replicate <- function(rep_id) {
  d <- counts[counts$rep == rep_id, ]
  d <- d[order(d$time.h), ]

  T_obs <- d$total
  S_obs <- d$spore

  transition <- d$time.h <  PHASE_BREAK_H
  equilibrium <- d$time.h >  PHASE_BREAK_H
  n_transition <- sum(transition)

  # The naive vegetative estimate, used only to set the priors, so the
  # non-positive values are dropped here rather than modelled away.
  v_naive <- T_obs - S_obs
  usable  <- v_naive > 0

  jdat <- list(
    n = nrow(d), n_transition = n_transition,
    Sobs = S_obs, Tobs = T_obs,
    s_alpha1 = gamma_moments(S_obs[transition])[1],
    s_beta1  = gamma_moments(S_obs[transition])[2],
    s_alpha2 = gamma_moments(S_obs[equilibrium])[1],
    s_beta2  = gamma_moments(S_obs[equilibrium])[2],
    v_alpha1 = gamma_moments(v_naive[usable & transition])[1],
    v_beta1  = gamma_moments(v_naive[usable & transition])[2],
    v_alpha2 = gamma_moments(v_naive[usable & equilibrium])[1],
    v_beta2  = gamma_moments(v_naive[usable & equilibrium])[2]
  )

  message(sprintf("  rep %s: %d timepoints (%d transition, %d equilibrium)",
                  rep_id, nrow(d), n_transition, sum(equilibrium)))

  mod <- rjags::jags.model(textConnection(MODEL), data = jdat,
                           n.chains = N_CHAINS, quiet = TRUE)
  update(mod, N_BURN)
  samples <- rjags::coda.samples(
    mod, c("V", "S", "sig_spore_obs", "sig_total_obs"), n.iter = N_ITER)

  f <- summary(samples)
  pick <- function(tab, prefix, col) {
    rows <- grep(sprintf("^%s\\[", prefix), rownames(tab))
    tab[rows, col][order(as.integer(gsub("\\D", "", rownames(tab)[rows])))]
  }

  data.frame(
    time    = d$time,
    s.est   = pick(f$statistics, "S", "Mean"),
    s.lower = pick(f$quantiles,  "S", "2.5%"),
    s.upper = pick(f$quantiles,  "S", "97.5%"),
    v.est   = pick(f$statistics, "V", "Mean"),
    v.lower = pick(f$quantiles,  "V", "2.5%"),
    v.upper = pick(f$quantiles,  "V", "97.5%")
  )
}

## ---- run -------------------------------------------------------------------
# The original ran one replicate at a time, by editing `trans.rep <- trans.A`
# and the output filename by hand. That is how replicate D ended up reading
# replicate C's file downstream.

message("Fitting the latent-state model for each replicate")

for (rep_id in REPS) {
  out <- fit_replicate(rep_id)
  write.table(out,
              file.path(REFIT_DIR, sprintf("transition.Bayes.rep%s.txt", rep_id)),
              sep = "\t", col.names = TRUE, row.names = FALSE)
}

## ---- compare against the archived fits -------------------------------------

cat("\nRefit vs. archived state estimates (median absolute % difference)\n\n")
for (rep_id in REPS) {
  a <- read.table(file.path(DATA_DIR, "latent_states",
                            sprintf("transition.Bayes.rep%s.txt", rep_id)), header = TRUE)
  b <- read.table(file.path(REFIT_DIR,
                            sprintf("transition.Bayes.rep%s.txt", rep_id)), header = TRUE)
  cat(sprintf("  rep %s: spore %.2f%%, vegetative %.2f%%\n", rep_id,
              100 * median(abs(b$s.est - a$s.est) / a$s.est),
              100 * median(abs(b$v.est - a$v.est) / a$v.est)))
}

message("\nRefits written to output/latent_states_refit/ (the archive is untouched).")
