################################################################################
# 91_latent_state_diagnostics.R
#
# Convergence diagnostics for the latent-state model, and a check that it
# reproduces the archived fits.
#
# The original chains were not kept, and data/latent_states/*.txt hold only
# posterior summaries, so R-hat and effective sample size cannot be recovered
# from them -- the model has to be rerun. This does that.
#
# It drives JAGS through its command-line interface rather than through rjags.
# rjags 4-17 requires JAGS 4.x and refuses to build against the JAGS 5 that
# homebrew installs; the CLI writes CODA files that the coda package reads, so
# nothing is lost by going this way.
#
# Input : data/spore.transition.txt
# Output: output/latent_state_diagnostics.csv
#         output/table_latent_diagnostics.md
#
# Requires: a JAGS binary on PATH. macOS: brew install jags.
################################################################################

.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the populationDynamics project root (or R/).")
source(.setup[1])
suppressMessages(library(dplyr))

JAGS <- Sys.which("jags")
if (JAGS == "") JAGS <- "/opt/homebrew/bin/jags"
if (!file.exists(JAGS)) stop("No JAGS binary found. macOS: brew install jags")
if (!requireNamespace("coda", quietly = TRUE)) stop("Package 'coda' is required.")

set.seed(20151112)
N_CHAINS <- 5L; N_BURN <- 1000L; N_ITER <- 5000L   # as in 90_latent_states_jags.R

## Reuse the data preparation and model text from the rjags script, so the two
## cannot drift apart. Everything before the rjags call is pure data handling.
src <- readLines(file.path(PROJ, "R", "90_latent_states_jags.R"))
eval(parse(text = paste(
  src[grep("^counts <- read.table", src):grep('^\\}"$', src)[1]], collapse = "\n")),
  globalenv())

WORK <- file.path(OUT_DIR, "latent_states_cli"); dir.create(WORK, FALSE, TRUE)
writeLines(MODEL, file.path(WORK, "model.bug"))

split_rhat <- function(chains, v) {          # Gelman et al., split chains
  x <- sapply(chains, function(c) as.numeric(c[, v]))
  n <- nrow(x); h <- floor(n / 2)
  s <- cbind(x[1:h, ], x[(n - h + 1):n, ])
  W <- mean(apply(s, 2, var)); B <- h * var(colMeans(s))
  sqrt(((h - 1) / h * W + B / h) / W)
}

run_rep <- function(rep_id) {
  d <- counts[counts$rep == rep_id, ]; d <- d[order(d$time.h), ]
  Tobs <- d$total; Sobs <- d$spore
  transition <- d$time.h < PHASE_BREAK_H; equilibrium <- d$time.h > PHASE_BREAK_H
  v_naive <- Tobs - Sobs; usable <- v_naive > 0
  n <- nrow(d); n_transition <- sum(transition)
  s_alpha1 <- gamma_moments(Sobs[transition])[1];  s_beta1 <- gamma_moments(Sobs[transition])[2]
  s_alpha2 <- gamma_moments(Sobs[equilibrium])[1]; s_beta2 <- gamma_moments(Sobs[equilibrium])[2]
  v_alpha1 <- gamma_moments(v_naive[usable & transition])[1]
  v_beta1  <- gamma_moments(v_naive[usable & transition])[2]
  v_alpha2 <- gamma_moments(v_naive[usable & equilibrium])[1]
  v_beta2  <- gamma_moments(v_naive[usable & equilibrium])[2]

  old <- setwd(WORK); on.exit(setwd(old), add = TRUE)
  dump(c("n","n_transition","Sobs","Tobs","s_alpha1","s_beta1","s_alpha2","s_beta2",
         "v_alpha1","v_beta1","v_alpha2","v_beta2"), file = "data.R")
  writeLines(c('model in "model.bug"', 'data in "data.R"',
    sprintf("compile, nchains(%d)", N_CHAINS), "initialize",
    sprintf("update %d", N_BURN),
    "monitor sig_spore_obs", "monitor sig_total_obs",
    "monitor S, thin(10)", "monitor V, thin(10)",
    sprintf("update %d", N_ITER),
    'coda sig_spore_obs, stem("g1_")', 'coda sig_total_obs, stem("g2_")',
    'coda S, stem("S_")', 'coda V, stem("V_")', "exit"), "run.cmd")
  system2(JAGS, "run.cmd", stdout = FALSE, stderr = FALSE)

  rd <- function(stem) lapply(seq_len(N_CHAINS), function(i)
    coda::read.coda(sprintf("%schain%d.txt", stem, i),
                    sprintf("%sindex.txt", stem), quiet = TRUE))
  out <- do.call(rbind, lapply(c(g1_ = "g1_", g2_ = "g2_"), function(st) {
    ch <- rd(st); v <- coda::varnames(ch[[1]])[1]
    data.frame(replicate = rep_id, parameter = v,
               rhat = split_rhat(ch, v),
               ess  = as.numeric(coda::effectiveSize(coda::as.mcmc.list(ch))[1]))
  }))
  S <- colMeans(do.call(rbind, lapply(rd("S_"), as.matrix)))
  V <- colMeans(do.call(rbind, lapply(rd("V_"), as.matrix)))
  arc <- read.table(file.path(DATA_DIR, "latent_states",
                              sprintf("transition.Bayes.rep%s.txt", rep_id)), header = TRUE)
  out$cor_spore <- cor(arc$s.est, S); out$cor_veg <- cor(arc$v.est, V)
  rownames(out) <- NULL
  out
}

res <- do.call(rbind, lapply(c("A","B","C","D"), function(r) {
  message("  replicate ", r); run_rep(r) }))
write.csv(res, file.path(OUT_DIR, "latent_state_diagnostics.csv"), row.names = FALSE)

md <- c("## Table. Latent-state model: convergence and reproduction", "",
  sprintf(paste("%d chains, %s burn-in, %s iterations, the settings of the original fit.",
                "R-hat is the split potential scale reduction factor and ESS the effective",
                "sample size. The last two columns correlate the refit posterior means",
                "against the archived fits in `data/latent_states/`."),
          N_CHAINS, format(N_BURN, big.mark = ","), format(N_ITER, big.mark = ",")), "",
  "| Replicate | Parameter | R-hat | ESS | r (spore) | r (vegetative) |",
  "| :--- | :--- | ---: | ---: | ---: | ---: |",
  sprintf("| %s | %s | %.3f | %.0f | %.4f | %.4f |",
          res$replicate, gsub("_", "\\\\_", res$parameter), res$rhat, res$ess,
          res$cor_spore, res$cor_veg))
writeLines(md, file.path(OUT_DIR, "table_latent_diagnostics.md"))
print(res, row.names = FALSE, digits = 4)
message("\nWrote output/latent_state_diagnostics.csv and table_latent_diagnostics.md")
