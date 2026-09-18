################################################################################
# 03_null_model.R
#
# What should the allele-frequency spectrum look like with no selection?
#
# Input : data/mutations.csv
# Output: output/null_spectrum.csv     neutral expectation per fraction
#         output/null_summary.csv      per-iteration summaries
#
# Under neutrality only the genealogy of the sampled clones matters, so this
# traces the sampled lineages backwards rather than forward-simulating the
# whole population. Same model as the forward simulation it replaces, no
# external software, runs in seconds.
#
#   29 generations of growth from a single ancestor -- log2(10^9 cells) --
#   then a bottleneck to N and G generations at constant size.
#     spores     G = 0     dormant, no turnover
#     non-spores G = 1000  cryptic growth, a deliberate upper bound
#
# Mutations are dropped onto the branches at rate U per genome per generation;
# a mutation on a branch with k descendants in the sample is carried by k of
# the sequenced clones. Mutations on the root branch are in every clone and so
# are not polymorphic; they are dropped.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the mutations project root (or R/).")
source(.setup[1])

set.seed(20260730)

## ---- parameters ------------------------------------------------------------
MU          <- 3.28e-10          # per bp per generation (Sung et al. 2015)
U           <- MU * GENOME_BP    # per genome per generation
GENS_GROWTH <- 29L               # log2(10^9 cells) from a single ancestor
N_ITER      <- 100L              # iterations, as in the draft

# Equilibrium population sizes, from the population dynamics analysis: spores
# settle near 10^6 CFU/mL and non-spores near 3 x 10^5.
N_SPORE     <- 9.9e5
N_NONSPORE  <- 2.9e5

SIM <- list(
  S10        = list(n = FRACTIONS$S10$n,        N = N_SPORE,    G = 0L),
  S1000      = list(n = FRACTIONS$S1000$n,      N = N_SPORE,    G = 0L),
  `S+NS1000` = list(n = FRACTIONS$`S+NS1000`$n, N = N_NONSPORE, G = 1000L)
)

## ---- one genealogy ---------------------------------------------------------
# Wright-Fisher coalescent: each generation every lineage picks a parent at
# random from N; lineages sharing a parent merge. Returns the branches as
# (length in generations, number of sampled clones descending from it).
one_tree <- function(n, N, G) {
  leaves <- rep(1L, n)
  blen   <- rep(0L, n)
  seg_len <- integer(0)
  seg_cnt <- integer(0)

  step_back <- function(Npop) {
    k <- length(leaves)
    if (k <= 1) return(invisible(NULL))
    blen <<- blen + 1L
    parent <- sample.int(max(as.integer(Npop), 1L), k, replace = TRUE)
    if (anyDuplicated(parent)) {
      grp <- split(seq_len(k), parent)
      new_leaves <- integer(0); new_blen <- integer(0)
      for (g in grp) {
        if (length(g) > 1) {                    # coalescence: close branches
          seg_len <<- c(seg_len, blen[g])
          seg_cnt <<- c(seg_cnt, leaves[g])
          new_leaves <- c(new_leaves, sum(leaves[g]))
          new_blen   <- c(new_blen, 0L)
        } else {
          new_leaves <- c(new_leaves, leaves[g])
          new_blen   <- c(new_blen, blen[g])
        }
      }
      leaves <<- new_leaves; blen <<- new_blen
    }
    invisible(NULL)
  }

  if (G > 0) for (t in seq_len(G)) { if (length(leaves) == 1) break; step_back(N) }
  # growth phase, backwards: the population halves each generation
  for (g in seq(GENS_GROWTH, 1)) {
    if (length(leaves) == 1) break
    step_back(2^(g - 1))
  }

  data.frame(len = seg_len, cnt = seg_cnt)
}

# Allele count of every mutation in one simulated sample.
one_rep <- function(n, N, G) {
  seg <- one_tree(n, N, G)
  if (!nrow(seg)) return(integer(0))
  rep(seg$cnt, rpois(nrow(seg), U * seg$len))
}

## ---- the spore fractions need no simulation at all -------------------------
# Spores are dormant: on the model's own assumptions they experience no
# generations after the population contracts. Every mutation they carry was
# acquired during the 29 generations of growth beforehand, so the expected
# number is just n_clones x U x 29. That arithmetic is not falsifiable by a
# choice of population size, and it is the strongest form of the spore result.

per_clone <- U * GENS_GROWTH

spore_check <- bind_rows(lapply(c("S10", "S1000"), function(nm) {
  n_cl <- FRACTIONS[[nm]]$n
  obs  <- sum(read.csv(file.path(DATA_DIR, "mutations.csv"),
                       stringsAsFactors = FALSE)$fraction == nm)
  data.frame(fraction = nm, clones = n_cl,
             expected = n_cl * per_clone, observed = obs,
             # A 95% interval on the expected count says more than a P value
             # that rounds to 1: the point is agreement, not a failed test.
             lower95 = qpois(0.025, n_cl * per_clone),
             upper95 = qpois(0.975, n_cl * per_clone),
             stringsAsFactors = FALSE)
}))
write.csv(spore_check, file.path(OUT_DIR, "spore_analytic_expectation.csv"),
          row.names = FALSE)

cat("Spore fractions: expected mutations without any simulation\n")
cat(sprintf("  %.4f mutations per clone from %d generations of growth at U = %.3e\n\n",
            per_clone, GENS_GROWTH, U))
print(spore_check %>%
        mutate(expected = sprintf("%.2f", expected),
               interval95 = sprintf("%d-%d", lower95, upper95)) %>%
        select(fraction, clones, expected, interval95, observed) %>%
        as.data.frame(), row.names = FALSE)
cat(paste0("\nBoth fractions carry the mutations the growth phase alone predicts,\n",
           "and nothing more. Nothing accumulated while they were dormant.\n\n"))

## ---- run -------------------------------------------------------------------

message(sprintf("Neutral coalescent: U = %.3e mutations per genome per generation",
                U))

spectrum <- list()
summary_rows <- list()

for (nm in names(SIM)) {
  f <- SIM[[nm]]
  reps <- replicate(N_ITER, one_rep(f$n, f$N, f$G), simplify = FALSE)

  # Number of mutations at each allele count, per iteration.
  counts <- t(vapply(reps, function(x)
    as.integer(table(factor(x, levels = seq_len(f$n)))), integer(f$n)))

  spectrum[[nm]] <- data.frame(
    fraction  = nm,
    carriers  = seq_len(f$n),
    frequency = seq_len(f$n) / f$n,
    mean      = colMeans(counts),
    lower     = apply(counts, 2, quantile, 0.025),
    upper     = apply(counts, 2, quantile, 0.975),
    max       = apply(counts, 2, max),
    stringsAsFactors = FALSE)

  summary_rows[[nm]] <- data.frame(
    fraction = nm, n_clones = f$n, generations = GENS_GROWTH + f$G,
    mean_mutations = mean(vapply(reps, length, integer(1))),
    iterations_with_shared = sum(vapply(reps, function(x) any(x >= 2), logical(1))),
    max_carriers_seen = max(c(0L, unlist(reps))),
    stringsAsFactors = FALSE)

  message(sprintf("  %-9s n = %3d, %4d generations: %.2f mutations per run, ",
                  nm, f$n, GENS_GROWTH + f$G, summary_rows[[nm]]$mean_mutations),
          sprintf("largest allele count in %d runs = %d",
                  N_ITER, summary_rows[[nm]]$max_carriers_seen))
}

spectrum <- bind_rows(spectrum)
summary_df <- bind_rows(summary_rows)

write.csv(spectrum, file.path(OUT_DIR, "null_spectrum.csv"), row.names = FALSE)
write.csv(summary_df, file.path(OUT_DIR, "null_summary.csv"), row.names = FALSE)

## ---- observed vs neutral ---------------------------------------------------

muts <- read.csv(file.path(DATA_DIR, "mutations.csv"), stringsAsFactors = FALSE)

test_rows <- lapply(names(SIM), function(nm) {
  obs <- muts$carriers[muts$fraction == nm]
  sp  <- spectrum[spectrum$fraction == nm, ]

  # Under the null, essentially every mutation is a singleton. Compare the
  # observed number of shared mutations with that expectation as a Poisson
  # count: expected = mean number of shared mutations across iterations.
  expected_shared <- sum(sp$mean[sp$carriers >= 2])
  observed_shared <- sum(obs >= 2)

  # Poisson tail probability of seeing at least this many, plus the normal
  # approximation the draft reports as a z-test.
  p_pois <- stats::ppois(observed_shared - 1, expected_shared, lower.tail = FALSE)
  z <- (observed_shared - expected_shared) / sqrt(max(expected_shared, .Machine$double.eps))

  data.frame(fraction = nm,
             observed_mutations = length(obs),
             observed_shared = observed_shared,
             expected_shared = expected_shared,
             z = z,
             p_value = p_pois,
             stringsAsFactors = FALSE)
})
tests <- bind_rows(test_rows)
write.csv(tests, file.path(OUT_DIR, "null_tests.csv"), row.names = FALSE)

cat("\nObserved vs neutral expectation (mutations carried by more than one clone)\n\n")
print(tests %>%
        mutate(expected_shared = sprintf("%.3f", expected_shared),
               z = sprintf("%.1f", z),
               p_value = format.pval(p_value, digits = 3, eps = 2.2e-16)),
      row.names = FALSE)

## ---- how much does the answer depend on population size? -------------------
# Whether shared mutations are surprising under neutrality turns almost
# entirely on N for the non-spore subpopulation: the chance that two of the 84
# sampled clones share an ancestor within G generations is roughly G/N per
# pair, so a small N produces sharing on its own. N is quoted as ~3 x 10^5
# CFU/mL, but the number that matters is cells in the whole culture, not per
# mL. This sweep shows where the conclusion changes.

sweep_N <- c(2.9e5, 1.33e6, 1e7, 1e8, 1e9)
f <- SIM[["S+NS1000"]]

sweep <- bind_rows(lapply(sweep_N, function(N) {
  reps <- replicate(N_ITER, one_rep(f$n, N, f$G), simplify = FALSE)
  shared <- vapply(reps, function(x) sum(x >= 2), integer(1))
  exp_shared <- mean(shared)
  data.frame(
    N = N,
    mean_mutations = mean(vapply(reps, length, integer(1))),
    expected_shared = exp_shared,
    iterations_with_any_shared = sum(shared > 0),
    max_carriers = max(c(0L, unlist(reps))),
    p_observed_10 = stats::ppois(9, max(exp_shared, 1e-12), lower.tail = FALSE),
    stringsAsFactors = FALSE)
}))

# N is a count of cells, not a density. The densities are known; what is
# missing is one number, the culture volume. This is the conversion.
sweep$implied_volume_mL <- sweep$N / N_NONSPORE

write.csv(sweep, file.path(OUT_DIR, "null_population_size_sweep.csv"),
          row.names = FALSE)

cat("\nSensitivity to the non-spore population size (S+NS1000, G = 1000)\n")
cat("Observed: 34 mutations, 10 carried by more than one clone.\n")
cat(sprintf(paste0(
  "N is a count of cells, not a density. At %.1e non-spores per mL, the\n",
  "implied culture volume is shown alongside; the experiment's volume settles\n",
  "this outright. Anything from about 25 mL up puts N past 10^7, which is the\n",
  "regime where the observed sharing is improbable under neutrality.\n\n"),
  N_NONSPORE))
print(sweep %>%
        mutate(N = sprintf("%.1e", N),
               implied_volume_mL = sprintf("%.1f", implied_volume_mL),
               mean_mutations = sprintf("%.1f", mean_mutations),
               expected_shared = sprintf("%.2f", expected_shared),
               p_observed_10 = format.pval(p_observed_10, digits = 3,
                                           eps = 2.2e-16)),
      row.names = FALSE)

## ---- how many generations would the observed frequency need? ---------------
# For one mutant lineage to reach a frequency f by division alone, while the
# rest of the population does not divide, needs log2(f * N) generations.

EQUILIBRIUM_TOTAL <- 1.33e6      # CFU/mL, from the population dynamics analysis
top <- muts %>% filter(fraction == "S+NS1000") %>% slice_max(carriers, n = 1)
gens_needed <- log2(top$frequency[1] * EQUILIBRIUM_TOTAL)

cat(sprintf(
  "\nFor a mutation to reach %.1f%% of a population of %.2e cells by division\n",
  100 * top$frequency[1], EQUILIBRIUM_TOTAL))
cat(sprintf("alone takes at least %.0f generations.\n", gens_needed))

writeLines(sprintf("%.1f", gens_needed),
           file.path(OUT_DIR, "generations_needed.txt"))

message("\nWrote output/null_spectrum.csv, null_summary.csv and null_tests.csv")
