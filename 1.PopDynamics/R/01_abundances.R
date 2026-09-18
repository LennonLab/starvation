################################################################################
# 01_abundances.R
#
# Turn colony counts into abundances, and attach the modelled spore and
# vegetative states.
#
# Input : data/spore.transition.txt          colony counts over ~900 days
#         data/spore.transition.ltde.txt     counts at the end of the experiment
#         data/latent_states/transition.Bayes.rep{A,B,C,D}.txt
# Output: data/pop_dynamics.csv              the analysis table
#         output/negative_vegetative.csv     where T - S goes negative
#
# Total and spore abundances are counted directly. Vegetative cells are not:
# they are what is left after subtracting spores from the total, and because
# both counts carry plating error that difference goes negative at a fair
# number of timepoints. That is why the vegetative series is modelled rather
# than subtracted -- see R/90_latent_states_jags.R for the state model, whose
# archived output is read here.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the populationDynamics project root (or R/).")
source(.setup[1])

## ---- 1. counts -> CFU/mL ---------------------------------------------------

counts <- read.table(file.path(DATA_DIR, "spore.transition.txt"), header = TRUE)
names(counts)[1] <- "time.h"

counts <- counts %>%
  mutate(
    time   = time.h / 24,                                        # days
    total  = total * 10^total_dil * PLATE_VOLUME_FACTOR,         # CFU/mL
    spore  = spore * 10^spor_dil  * PLATE_VOLUME_FACTOR,
    veg_naive = total - spore,                                   # before modelling
    phase  = ifelse(time.h < PHASE_BREAK_H, "transition", "equilibrium")
  )

# Sampling stopped between 13292 h and 14108 h (554-588 d) while the technician
# was away; the gap is real, not missing data.

neg <- counts %>%
  filter(veg_naive <= 0) %>%
  select(time.h, time, rep, total, spore, veg_naive)

write.csv(neg, file.path(OUT_DIR, "negative_vegetative.csv"), row.names = FALSE)

# Strictly negative counts reproduce the ones noted in the original script
# (A=23, B=31, C=49, D=48). Exact zeros are just as unusable on a log axis.
n_neg  <- tabulate(factor(counts$rep[counts$veg_naive <  0], levels = REPS))
n_zero <- tabulate(factor(counts$rep[counts$veg_naive == 0], levels = REPS))

message(sprintf(
  "Naive T - S is negative at %d and zero at %d of %d timepoints (negatives %s);\n  that is why V is modelled rather than subtracted.",
  sum(n_neg), sum(n_zero), nrow(counts),
  paste(sprintf("%s=%d", REPS, n_neg), collapse = ", ")))

## ---- 2. modelled spore and vegetative states -------------------------------

# One file per replicate population, each holding the posterior mean and 95%
# interval for the spore (s) and vegetative (v) state at every timepoint.
#
# NOTE: the original script read transition.Bayes.repC.txt twice -- once as C
# and once as D -- so replicate D's own estimates never reached the figures.
# The archived repD.txt is fine; only the reading was wrong. Each replicate is
# read from its own file here.
# The state files record time in days written to 15 significant digits, so
# 186.833333333333 there is not bit-identical to 4484/24 here. Join on whole
# hours, which both sides represent exactly.
read_states <- function(rep) {
  f <- file.path(DATA_DIR, "latent_states",
                 sprintf("transition.Bayes.rep%s.txt", rep))
  read.table(f, header = TRUE) %>%
    mutate(rep = rep, time.h = round(time * 24)) %>%
    select(-time) %>%
    rename(spore_est = s.est, spore_lo = s.lower, spore_hi = s.upper,
           veg_est   = v.est, veg_lo   = v.lower, veg_hi   = v.upper)
}

states <- bind_rows(lapply(REPS, read_states)) %>%
  # back onto the CFU/mL scale the model was taken off
  mutate(across(c(spore_est, spore_lo, spore_hi,
                  veg_est, veg_lo, veg_hi), ~ .x * MODEL_SCALE))

## ---- 3. join ---------------------------------------------------------------

pop <- counts %>%
  select(time.h, time, rep, phase, total, spore_obs = spore, veg_naive) %>%
  left_join(states, by = c("time.h", "rep"))

if (anyNA(pop$spore_est)) {
  stop("Some timepoints have no modelled state -- check that the count and ",
       "state files cover the same times.")
}

write.csv(pop, file.path(DATA_DIR, "pop_dynamics.csv"), row.names = FALSE)

## ---- 4. end-of-experiment counts -------------------------------------------
# Read at ~1580 and ~1910 days, well past the end of the main series. Not part
# of the figure, but kept alongside it.

end <- read.table(file.path(DATA_DIR, "spore.transition.ltde.txt"), header = TRUE) %>%
  mutate(time  = time.h / 24,
         total = total * 10^total_dil * PLATE_VOLUME_FACTOR,
         spore = spore * 10^spor_dil  * PLATE_VOLUME_FACTOR,
         veg_naive = total - spore)

write.csv(end, file.path(DATA_DIR, "pop_dynamics_end.csv"), row.names = FALSE)

message(sprintf("Wrote data/pop_dynamics.csv (%d rows, %d timepoints x %d replicates)",
                nrow(pop), length(unique(pop$time)), length(REPS)))
