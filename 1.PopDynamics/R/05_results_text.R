################################################################################
# 05_results_text.R
#
# Write the results paragraphs with every number pulled from the analysis, so
# the prose and the tables cannot drift apart.
#
# Input : output/population_descriptors.csv, output/sigmoidal_parms.csv
# Output: output/results_text.md
#
# The wording follows the manuscript draft. Numbers that differ from that draft
# are listed at the end of the file, with the reason.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the populationDynamics project root (or R/).")
source(.setup[1])

need <- file.path(OUT_DIR, c("population_descriptors.csv", "sigmoidal_parms.csv"))
if (!all(file.exists(need))) stop("Run R/02_sigmoidal_fits.R and R/04_report.R first.")

desc  <- read.csv(file.path(OUT_DIR, "population_descriptors.csv"))
parms <- read.csv(file.path(OUT_DIR, "sigmoidal_parms.csv"))

## ---- number formatting -----------------------------------------------------

# "1.33 x 10^6", the way it is written in the text.
sci <- function(x, digits = 2) {
  e <- floor(log10(abs(x)))
  sprintf("%.*f x 10^%d^", digits, x / 10^e, e)
}

pm     <- function(x, digits = 1) sprintf("%.*f ± %.*f", digits, mean(x), digits, sem(x))
pm_sci <- function(x, digits = 2) sprintf("%s ± %s", sci(mean(x), digits), sci(sem(x), 1))

pooled <- function(series, parameter) {
  parms$estimate[parms$series == series & parms$fit == "pooled" &
                   parms$parameter == parameter]
}

# How far a series travels between its two asymptotes, as a fold change.
fold_change <- function(series) {
  10^abs(pooled(series, "a") - pooled(series, "b"))
}

# Rounded to the nearest hundred, since the sentence says "nearly".
n_days <- signif(max(read.csv(file.path(DATA_DIR, "pop_dynamics.csv"))$time), 2)

## ---- the paragraphs --------------------------------------------------------

results <- c(
  "## Populations persist under resource limitation",
  "",
  "### 1. Energy limitation alters population dynamics",
  "",
  sprintf(paste(
    "Energy limitation led to a large drop in total population size (*T*).",
    "There was a %s%% reduction over the first 10 days, a %s-fold fall.",
    "Following this bottleneck, there was an extended period of equilibrium",
    "where total cell densities remained relatively unchanged for nearly %d",
    "days, converging on a population size of %s cell/mL (Fig. 1, Table 2)."),
    pm(desc$crash_percent, 2), pm(desc$crash_fold, 0), n_days,
    pm_sci(desc$equilibrium_total)),
  "",
  "### 2. Energy limitation creates population heterogeneity",
  "",
  sprintf(paste(
    "The population dynamics involved a large transition from mostly",
    "vegetative cell types to spores. In other words, energy limitation",
    "generated a shift in the equilibrium abundances of the two",
    "sub-populations. This is not entirely surprising, since the *Bacillus*",
    "stress response is known to let individuals form endospores under extreme",
    "energy limitation, and it is also known that not all individuals commit to",
    "sporulation synchronously. Despite endosporulation being understood as a",
    "responsive process, there appears to be some degree of bet-hedging,",
    "albeit consistent across populations. During the first 24 h, sporulation",
    "efficiency was low (%s%%), but %s%% of remaining individuals were spores",
    "by day 10 (Fig. 1)."),
    pm(desc$spore_eff_day1, 3), pm(desc$spore_eff_day10, 1)),
  "",
  sprintf(paste(
    "This was not due entirely to the death of vegetative cells: the absolute",
    "abundance of spores rose %.0f-fold over the experiment, from %s to %s",
    "CFU/mL, while non-spores fell %.0f-fold. Maximum-likelihood fitting of a",
    "sigmoidal function put the midpoint of this metabolic transition at %s d",
    "following the onset of energy limitation for the spore sub-population and",
    "%s d for the non-spore sub-population. By the end of the experiment,",
    "%s%% of cells were spores and %s%% were vegetative. A non-trivial",
    "fraction of a *B. subtilis* population can therefore withstand extended",
    "resource limitation as vegetative cells, whether actively dividing on",
    "resources turned over by dead cells and cannibalism, or in some non-spore",
    "form of dormancy."),
    fold_change("spore"), sci(10^pooled("spore", "b")), sci(10^pooled("spore", "a")),
    fold_change("veg"),
    pm(desc$midpoint_spore_d, 2), pm(desc$midpoint_veg_d, 2),
    pm(desc$final_spore_pct, 0), pm(desc$final_veg_pct, 0)),
  "",
  "---",
  "",
  "### Where these numbers differ from the earlier draft",
  "",
  paste(
    "Every value above is generated from the analysis rather than typed in, so",
    "it moves if the data or the fits do. Four differences from the draft are",
    "worth knowing about."),
  "",
  sprintf(paste(
    "- **Averages across replicates.** The original script used `mean(a, b, c,",
    "d)` to average the four replicate populations. R takes only the first",
    "argument as data and passes the rest to `trim` and `na.rm`, so those",
    'means were replicate A alone. The bottleneck is %s%% rather than 99.6%%,',
    "and day-10 sporulation efficiency %s%% rather than 79.3%%."),
    pm(desc$crash_percent, 2), pm(desc$spore_eff_day10, 1)),
  paste(
    "- **Replicate D.** The original read replicate C's state estimates in",
    "place of D's, so D never contributed its own spore and non-spore numbers.",
    "Fixed here; D's own estimates are used."),
  paste(
    "- **Uncertainty on the percentages.** The draft reports 99.6 +/- 0.0005%",
    "and 0.13 +/- 0.021%, mixing a percentage with the standard error of the",
    "underlying proportion. Both are on the same scale here."),
  sprintf(paste(
    "- **Transition midpoint.** The draft gives 4.4 ± 2.27 d. The midpoint",
    "itself reproduces (%s d for the spore sub-population), but the ± 2.27",
    "could not be reproduced from these fits by any route tried -- the",
    "standard error across replicates is %.2f d. Worth checking against",
    "whatever produced the draft before reusing it."),
    pm(desc$midpoint_spore_d, 2), sem(desc$midpoint_spore_d)),
  ""
)

writeLines(results, file.path(OUT_DIR, "results_text.md"))
message("Wrote output/results_text.md")

cat("\n", paste(results, collapse = "\n"), "\n", sep = "")
