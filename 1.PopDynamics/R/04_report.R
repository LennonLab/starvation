################################################################################
# 04_report.R
#
# Population descriptors and the manuscript tables.
#
# Input : data/pop_dynamics.csv, output/sigmoidal_parms.csv
# Output: output/population_descriptors.csv
#         output/table1_sigmoidal_parms.md
#         output/table2_descriptors.md
#         output/methods_stats.md
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the populationDynamics project root (or R/).")
source(.setup[1])

pop   <- read.csv(file.path(DATA_DIR, "pop_dynamics.csv"))
parms <- read.csv(file.path(OUT_DIR, "sigmoidal_parms.csv"))

# The bottleneck is over by day 10; everything after it is the equilibrium.
EQUILIBRIUM_FROM <- 10
# "End of the experiment" for the final composition: the last stretch of it.
FINAL_WINDOW_D   <- 60

md_table <- function(df) {
  numeric_only <- vapply(df, function(col)
    all(grepl("^[-–<>0-9.%]+$", as.character(col))), logical(1))
  align <- ifelse(numeric_only, "r", "l")
  sep <- ifelse(align == "l", ":---", "---:")
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("| ", paste(sep, collapse = " | "), " |"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

write_md <- function(lines, file, title, caption) {
  writeLines(c(paste("##", title), "", caption, "", lines, ""),
             file.path(OUT_DIR, file))
  message("Wrote output/", file)
}

## ---- population descriptors ------------------------------------------------

# Value of one column for one replicate at one time, to the nearest reading.
at_time <- function(rep, day, column) {
  d <- pop[pop$rep == rep, ]
  d[[column]][which.min(abs(d$time - day))]
}

# The midpoint of a fitted transition, in days: the curve is fit against
# log10(time), so the parameter m comes back as 10^m.
midpoint_days <- function(series, fit) {
  10^parms$estimate[parms$series == series & parms$fit == fit &
                      parms$parameter == "m"]
}

desc <- bind_rows(lapply(REPS, function(r) {
  d   <- pop[pop$rep == r, ]
  t0  <- at_time(r, 0,  "total")
  t10 <- at_time(r, 10, "total")

  eq    <- d[d$time >= EQUILIBRIUM_FROM, ]
  final <- d[d$time > max(d$time) - FINAL_WINDOW_D, ]

  data.frame(
    rep = r,
    total_day0  = t0,
    total_day10 = t10,
    # The bottleneck: how far the population falls in the first ten days.
    crash_fold    = t0 / t10,
    crash_percent = 100 * (t0 - t10) / t0,
    # Where total abundance settles once the bottleneck is over.
    equilibrium_total = mean(eq$total),
    # Sporulation efficiency: the fraction of the population that is spores.
    spore_eff_day1  = 100 * at_time(r, 1,  "spore_est") / at_time(r, 1,  "total"),
    spore_eff_day10 = 100 * at_time(r, 10, "spore_est") / at_time(r, 10, "total"),
    # When the vegetative-to-spore transition happens, from the fitted curves.
    midpoint_spore_d = midpoint_days("spore", r),
    midpoint_veg_d   = midpoint_days("veg", r),
    # Composition over the last stretch of the experiment.
    final_spore_pct = 100 * mean(final$spore_est / (final$spore_est + final$veg_est)),
    final_veg_pct   = 100 * mean(final$veg_est   / (final$spore_est + final$veg_est))
  )
}))

write.csv(desc, file.path(OUT_DIR, "population_descriptors.csv"), row.names = FALSE)

summarise_col <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean(x), sem(x))
}

ROWS <- list(
  list(name = "Total at day 0 (CFU/mL)",                 col = "total_day0",        fmt = "sci"),
  list(name = "Total at day 10 (CFU/mL)",                col = "total_day10",       fmt = "sci"),
  list(name = "Bottleneck (fold decline, days 0–10)",    col = "crash_fold",        digits = 0),
  list(name = "Bottleneck (% decline)",                  col = "crash_percent",     digits = 2),
  list(name = "Equilibrium total (CFU/mL, after day 10)", col = "equilibrium_total", fmt = "sci"),
  list(name = "Sporulation efficiency at day 1 (%)",     col = "spore_eff_day1",    digits = 3),
  list(name = "Sporulation efficiency at day 10 (%)",    col = "spore_eff_day10",   digits = 1),
  list(name = "Spore transition midpoint (d)",           col = "midpoint_spore_d",  digits = 2),
  list(name = "Non-spore transition midpoint (d)",       col = "midpoint_veg_d",    digits = 2),
  list(name = "Spores at end of experiment (%)",         col = "final_spore_pct",   digits = 1),
  list(name = "Non-spores at end of experiment (%)",     col = "final_veg_pct",     digits = 1)
)

cell <- function(x, row) {
  if (identical(row$fmt, "sci")) sprintf("%.2e", x)
  else sprintf(paste0("%.", row$digits, "f"), x)
}

summary_cell <- function(row) {
  x <- desc[[row$col]]
  if (identical(row$fmt, "sci")) sprintf("%.2e ± %.1e", mean(x), sem(x))
  else summarise_col(x, row$digits)
}

t2 <- data.frame(
  Quantity = vapply(ROWS, `[[`, character(1), "name"),
  `Mean ± SEM` = vapply(ROWS, summary_cell, character(1)),
  A = vapply(ROWS, function(r) cell(desc[[r$col]][1], r), character(1)),
  B = vapply(ROWS, function(r) cell(desc[[r$col]][2], r), character(1)),
  C = vapply(ROWS, function(r) cell(desc[[r$col]][3], r), character(1)),
  D = vapply(ROWS, function(r) cell(desc[[r$col]][4], r), character(1)),
  check.names = FALSE
)

write_md(
  md_table(t2), "table2_descriptors.md",
  "Table 2. Population descriptors",
  paste0(
    "Per replicate population and averaged across the four (mean ± SEM). ",
    "The bottleneck is the decline in total abundance over the first ten ",
    "days, after which the population settles. Sporulation efficiency is the ",
    "modelled spore abundance as a percentage of the counted total. Transition ",
    "midpoints are 10^m from the fitted sigmoidal curves. End-of-experiment ",
    "composition is averaged over the last 60 days."))

## ---- sigmoidal parameter table ---------------------------------------------

wide <- parms %>%
  mutate(value = sprintf("%.3f (%.3f)", estimate, se)) %>%
  select(label, fit, parameter, value) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  arrange(match(label, c("Total (S + V)", "Spore (S)", "Vegetative (V)")),
          match(fit, c("pooled", REPS)))

names(wide)[1:2] <- c("Series", "Fit")

write_md(
  md_table(as.data.frame(wide)), "table1_sigmoidal_parms.md",
  "Table 1. Sigmoidal fits to the abundance series",
  paste0(
    "Maximum-likelihood estimates (standard errors) for ",
    "log10(N) = b + (a − b) / (1 + exp((m − log10(t)) / w)), fit to all four ",
    "replicate populations together (`pooled`) and to each on its own. ",
    "*a* is the level the population settles at and *b* the level it starts ",
    "from, both log10 CFU/mL; *m* is the midpoint in log10 days, *w* the ",
    "width of the transition, and *z* the residual standard deviation."))

## ---- methods ---------------------------------------------------------------

n_times <- length(unique(pop$time))
n_days  <- round(max(pop$time))
n_neg   <- sum(pop$veg_naive <= 0)
pkg <- function(p) as.character(utils::packageVersion(p))

methods <- c(
  "### Population sampling",
  "",
  sprintf(paste(
    "Four replicate populations (A–D) were sampled at %d timepoints over %d",
    "days. At each timepoint, total and spore abundances were measured as",
    "colony-forming units: total counts from a plated dilution series, spore",
    "counts from the same sample after treatment that leaves only spores",
    "viable. Counts were converted to CFU per mL by multiplying by the",
    "dilution factor and by ten, the plated volume being 100 µL. Sampling",
    "paused between days 554 and 588."), n_times, n_days),
  "",
  "### Vegetative abundance",
  "",
  sprintf(paste(
    "Vegetative cells were not counted directly: they are the difference",
    "between the total and the spore count, and because both counts carry",
    "plating error that difference is non-positive at %d of %d",
    "replicate-timepoints. Vegetative and spore abundances were therefore",
    "treated as latent states inferred jointly. For timepoint *i*, the spore",
    "state S[i] and vegetative state V[i] were given gamma priors, and the",
    "observations modelled as Sobs[i] ~ Normal(S[i], tau_1) and Tobs[i] ~",
    "Normal(V[i] + S[i], tau_2), so that the counted total constrains the sum",
    "and the counted spores constrain one component. Gamma shape and rate",
    "were set by moment-matching the observed abundances, separately for the",
    "transition phase (before day 10) and the equilibrium phase after it,",
    "since the two differ by orders of magnitude. Observation standard",
    "deviations were given uniform(0, 1000) priors. The model was fit per",
    "replicate in JAGS with five chains, 1,000 burn-in iterations and 5,000",
    "samples; posterior means are used as the state estimates."),
    n_neg, nrow(pop)),
  "",
  "### Abundance trajectories",
  "",
  paste(
    "Each series was described with a sigmoidal curve on log-log axes,",
    "log10(N) = b + (a − b) / (1 + exp((m − log10(t)) / w)), where *b* is the",
    "level the population starts from, *a* the level it settles at, *m* the",
    "midpoint in log10 days and *w* the width of the transition. Curves were",
    "fit by maximum likelihood assuming normal error on log10(N), to the four",
    "replicates pooled and to each replicate separately. Day-zero readings",
    "were retained: log10(0) sends the curve to exactly *b*, so they inform",
    "the starting asymptote and nothing else."),
  "",
  "### Software",
  "",
  sprintf(paste(
    "Analyses were run in %s with bbmle %s, ggplot2 %s and patchwork %s.",
    "The latent-state model was fit in JAGS via rjags. Code and data are in",
    "the populationDynamics project."),
    R.version.string, pkg("bbmle"), pkg("ggplot2"), pkg("patchwork")),
  ""
)

writeLines(c("## Methods: statistics", "", methods),
           file.path(OUT_DIR, "methods_stats.md"))
message("Wrote output/methods_stats.md")

## ---- report ----------------------------------------------------------------

cat("\nPopulation descriptors (mean ± SEM across replicates A-D)\n\n")
cat(sprintf("  bottleneck, days 0-10      %.0f-fold (%.2f%% ± %.2f%%)\n",
            mean(desc$crash_fold), mean(desc$crash_percent), sem(desc$crash_percent)))
cat(sprintf("  sporulation eff., day 1    %.3f%% ± %.3f%%\n",
            mean(desc$spore_eff_day1), sem(desc$spore_eff_day1)))
cat(sprintf("  sporulation eff., day 10   %.1f%% ± %.1f%%\n",
            mean(desc$spore_eff_day10), sem(desc$spore_eff_day10)))
