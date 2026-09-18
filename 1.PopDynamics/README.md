# SporeMut population dynamics

Recreation of the population dynamics analysis from the SporeMut project: how
four replicate *Bacillus subtilis* populations collapse and then persist over
~900 days, and how they split between spores and vegetative cells.

## Run it

```
Rscript run_all.R
```

A few seconds. Needs `dplyr`, `tidyr`, `ggplot2`, `bbmle` and `patchwork`.
Everything lands in `output/`.

## What the experiment measures, and what it doesn't

Two things are counted at each timepoint: **total** abundance and **spore**
abundance, both as colonies from a plated dilution series.

**Vegetative cells are never counted.** They are what is left after
subtracting spores from the total — and because both counts carry plating
error, that difference is negative at 151 of 1844 replicate-timepoints and
exactly zero at another 23. You cannot plot that on a log axis, let alone fit
it. So spores and vegetative cells are treated as *latent states* inferred
jointly from the two counts: the total constrains their sum, the spore count
constrains one part, and gamma priors keep both positive. That model is
`R/90_latent_states_jags.R`; its archived output is what the rest reads.

`output/negative_vegetative.csv` lists every timepoint where the naive
subtraction fails, if you want to see the problem directly.

## Layout

```
data/
  spore.transition.txt        colony counts, 4 replicates x 463 timepoints
  spore.transition.ltde.txt   counts at ~1580 and ~1910 days, after the series
  latent_states/              archived per-replicate state estimates (JAGS)
  pop_dynamics.csv            built by R/01: the analysis table
R/
  00_setup.R              paths, packages, theme, constants
  01_abundances.R         counts -> CFU/mL, joined to the latent states
  02_sigmoidal_fits.R     abundance series -> sigmoidal curve fits
  03_figures.R            the two-panel figure
  04_report.R             descriptors, tables, methods text
  05_results_text.R       the results paragraphs
  90_latent_states_jags.R re-run the state model  (needs JAGS, not in run_all)
run_all.R
```

## Output

| file | what it is |
| --- | --- |
| `fig_spore_vegetative.pdf` / `.png` | spore and non-spore abundance — the manuscript panel |
| `fig_total.pdf` / `.png` | total abundance (S + V), counted directly |
| `fig_pop_dynamics.pdf` / `.png` | both stacked, the archived two-panel layout |
| `results_text.md` | the results paragraphs, numbers filled in from the analysis |
| `table1_sigmoidal_parms.md` | curve parameters, pooled and per replicate |
| `table2_descriptors.md` | bottleneck, equilibrium, sporulation efficiency, transition midpoints |
| `methods_stats.md` | the statistics paragraph for the methods section |
| `sigmoidal_parms.csv`, `sigmoidal_curves.csv`, `population_descriptors.csv` | the same numbers unrounded |
| `negative_vegetative.csv` | timepoints where total − spore is non-positive |

`results_text.md` writes the prose with every number pulled from the analysis
rather than typed in, so the text and the tables cannot drift apart. It ends
with a short list of where those numbers differ from the earlier draft, and
why.

## The analysis

**Counts to abundances.** CFU/mL = colonies × 10^dilution × 10, the ten
converting a 100 µL plated volume to 1 mL. The state model was fit on counts
divided by 10^4 to keep the numbers off the shoulders of the gamma priors;
`R/01_abundances.R` multiplies back.

**Latent states.** Per replicate, for timepoint *i*:

```
S[i] ~ Gamma(s_alpha, s_beta)          spore state
V[i] ~ Gamma(v_alpha, v_beta)          vegetative state
Sobs[i] ~ Normal(S[i],        tau_spore_obs)
Tobs[i] ~ Normal(V[i] + S[i], tau_total_obs)
```

Gamma shape and rate come from moment-matching the observed abundances,
separately for the transition phase (before 240 h) and the equilibrium phase
after it, since the two differ by orders of magnitude. Five chains, 1,000
burn-in, 5,000 samples; posterior means are the state estimates.

**Trajectories.** Each series is fit with a sigmoidal curve on log-log axes,

```
log10(N) = b + (a - b) / (1 + exp((m - log10(t)) / w))
```

by maximum likelihood, pooled across replicates and per replicate. Pooled:
total falls 116-fold from 1.5×10⁸ to 1.3×10⁶ CFU/mL, vegetative cells fall
201-fold, and spores rise 12-fold to about 10⁶.

## Things worth knowing

**Replicate D was reading replicate C's file.** The original script had

```r
data.B.D <- read.table(".../transition.Bayes.repC.txt")   # then labelled "D"
```

so replicate D's own state estimates never reached the figures or the fitted
parameters. The archived `transition.Bayes.repD.txt` is fine — only the
reading was wrong. Fixed here, and the footprint is visible: in the archived
`sigmoidal.parms.txt`, every spore and vegetative parameter for D is identical
to C. Refitting with D's own data reproduces all 60 archived parameters to
within 0.06 **except** D's spore and vegetative rows, which is exactly the set
the bug touched.

**Two more slips in the descriptor block, both fixed.** `crash.D.fold` divided
replicate C's day-zero total by D's day-ten total (333-fold; the correct
value is 489-fold). And `mean(a, b, c, d)` was used to average four
replicates, which in R returns `a` and passes the rest to `trim` and `na.rm` —
so every "mean across replicates" in the original was just replicate A. The
corrected averages are close but not identical: 99.70% vs 99.6% for the
bottleneck, 81.5% vs 79.3% for day-10 sporulation efficiency. Replicates A, B
and C reproduce the original per-replicate numbers exactly.

**The optimiser needed bounds.** Unconstrained, `mle2` wandered into negative
standard deviations (hundreds of `NaNs produced` warnings) and, for the total
series of replicate B, off to a settle level of 10³⁵ CFU/mL. Wide bounds on
all five parameters fix both without moving any of the good fits.

**`longtermdormancy_20151112_nocomments.csv` is not used.** It sits in the
source data folder but no script in `code/1.PopulationDynamics/` reads it —
it is a different experiment (the KBS long-term dormancy strains), so it is
not copied here.

**The end-of-experiment counts are not plotted.** `spore.transition.ltde.txt`
holds readings at ~1580 and ~1910 days. The original computed abundances from
them but never put them on the figure; they are carried through to
`data/pop_dynamics_end.csv` in case they should be.

**Re-running the state model needs JAGS**, which is not installed here, so
`R/90_latent_states_jags.R` exits with install instructions rather than
running. It writes to `output/latent_states_refit/` and compares against the
archive rather than overwriting it. It also loops over all four replicates —
the original was run one at a time by editing `trans.rep <- trans.A` and the
output filename by hand, which is how the C/D mix-up happened.

**Two tidy-ups to the model text**, neither of which changes the posterior:
`sig.v`, `sig.s` and their precisions were declared but never used by the
likelihood; and the two observation precisions were named the other way round
from the quantity they were applied to (`tau.to` on `Sobs`, `tau.so` on
`Tobs`), which only mislabelled the output.

**Which script was newest.** Of the four in `code/1.PopulationDynamics/`,
`1.PopDynamics_Log_2panel.Rmd` produces the two-panel figure and is the one
recreated here; the others are earlier single-panel and variance-stabilised
variants.

## Source

Recreated from `SporeMut/code/1.PopulationDynamics/1.PopDynamics_Log_2panel.Rmd`
with data from `data/1.PopDynamics/` and archived state estimates from
`output/1.PopDynamics/`.
