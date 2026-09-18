# SporeMut growth curves

Recreation of the growth-rate analysis from the SporeMut project: how the
maximum growth rate, yield, and lag time of ten evolved *Bacillus subtilis*
clones compare to their ancestor.

## Run it

```
Rscript run_all.R
```

About half a minute. Needs `dplyr`, `tidyr`, `ggplot2`, `ggridges`, `loo`, and
`lme4`. Everything lands in `output/`.

## Layout

```
data/
  comp_data_annotated.csv            analysis input: one row per growth curve,
                                     with the fitted Gompertz parameters, their
                                     SEs, and the plate run each curve came from
  comp_data.csv                      the original file, kept for reference
  treatments_original_corrected.csv  strain -> mutation / cell type
  plate_runs.csv                     which strain was run on which plate
  reference_jags_umax_relative.csv   posterior quantiles from the original run
  gompertz_fits/                     archived per-run curve fits (.txt) and
                                     diagnostic plots (.pdf)
  raw/                               Synergy MX plate-reader exports
bin/                                 the lab's Gompertz fitting machinery,
                                     unchanged from the original repo
R/
  00_setup.R                paths, packages, theme, strain metadata, brackets
  01_bayes_fitness.R        curve fits -> posterior draws
  02_figures.R              posterior draws -> figures
  03_validate.R             check against the original JAGS run
  04_group_models.R         do the strain groupings explain the parameters?
  05_report.R               manuscript tables and methods text
  90_gompertz_fits.R        raw plates -> curve fits    (slow, not in run_all)
  91_build_comp_data.R      curve fits -> analysis input
run_all.R
```

## Output

Three figures per parameter (`umax`, `yield`, `lag`), each annotated with the
nested group brackets — mutation (*sinR*, *slrC*), cell type (Spore, Total),
origin (Ancestor, Evolved):

| file | what it shows |
| --- | --- |
| `fig_<p>_relative_dist.pdf` | posterior density relative to the ancestor — the published panel |
| `fig_<p>_absolute_dist.pdf` | the same densities on the measured scale, ancestor included |
| `fig_<p>_mean_sd.pdf` | per-strain mean ± SD over the six curve fits, raw curves behind |

Strain names are deliberately off the figures: the panels are about group
membership, and the per-strain numbers are in Table 3 and the data files. Pass
`show_clones = TRUE` to `ridge_plot()` or `mean_sd_plot()` in
`R/02_figures.R` to put the labels back.

The brackets carry end caps so each group's extent is unambiguous. The
sporulation mutants are left unlabelled in the mutation row, as in the original
manuscript panel — the cell-type row already marks them as Spore. Tick marks at
strain positions are off, since the brackets identify the groups; the value
axes keep theirs.

Tables and text, as markdown ready to paste into a manuscript:

| file | what it is |
| --- | --- |
| `table1_group_models.md` | each grouping: between-strain SD remaining, and a likelihood-ratio test |
| `table2_group_ratios.md` | ratios of group means with 95% credible intervals |
| `table3_strain_estimates.md` | per-strain posterior estimates for all three parameters |
| `methods_stats.md` | the statistics paragraph for the methods section |
| `model_comparison.csv`, `model_contrasts.csv`, `posterior_summary.csv` | the same numbers unrounded, plus WAIC |

## The analysis

**Curve fitting.** Each well's OD600 trace is fit with a modified Gompertz
model, giving an intercept (`b0`), maximum yield (`A`), maximum growth rate
(`umax`), and lag time (`L`), each with a standard error. This is the original
Lennon lab code in `bin/`, driven by `R/90_gompertz_fits.R`.

**Curation.** Curves that failed visual QC were dropped, and the six
lowest-RMSE curves per strain kept. `R/91_build_comp_data.R` replays this.

**Bayesian estimation.** Each replicate carries its own uncertainty from the
curve fit, so replicates are combined by inverse-variance weight rather than
averaged. For each parameter and each strain *j*:

```
x[i,j] ~ LogNormal( mu[j], tau[j] * w[i,j] )     i = 1..6 curves
mu[j]  ~ Normal(0, precision 0.001)
tau[j] ~ Gamma(0.01, 0.01)
```

`w[i,j]` is the inverse-variance weight of curve *i* within strain *j*,
normalised to sum to 1. The log-normal keeps estimates positive during
sampling. `exp(mu[j])` is the strain's parameter; dividing by the ancestor's
draw gives relative fitness. Four chains, 10,000 draws each after burn-in.

**Group models** (`R/04_group_models.R`) ask whether the strain groupings
explain variation *between strains*:

```
y[i,j]   = log(parameter)
y[i,j]   ~ Normal( theta[j], sigma2_e / w[i,j] )   curve level
theta[j] ~ Normal( x[j]'beta, sigma2_s )           strain level
```

with `x[j]` a cell-means indicator for the grouping under test. Candidates:

| model | grouping | strains |
| --- | --- | --- |
| 1 | global mean | each scope |
| 2 | ancestor vs. evolved | all 11 |
| 3 | mutation vs. no mutation | — same partition as model 2, see below |
| 4 | spore vs. total | 10 evolved |
| 5 | *sinR* vs. *slrC* | 8 non-sporulation mutants |

Models 4 and 5 apply to subsets, so each is compared against its own
global-mean model within that subset; nothing is comparable across scopes.

Read `sigma_strain` (between-strain SD left after the grouping) and the group
ratios with their credible intervals. WAIC is reported but is close to blind
here: it scores prediction of one more curve from a strain already in the
model, by which point the strain effect has absorbed the group difference, and
`loo` flags it as unreliable with six curves per strain.

The Gibbs sampler and `lme4` fit the same model by different routes and agree
(ancestor/evolved ratio for umax: 1.466 vs 1.469).

## Things worth knowing

**The data files say `ywcC`; the figures and tables say `slrC`.** *ywcC* is the
legacy synonym for BSU_38220, whose current symbol is *slrC*. The mapping is
`GENE_SYNONYMS` in `R/00_setup.R` and is display-only — `treatments_original_corrected.csv`
and `comp_data.csv` are untouched, so the original name stays on record and the
archive checks still pass. Cite as **slrC (ywcC)** at first mention, since the
lab's records and the 2016-era literature use *ywcC*.

**Model 3 is model 2 here.** Every evolved isolate carries a mutation and the
ancestor carries none, so the two groupings are the same partition. They would
differ only with the unmutated spore isolates (S1, S6, S11, S22, S51, S95)
included; those are in `comp_data_annotated.csv` and are out of scope for the
growth-rate figures.

**The ancestor is one strain.** Its six curves are technical replicates, so in
model 2 the ancestor's group mean and that strain's own effect are the same
quantity, and the between-strain SD is informed only by the ten evolved clones.

**The `Curve` names.** These are raw well labels from the plate reader:
`M21_1` is strain M21, well 1. Each strain was run on two plates whose wells
are both numbered from 1, so a name like `M21_1` appears twice in the original
`comp_data.csv` with no way to tell the runs apart. `comp_data_annotated.csv`
adds `run`, `well`, and a unique `curve_id` (`M21-20230824-1`). The odd name
`S95_4.1` is R renaming a duplicated column header: `data/raw/20230824_...csv`
has two columns labelled `S95_4`, a plate-map typo most likely meant to be
`S95_3`.

**The sampler no longer needs JAGS.** The original ran the estimation model
through `rjags`. Because the model is conjugate, `R/01_bayes_fitness.R` also
carries an exact Gibbs sampler in ~15 lines of base R and defaults to it when
JAGS is not installed. `R/03_validate.R` compares against the archived JAGS
posterior: medians agree to 0.15%, the 95% interval bounds to 1.2%. To use
JAGS instead, install it plus `rjags` and the `"auto"` engine picks it up.

**Refitting reproduces the archive.** `R/90_gompertz_fits.R` writes to
`output/gompertz_refit/` rather than over `data/gompertz_fits/`, and compares
the two: 27 of 32 fit files come back identical. The five that differ do so
only on curves that were rejected during QC anyway (e.g. `M17_1`, a runaway fit
at umax 0.47 vs 0.50). Refitting renumbers rows, which would invalidate the
QC row positions recorded in `R/91_build_comp_data.R` — that is why the archive
is what the analysis reads.

**One inherited bug, fixed in the annotated file.** The original read
`S1_new.fit.parms.txt` where it meant `S6_new`, so the six S6 rows in
`comp_data.csv` are a second copy of S1. `comp_data_annotated.csv` reads the S6
fits that were there all along. The M clones are untouched, so the growth-rate
figures are unaffected; `R/91_build_comp_data.R` still verifies that the
original (buggy) file reproduces exactly from the archived fits.

**Time units.** The raw exports give time as `HH:MM:SS`; the original parses it
as `HH.MM`, so lag and growth rate are in those units rather than decimal
hours. Kept as is so the numbers match the manuscript.

**Everything is a script.** The original lived in a 1,800-line R Markdown file
(`code/4.GrowthCurves/test/GrowthCurve_CK.Rmd` in the SporeMut repo) that mixed
several parallel analyses — frequentist mixed models, `brms`, and the JAGS
version — with later chunks overwriting earlier variables. Only the JAGS
analysis reached the manuscript, so only that path is kept here.

**Updated for current R.** `stat(ecdf)` became `after_stat(ecdf)`, `size=`
became `linewidth=`, and `legend.position = c(.9, .8)` became
`legend.position = "inside"` plus `legend.position.inside`. Verified against
R 4.6.0 with ggplot2 4.0.3.

## Source

Recreated from `SporeMut/code/4.GrowthCurves/` — the `test/GrowthCurve_CK.Rmd`
chunks "Models with rjags" and "Yield and lag - Don's Bayes".
