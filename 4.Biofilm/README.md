# SporeMut biofilm assay

Recreation of the biofilm analysis from the SporeMut project: how much biofilm
ten evolved *Bacillus subtilis* clones produce relative to their ancestor.

## Run it

```
Rscript run_all.R
```

A few seconds. Needs `dplyr`, `tidyr`, `ggplot2`, `ggridges`, `loo`, and
`lme4`. Everything lands in `output/`.

## The plate layout, and where it disagrees

The layout for the 2020 reading is not in a file of its own — it is written as
notes underneath the grids in `OneDrive_3_7-15-2026/Biofilm/Behringer_Plate1.xlsx`
and `Behringer_Plate2.xlsx`. That is the only record of it anywhere in the
project folder, and it is transcribed to `data/plate_maps.csv`.

It disagrees with `data/biofil.csv` in three places:

| plate | column | Behringer note | biofil.csv |
| --- | ---: | :--- | :--- |
| 1 | 10 | **M66** | M79 |
| 1 | 11 | **B. subtilis 168 delta 6** | Ancestor |
| 2 | 9 | **S52** | S51 |

**M66 vs M79 is the one that matters.** That column holds the highest readings
on the plate (1.24–1.94), and under the `biofil.csv` labelling it is *m79* —
one of the three clones in the *slrC*/*epsA–slrR* group whose elevated biofilm
drives the group difference. `M66` appears **exactly once in the whole 470-file
folder**: in that note. `M79` appears throughout — in `treatments.csv`,
`comp_data.csv`, the growth-curve plates and the mutation matrices. On that
evidence the note is almost certainly a slip, and `biofil.csv` is right. But it
is the only independent record of the layout, so it should be confirmed rather
than assumed.

The same pattern covers **S52 vs S51**: `S51` appears 42 times in the project,
`S52` only in that note.

### 168 Δ6 vs "Ancestor" — the assay used the wrong ancestor

This one is not a slip. Δ6 is a specific domesticated 168 derivative, and the
recollection from the lab is that the original plate was run against the wrong
ancestor and the ancestor alone was re-read later. The 2023 file fits exactly:
its experiment is named **`SporeMutWT_Biofilm_230615`** — wild type only, no
strain panel, a lid comparison to settle the protocol.

**So everything expressed relative to the ancestor is provisional**, because
the denominator is probably the wrong strain. That is the relative ridge
figure, the `Relative biofilm` column of Table 3, and model 2
(ancestor vs. evolved). The scaling is a single divisor: if the true ancestor
sits a factor *k* away, every relative value moves by 1/*k*.

**What does not depend on it**, and can be used as it stands:

- the absolute posteriors and `fig_biofilm_absolute_dist.pdf`
- the lineage-group comparison, Table 4 and `fig_biofilm_lineage_groups.pdf` —
  *sinR* vs *slrC*/*epsA–slrR* vs spore never touches the ancestor
- models 4 and 5, which are fit to evolved clones only

`R/00_setup.R` prints this at every run and `ANCESTOR_IS_PROVISIONAL` there is
the flag to clear once the replacement is merged.

**Why the 2023 ancestor is not spliced in.** It is a separate run on a
different protocol — 540/600 nm against 550 nm, and with a lid — with no shared
reference well to calibrate between the two. Dividing 2020 clone readings by a
2023 ancestor reading would put a run-to-run difference straight into every
relative value. Merging it properly needs either a strain measured on both
plates or an ancestor re-read under the 2020 protocol.

## Why the 2023 file does not supersede this one

`OneDrive_3_7-15-2026/20230611/20230611_SporeMut_Biofilm.xlsx` is newer than
the 2020 reading, has no plate map, and looks at first like the assay this
project should be built on. It is not the same experiment. Four independent
signs, and no plate map exists to find:

- The Synergy export records its own experiment file as
  **`SporeMutWT_Biofilm_230615_105535.xpt`** — WT is in the name.
- Its two sheets are **`withLid` and `withoutLid`**: a methods comparison, not
  a strain comparison.
- The growth-curve file sitting beside it, from the same day, has a sheet named
  **"final-all samples are WT"**.
- The readings are far too uniform for a strain panel. Across 88 wells the 2023
  plate spans 0.083–0.503 (CV 52%); the 2020 strain plate spans 0.066–1.944
  (CV 134%), because it contains both m79 at 1.9 and m17 at 0.07. A 6-fold
  spread is ordinary well-to-well variation for one strain; a 30-fold spread is
  a panel.

So it is a wild-type lid test, and `data/biofil.csv` remains the biofilm
dataset. If a 2023 strain-panel reading exists it is somewhere else, and would
still need its layout.

## Biofilm by lineage group

An early manuscript draft splits the likely-vegetative clones by which biofilm
regulator they carry: *"clones containing a mutation in sinR did not also
contain mutations in ywcC or upstream of slrR"*. `R/06_lineage_groups.R` runs
that comparison on `data/biofil.csv`.

| group | strains | wells | mean OD₅₅₀ |
| --- | ---: | ---: | ---: |
| *sinR* | 5 | 40 | 0.042 ± 0.005 |
| spore | 2 | 16 | 0.144 ± 0.038 |
| *slrC*/*epsA–slrR* | 3 | 24 | 0.768 ± 0.118 |

*sinR* mutants form less biofilm, *slrC*/*epsA–slrR* mutants a great deal more
— P = 1.6 × 10⁻¹² between them, and 1.5 × 10⁻⁷ for *slrC*/*epsA–slrR* against
spores. **The *sinR*-versus-spore comparison is not significant** (P = 0.29),
and the reason is visible in the data: the two spore clones differ 16-fold from
each other, m23 at 0.272 against m26 at 0.017. Two strains with that spread
cannot support a comparison.

The draft's own P values are from an earlier round of the assay and are not
reproduced.

**On scope.** `data/biofil.csv` is plate 1 of the two-plate raw reading.
Plate 2 holds S1, S6, S11, S22, S51 and S95 in alternating columns with paired
media controls — a different layout, and not part of this analysis.
`R/06_lineage_groups.R` reads it only as a sensitivity check, because it is the
one thing that moves the weak comparison: adding those six to the spore group
takes *sinR* vs. spore from P = 0.29 to P = 5.9 × 10⁻⁴. The script also
verifies plate 1 of the raw file against `biofil.csv` (they agree to 10⁻¹⁶),
which pins the blank correction as the **median** — plate 1 has two
contaminated blank wells, 0.262 and 0.179 against a baseline near 0.064.

## Layout

```
data/
  biofil.csv    the assay: one row per well, with raw OD550, the plate blank,
                the blank-corrected reading, and each strain's treatment and
                mutation
R/
  00_setup.R            paths, packages, theme, strain metadata, brackets
  01_bayes_biofilm.R    assay -> posterior draws
  02_figures.R          posterior draws -> figures
  03_validate.R         check the sampler against the exact posterior
  04_group_models.R     do the strain groupings explain biofilm production?
  05_report.R           manuscript tables and methods text
run_all.R
```

## Output

Three figures, each annotated with nested group brackets — mutation (*sinR*,
*ywcC*), cell type (Spore, Vegetative), origin (Ancestor, Evolved):

| file | what it shows |
| --- | --- |
| `fig_biofilm_relative_dist.pdf` | posterior density relative to the ancestor — the published panel |
| `fig_biofilm_absolute_dist.pdf` | the same densities on the measured scale, ancestor included |
| `fig_biofilm_mean_sd.pdf` | per-strain geometric mean ± SD over the eight wells, wells behind |

Strain names are off the figures: the panels are about group membership, and
the per-strain numbers are in Table 3. Pass `show_clones = TRUE` to
`ridge_plot()` or `mean_sd_plot()` in `R/02_figures.R` to put them back.

Tables and text, as markdown ready to paste into a manuscript:

| file | what it is |
| --- | --- |
| `table1_group_models.md` | each grouping: between-strain SD remaining, and a likelihood-ratio test |
| `table2_group_ratios.md` | ratios of group means with 95% credible intervals |
| `table3_strain_estimates.md` | per-strain posterior estimates |
| `methods_stats.md` | the statistics paragraph for the methods section |
| `model_comparison.csv`, `model_contrasts.csv`, `posterior_summary.csv` | the same numbers unrounded, plus WAIC |

## The analysis

**Strain-level estimates.** Corrected OD550 spans two orders of magnitude
across strains, so it is modelled on the log scale. For strain *j*:

```
x[i,j] ~ LogNormal( mu[j], tau[j] )     i = 1..8 wells
mu[j]  ~ Normal(0, precision 0.001)
tau[j] ~ Gamma(0.01, 0.01)
```

Unlike the growth-curve assay there is no per-well standard error to weight by
— each well is a single reading — so replicates enter unweighted and each
strain gets its own variance. Four chains, 10,000 draws each after burn-in.
Relative biofilm is `exp(mu[j])` divided by the ancestor's, computed draw by
draw.

**Group models** (`R/04_group_models.R`) ask whether the strain groupings
explain variation *between strains*:

```
y[i,j]   = log(OD550 corrected)
y[i,j]   ~ Normal( theta[j], sigma2_e )    well level
theta[j] ~ Normal( x[j]'beta, sigma2_s )   strain level
```

| model | grouping | strains |
| --- | --- | --- |
| 1 | global mean | each scope |
| 2 | ancestor vs. evolved | all 11 |
| 3 | mutation vs. no mutation | — same partition as model 2 |
| 4 | spore vs. vegetative | 10 evolved |
| 5 | *sinR* vs. *slrC* | 8 non-sporulation mutants |

Read `sigma_strain` (between-strain SD left after the grouping) and the group
ratios. WAIC is reported but is close to blind here: with a strain effect in
the model it scores prediction of one more well from a strain already seen.

## Things worth knowing

**Everything is on a log axis.** Biofilm OD runs from 0.014 (m26) to 1.28
(m79), the model is fit on the log scale, and one strain's mean minus its SD
is negative. A linear axis would be unreadable and could not carry that error
bar. The mean ± SD panel therefore shows the geometric mean × / ÷ one
geometric standard deviation, with ticks labelled in OD units.

**Model 3 is model 2 here.** Every evolved isolate carries a mutation and the
ancestor carries none, so the two groupings are the same partition. Unlike the
growth-curve project there are no unmutated spore isolates in this assay, so
there is no version of the data that separates them.

**The ancestor is one strain.** Its eight wells are technical replicates, so
in model 2 the ancestor's group mean and that strain's own effect are the same
quantity, and the between-strain SD is informed only by the ten evolved clones.

**The data file says `ywcC/slrR`; the figures say `slrC`.** *ywcC* is the
legacy synonym for BSU_38220, whose current symbol is *slrC* — it sits beside
*slrA* and its product is the regulator of it. The mapping lives in
`MUTATION_LABEL` in `R/00_setup.R` and is display-only; `data/biofil.csv` is
untouched, so the original label stays on record. Cite as **slrC (ywcC)** at
first mention, since the lab's records and the 2016-era literature use *ywcC*.

**What `/slrR` means — resolved.** It is not the *slrR* gene. An early
manuscript draft names this group "*ywcC/epsA-slrR* mutants" and describes the
split as *"clones containing a mutation in sinR did not also contain mutations
in ywcC or upstream of slrR"*. "Upstream of *slrR*" is the intergenic
*epsA–slrR* variant at 3,529,981 — the marker carried by 47 of 84 sequenced
clones and excluded from Figure 2. So the draft's two subgroups are the two
clades the sequencing finds, and the mutual exclusivity was already known.

## `data/biofil.csv` is only half the experiment

`Biofilm_06_16_20.txt` holds **two plates**. Plate 1 is the eleven strains
`biofil.csv` carries; **plate 2 is S1, S6, S11, S22, S51 and S95** — the six
spore isolates — each with a paired blank. Only plate 1 was ever parsed.

That matters for the draft's comparison, which uses "known/likely spores":
the six S isolates plus m23 and m26. With plate 1 alone the spore group is
m23 and m26 only, whose means differ 16-fold (0.272 against 0.017), and
**sinR vs. spore is not significant at all** (P = 0.29). With both plates it is
(P = 5.9 × 10⁻⁴).

`R/06_lineage_groups.R` parses both plates and runs the three-group
comparison. It verifies plate 1 against `biofil.csv` first — they agree to
1 × 10⁻¹⁶, which also pins down the blank correction: the **median**, not the
mean, because plate 1 has two contaminated blank wells (0.262 and 0.179
against a baseline near 0.064).

| group | strains | wells | mean OD550 |
| --- | ---: | ---: | ---: |
| *sinR* | 5 | 40 | 0.041 ± 0.005 |
| known/likely spores | 8 | 64 | 0.227 ± 0.014 |
| *slrC*/*epsA–slrR* | 3 | 24 | 0.768 ± 0.118 |

All three pairwise differences are significant and in the directions the draft
reports — *sinR* mutants form less biofilm than spores, *slrC*/*epsA–slrR*
mutants more. **The exact P values do not reproduce**: the draft gives
6.30 × 10⁻⁶, 6.80 × 10⁻¹¹ and 4.7 × 10⁻⁵ where Holm-adjusted pooled-SD t-tests
on these data give 4.7 × 10⁻¹⁹, 5.9 × 10⁻⁴ and 4.3 × 10⁻¹⁴. Restricting the
spore group to the S isolates alone does not close the gap either. The
qualitative result is robust; the test the draft actually ran is not
recoverable from what survives.

**The sampler no longer needs JAGS.** The original ran this model through
`rjags`. Because it is conjugate, `R/01_bayes_biofilm.R` also carries an exact
Gibbs sampler in base R and defaults to it when JAGS is not installed. For a
single strain the precision can be integrated out analytically, leaving a
one-dimensional posterior that `R/03_validate.R` evaluates on a grid — so the
sampler is checked against the exact answer, not against another sampler.
Medians agree to 0.2%, interval bounds to 1.1%.

**Everything is a script.** The original lived in
`code/3.Biofilm/3.Biofilm_CK.Rmd`, which mixed `brms`, `lme4` and JAGS
analyses with later chunks overwriting earlier variables. Only the JAGS path
produced the archived figures, so only that path is kept here.

**Updated for current R.** `stat(ecdf)` became `after_stat(ecdf)`, `size=`
became `linewidth=`, and `legend.position = c(.85, .8)` became
`legend.position = "inside"` plus `legend.position.inside`.

**Shared with the growthCurves project.** `R/00_setup.R` (theme, brackets,
strain metadata) and the Gibbs samplers are near-duplicates of the ones in
`~/Desktop/growthCurves`. Kept self-contained so each project stands alone; if
more assays follow, they are worth factoring into a small package.

## Source

Recreated from `SporeMut/code/3.Biofilm/3.Biofilm_CK.Rmd` — the chunk
"Jags model - Don".
