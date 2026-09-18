# Mutation accumulation in heterogenous population of dormant bacteria

Four replicate populations of *Bacillus subtilis* were maintained for 1000 days
without added resources. This repository holds the analysis behind the
manuscript: population dynamics, whole-genome sequencing of clones from three
sampled fractions, and phenotyping of the evolved clones.

Each numbered folder is a self-contained project that reproduces its own outputs
from a single command. They run in order — later ones read earlier ones' output —
but each can be re-run on its own.

```bash
cd 1.PopDynamics && Rscript run_all.R
```

| | what it does | main output |
| --- | --- | --- |
| [`1.PopDynamics`](1.PopDynamics) | total, spore and vegetative abundance over 898 days | `fig_spore_vegetative.pdf`, population descriptors |
| [`2.Mutations`](2.Mutations) | variant matrices → mutation table, circos figure, neutral model, parallelism tests | `fig2a_circos.pdf`, `tableS1_all_mutations.md` |
| [`3.GrowthCurves`](3.GrowthCurves) | Gompertz fits → Bayesian strain estimates for µmax, yield, lag | posterior summaries, group models |
| [`4.Biofilm`](4.Biofilm) | crystal-violet assay → Bayesian strain estimates | posterior summaries, lineage-group comparison |
| [`5.MsFigures`](5.MsFigures) | assembles Figures 1–2 and the supplementary figures | `figure1_*.pdf`, `figure2_*.pdf`, `figureS1`–`S4` |

Every project has its own README with the detail: what the data are, what was
decided and why, and what is still uncertain.

The manuscript itself is drafted in Overleaf and is not kept here. Its LaTeX
tables are generated from these projects rather than transcribed, so they have
to be rebuilt and re-uploaded whenever an analysis changes.

## What the analysis shows

Total abundance fell 357-fold over the first ten days, then held near
1.33 × 10⁶ CFU/mL for the next 900 days. The population became mostly but not
entirely spores — 76% at the end, leaving a quarter of cells in a non-spore
state.

Spores accumulated nothing. Four mutations among 96 endpoint spore clones and
two among 52 early spore clones, in both cases what the 29 generations of growth
*before* resource limitation already predict, and every one a singleton.

The whole-population sample was different: 48 of 84 clones carried at least one
of 34 mutations. Two genes took more independent mutations than the genome-wide
rate allows — *sinR* (5 events, P < 5 × 10⁻⁵) and *slrC* (2 events, P = 0.009),
both biofilm regulators. That test is the one the argument rests on: it needs no
coalescent, no assumption about population size, and no assumption that the
sample is well mixed.

## Design

The four analysis projects share a structure: `data/` (inputs, unmodified),
`R/` (numbered scripts), `output/` (everything generated), `run_all.R`.

Nothing generated is edited by hand. `5.MsFigures` builds its panels by loading
the analysis projects' own figure code rather than reimplementing it, and the
manuscript's tables are generated from these projects' own table files — so a
number in the manuscript cannot differ from the number the analysis produced.

Re-running a project overwrites its `output/`. Large posterior draws
(`output/*.rds`) are not tracked; they regenerate in seconds.

## Data provenance

This is a reconstruction of an older project, and several things were resolved
along the way that a reader should know about. Each is documented where it
matters, in the relevant project's README:

- **The endpoint spore fraction** is the GSF1925 sequencing run (96 heat-killed
  libraries), not one of the ambiguously named `*Bacillus*.xlsx` matrices. Its
  full variant matrix was not retained; the four mutations survive in a summary
  table. See `2.Mutations/README.md`.
- **The Circos tracks are ±10 kb windows, not point positions.** Reading their
  start column as a position shifts every call 10 kb. This is asserted at every
  run against the two tracks whose variants also survive in a matrix.
- **BSU_22920 is *ypeB*, not *sleC*** — the RefSeq symbol is wrong, and gene
  names are resolved by locus tag with documented overrides in
  `2.Mutations/data/gene_name_overrides.csv`.
- **The biofilm reference well is provisional.** It is labelled "Ancestor" in the
  assay file, but the original plate record names it *B. subtilis* 168 Δ6, which
  is probably not this experiment's ancestor. Absolute values and every
  comparison among the evolved clones are unaffected; anything relative to the
  ancestor is not. See `4.Biofilm/README.md`.
- **The neutral-simulation result is conditional** on the assumed size of the
  non-spore population, which the data do not pin down. It is reported as
  supplementary for that reason.

## Requires

R (≥ 4.4) with `ggplot2`, `ggridges`, `patchwork`, `cowplot`, `gridGraphics`,
`dplyr`, `tidyr`, `readxl`, `circlize`, `bbmle`, `lme4`, `loo`, `viridisLite`.
The population-dynamics latent-state model uses JAGS via `rjags`; everything
else is base R or the packages above.

## Citation

<!-- Add the manuscript citation and a DOI for this repository once available. -->
