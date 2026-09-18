# Manuscript figures

Two figures, assembled from the four analysis projects.

```bash
Rscript run_all.R
```

Run the four projects first — this one reads their saved output and does not
recompute anything.

| | |
| --- | --- |
| `output/figure1_dynamics_mutations.pdf` | what the population did, and what changed in its genomes |
| `output/figure2_phenotypes.pdf` | phenotypes of the evolved clones, by mutation |
| `output/figureS1`…`figureS4` | supplementary figures |
| `output/figure1_caption_notes.md` | Figure 1's caption, including the parallelism numbers |
| `output/figure2_caption_notes.md` | Figure 2's caption, including the ancestor caveat |
| `output/supplementary_index.md` | what each supplementary figure is and where it comes from |

## Figure 1

| panel | what | from |
| --- | --- | --- |
| A | spore and non-spore abundance over 1000 days, with the counted total over them | `1.PopDynamics/R/03_figures.R` |
| B | mutations around the genome, one ring per fraction | `2.Mutations/R/02_circos.R` |

The link to Figure 2 is on the circos itself: genes taking more independent
events than the genome-wide rate allows are named with their event count and
set in bold — **sinR (5)**, **slrC (2)**. Those are the two genes phenotyped in
Figure 2.

This was a third panel at one point. It is two numbers, and two numbers do not
need a panel, so the statistics moved to `output/figure1_caption_notes.md`:
*sinR* 5 events against 0.0024 expected, 2100-fold, P < 5e-05; *slrC* 2 against
0.0047, 420-fold, P = 0.0086. The P values are min-P permutation family-wise
values, the right correction when the gene was nominated by the data rather
than chosen in advance. `N_PERM` is read out of
`2.Mutations/R/07_parallelism.R`, so the bound quoted for *sinR* cannot drift
from the number of permutations actually run.

Because the circos marks those genes, `2.Mutations/run_all.R` runs
`07_parallelism.R` before `02_circos.R`. Nothing else depended on that order.

Panel A's total is the counted S + V, drawn over the modelled spore/non-spore
split and carried in the same legend as a red line rather than as a floating
label.

## Figure 2

| panel | what | from |
| --- | --- | --- |
| A | maximum growth rate | `3.GrowthCurves/R/02_figures.R` |
| B | maximum yield | `3.GrowthCurves/R/02_figures.R` |
| C | lag time | `3.GrowthCurves/R/02_figures.R` |
| D | biofilm | `4.Biofilm/R/02_figures.R` |

Group brackets are drawn once per row, on the left-hand panel, and apply to
both panels in that row — all four share one strain order.

### Why the measured scale, not relative to the ancestor

`PHENOTYPE_SCALE` in `R/00_setup.R` is `"absolute"`. Both scales exist in the
projects, and the growth parameters would be fine either way, but the biofilm
ancestor is not: `4.Biofilm/data/biofil.csv` calls the reference well
"Ancestor" where the original plate record names it *B. subtilis* 168 Δ6,
which is probably not this experiment's ancestor. Every relative biofilm value
inherits that, while the absolute posteriors do not. Showing all four panels
absolute keeps them on one footing and keeps the figure clear of the problem.

Set `PHENOTYPE_SCALE <- "relative"` to switch; nothing else needs changing.

Each panel carries a dashed line at the ancestor's posterior mean — the same
quantity the relative panels divide by, drawn rather than divided out.

`SHOW_ANCESTOR_FOOTNOTE` in `R/00_setup.R` is `FALSE`, so the ancestor caveat
is **not** printed on the figure. It is still written to
`output/figure2_caption_notes.md` every run, which says there that it is off
the figure. Set it `TRUE` to print it on the biofilm panel; it would go on that
panel rather than on the shared brackets, since the brackets are read by the
lag-time panel beside it, whose ancestor is not in doubt.

## Supplementary figures

`R/03_supplementary.R` collects them so that everything going to the manuscript
is in one folder. Nothing is redrawn: each is the plot object its own project
builds, pulled out through `load_project()` and re-saved under a supplementary
name. `output/supplementary_index.md` lists them with a sentence each.

| | figure | from |
| --- | --- | --- |
| S1 | allele-frequency spectra against the neutral expectation | `2.Mutations/R/04_figures.R` |
| S2 | clone genotypes and the lineage structure | `2.Mutations/R/06_lineage_structure.R` |
| S3 | biofilm by lineage group | `4.Biofilm/R/06_lineage_groups.R` |
| S4 | the 2023 wild-type biofilm re-read (diagnostic) | `4.Biofilm/R/07_ancestor_2023.R` |

To add one, append an entry to `SUPPLEMENTARY` in `R/03_supplementary.R` giving
the project, the script, and the name that script gives its finished plot. A
missing object is a warning and a skip, not a failure, so a half-run project
does not stop the rest.

## How the panels get here

Nothing about a panel is restated in this project. `load_project()` in
`R/00_setup.R` sources the script that owns the panel, with the saving calls
(`ggsave`, `pdf`, `png`, `dev.off`) stubbed out, and hands back the
environment holding the objects. The manuscript figure and the project figure
are therefore the same code, and a change to a project's ordering, brackets or
axis rule reaches the manuscript figure the next time this is run.

`source` is stubbed too. The project scripts locate and source their own
`00_setup.R` with the default `local = FALSE`, which evaluates it in the
global environment; two projects loaded in one session would then overwrite
each other's `PLOT_ORDER` and `DATA_DIR`. Redirecting it keeps each project's
constants with that project.

What this project decides, and the only things to change here:

- layout, panel lettering, figure dimensions
- `PHENOTYPE_SCALE` — which scale Figure 2 uses
- `SHOW_ANCESTOR_FOOTNOTE` — whether the ancestor caveat is printed on Figure 2
- `TOTAL_COLOUR` — the counted total in Figure 1A
- which figures go in the supplement (`SUPPLEMENTARY`)
- `CELL_GREYS` — the two greys in Figure 1A are pushed further apart than the
  standalone panel uses, because the panel is smaller here and the two cell
  types are told apart by grey level alone
- brackets once per row rather than on every panel
- the ancestor footnote

## Requires

`ggplot2`, `patchwork`, `cowplot`, `gridGraphics`, plus whatever the four
projects need. `cowplot::as_grob` replays the circos, which `circlize` draws
with base graphics, onto a grid device — so panel B stays vector rather than
being rasterised into the figure.
