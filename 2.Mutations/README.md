# SporeMut mutation accumulation

Recreation of Figure 2 and the mutation-accumulation results: which mutations
were found in each sequenced fraction of a 1000-day *Bacillus subtilis*
starvation experiment, and whether the variation needs selection to explain it.

## Run it

```
Rscript run_all.R
```

A few seconds. Needs `dplyr`, `tidyr`, `ggplot2`, `readxl` and `circlize`.

## The three fractions

| fraction | what it is | clones sequenced | source |
| --- | --- | --- | --- |
| S10 | early spore fraction | 52 | `NewControlBacillus.goodcoverage.xlsx` |
| S1000 | endpoint spore fraction | 96 | the GSF1925 run — `S1000_endpoint_spore.tsv` (no matrix survives) |
| S+NS1000 | endpoint total fraction | 84 | `ControlBacillus.filtered.xlsx` |

Each matrix has one row per variant and two blocks of one column per clone:
the called allele, then a 0/1 carrier indicator. The trailing `Total` column
is the row sum of the indicator block — the number of clones carrying that
mutation, which is what the rings are coloured by. **The column layouts differ
between the two workbooks**, so `R/01_mutations.R` names the positions per
file rather than reusing indices; it also checks that the carrier block sums
to `Total` and stops if it doesn't.

### Which file is the endpoint spore fraction

The endpoint spore fraction is **not** one of the three `*Bacillus*.xlsx`
matrices. It is the **GSF1925** sequencing run, which is why it was hard to
place: it arrived under the run identifier instead of a matching filename.
`Sequences/GSF1925-demultiplexing-summary.xlsx` lists 96 libraries,
`GSF1925-BS_Heatkilled_1` … `_96`. Heat-killing selects the spore fraction,
so this is S1000, n = 96 — the same 96 the results text reports four
mutations in.

The three `.xlsx` matrices resolve as:

| file | clones | fraction |
| --- | --- | --- |
| `ControlBacillus.filtered.xlsx` | 84 | S+NS1000, endpoint total — authoritative |
| `NewControlBacillus.goodcoverage.xlsx` | 52 | S10, early spores |
| `EvolvedBacillus.filtered.xlsx` | 65 | the failed first run — not used |

One coincidence is worth naming, because it invites a wrong conclusion. The
84 clone columns in `ControlBacillus.filtered.xlsx` are numbered 1–96 with
twelve gaps (8, 20, 57, 61, 65, 68, 69, 73, 82, 83, 90, 91), which makes it
look like 84 survivors of the 96 GSF1925 libraries — and would mean the total
fraction is really the spore fraction. It isn't. `ControlBacillus` was created
2016-10-19 and GSF1925 was demultiplexed 2018-10-10, two years later, so one
cannot be a subset of the other. Both sets were simply picked into 96-well
plates and numbered from the plate. Any subset of 1–96 is trivially "contained
in" 1–96, so the numbering carries no information about provenance.

The 96-clone variant matrix itself was never retained. What survives is a
summary table — `Simulations/STHK/STHK_Bootstrap_29.xlsx`, sheet `Sheet3` —
giving position, substitution, gene and product for each of the four
mutations at 1/96. That is the source of `data/S1000_endpoint_spore.tsv`.
Two of its four genes, *yetA*/*ylmG* and *levB*–*aspP*, are exactly the spore
clones S1 and S95 in the growth-curve `treatments.csv`, which is an
independent file — so the identification cross-checks.

> Sheet3's other rows are **not** usable. Its `LTS-M` block has denominators
> of 65, i.e. it was built from `EvolvedBacillus.filtered.xlsx`, the failed
> run. Only the `LTS-S` rows are taken from it.

### The Circos tracks are windows, not positions

`*.mutations.txt` and `*.circos.txt` are Circos *highlight* tracks: each row
is a ±10 kb window, so `start` is the mutation position **minus 10,000**.
Reading `start` as the position shifted all four endpoint spore calls 10 kb
down and put the *levB*–*aspP* call inside *epsB*, 860 bp from the excluded
*epsA*–*slrR* marker — a coincidence that looked like a second hit in the eps
operon and isn't. At its true position, 3,539,121, it is 9,140 bp away and in
no gene.

The convention is verified rather than assumed: `check_circos_offset()` in
`R/01_mutations.R` tests the two tracks whose variants also survive in a
matrix. All 34 total-fraction rows and both S10 rows match at +10,000 and
none at +0. If that ever stops holding, the run stops.

## Layout

```
data/
  ControlBacillus.filtered.xlsx          endpoint total fraction, 84 clones
  NewControlBacillus.goodcoverage.xlsx   early spore fraction, 52 clones
  S1000_endpoint_spore.tsv               endpoint spore fraction, 4 singletons
  S1000_endpoint_spore.circos.txt        the same four, as a Circos window track
  LT_NoTreat.mutations.txt               Circos window tracks kept only to
  ST_Heat.mutations.txt                    check the +/-10 kb convention
  EvolvedBacillus.filtered.xlsx          a separate 65-clone set (see below)
  mutations.csv                          built by R/01: the analysis table
R/
  00_setup.R        paths, genome constants, allele-frequency palette
  01_mutations.R    variant matrices -> tidy mutation table
  02_circos.R       Fig 2a: mutations around the genome
  03_null_model.R   neutral coalescent expectation and its sensitivity
  04_figures.R      Fig 2b: allele-frequency spectra
  06_lineage_structure.R  Fig 2c: co-occurrence and clone backgrounds
  07_parallelism.R  genes hit more often than chance allows
  05_report.R       tables and results text
run_all.R
```

## Output

| file | what it is |
| --- | --- |
| `fig2a_circos.pdf` / `.png` | three rings around the genome, coloured and sized by allele frequency |
| `fig2b_spectrum.pdf` / `.png` | observed vs neutral allele-frequency spectrum, one panel per fraction |
| `fig2c_lineages.pdf` / `.png` | clone-by-mutation matrix, clones grouped by genetic background |
| `results_text.md` | the four results paragraphs, numbers filled in from the analysis |
| `table1_mutations_by_fraction.md` | mutations, singletons and highest frequency per fraction |
| `table2_repeated_genes.md` | genes with more than one independent mutation |
| `tableS1_all_mutations.md` / `.csv` | every mutation, with position, gene and carrier count |
| `null_spectrum.csv`, `null_tests.csv`, `null_population_size_sweep.csv` | the neutral model's output |
| `clone_backgrounds.csv`, `cooccurrence_tests.csv` | which background each clone belongs to, and how the shared mutations nest |
| `mutation_load_by_clade.csv`, `clones_without_mutations.txt` | private mutations per clone by clade; the ten clones carrying nothing |
| `parallelism_tests.csv`, `clade_distribution_test.csv` | gene-level Poisson tests; how mutations split between the clades |
| `clone_mutation_counts.csv` | mutations per clone |

## Figure 2a

Rings run outside in — endpoint total, endpoint spore, early spore — as in the
published version. Two changes from the Circos original:

- **Colour and height both encode allele frequency**, on magma rather than the
  original rainbow, which had no perceptual order — a reader could not tell
  from the colours which of two tiles was the more frequent. The original
  mapping is kept in `FREQ_COLOURS_ORIGINAL`, since it is what decodes the
  archived track files. Carrier count also sets bar height, on a square-root
  scale so singletons stay visible next to a mutation in eight clones.
- **Each ring carries its own background tint**, repeated as a swatch beside
  its label in the centre, so there is no ambiguity about which circle is
  which fraction.
- **Gene labels come from the data**, not from a labels file: every gene that
  was hit more than once or is carried by more than one clone is labelled, so
  each tile drawn in a colour other than singleton blue can be read by name.

Tiles are widened to 20 kb so a single base is visible on a 4.2 Mb circle. That
span is cosmetic.

## The spore result needs no simulation

Spores experience no generations once the population contracts, so every
mutation they carry was acquired during the 29 generations of growth
beforehand. At 3.28 × 10⁻¹⁰ per base per generation that is 0.0401 mutations
per clone: **3.8 expected against 4 observed** in the 96 endpoint spore clones,
**2.1 against 2** in the 52 early spore clones. Poisson P = 1.0 for both.
`R/03_null_model.R` reports this before it runs anything.

That arithmetic is not falsifiable by a choice of population size, generation
count or population structure, which is what makes it the version to quote.
The coalescent reaches the same conclusion and is now supporting detail.

## Figure 2b and the neutral model

Under neutrality only the genealogy of the sampled clones matters, so
`R/03_null_model.R` traces the sampled lineages backwards through a
Wright–Fisher coalescent rather than forward-simulating the population — the
same model as the forward simulation it replaces, in base R, in seconds.
29 generations of growth from a single ancestor, then G generations at
constant size: G = 0 for spores, G = 1000 for non-spores. Mutation rate
3.28 × 10⁻¹⁰ per base per generation (Sung et al. 2015).

**This belongs in the supplement.** At 1000 stationary generations the model
predicts ~110 mutations in the total fraction where 34 were found. That is a
problem with the assumed number of generations, not with population structure,
and refitting with two demes would not touch it — the mutations also
distribute 17/17 across the two clades against 19 expected (P = 0.49), so the
data are not asking for a structured model. The spore claim is carried
analytically above and the selection claim by parallelism below; the site
frequency spectrum is then a caveat rather than a contradiction.

## Things worth knowing

**Two different sequencing sets are conflated in the draft.** The draft's
section 2 — "Of the 65 clones that were sequenced, 21 (32%) had acquired a
mutation… 16 had only one mutation, two clones had three mutations, and three
clones had five" — matches `EvolvedBacillus.filtered.xlsx`: 65 clones, 22
(34%) with a mutation, and 17, 2 and 3 clones with one, three and five. But
that file shares **no positions and no genes** with the 84-clone
`ControlBacillus.filtered.xlsx` (*pbpD*, *cspB*, *yacL*… against *sinR*,
*recX*, *ywcC*). Per the project lead, the 84-clone run is the correct one and
the 65-clone file is a failed pooling; this recreation uses the 84-clone matrix
throughout, which is also what Fig. 2a is drawn from and what the draft's own
section 4 describes. The 65-clone file is kept in `data/` for reference only.

**"There were also 34 singleton mutations" should be 24.** There are 34
mutations in the total fraction, 10 of them in more than one clone.

**The neutral-model conclusion turns on one unmeasured number.** The draft
reports that mutations in the total fraction are at frequencies "far greater
than neutrally expected (P < 2.2 × 10⁻¹⁶)". That holds only for a large
non-spore population. At N = 2.9 × 10⁵ — the equilibrium non-spore density in
CFU per mL — 1000 generations of drift produce about 7.4 shared mutations by
chance against the 10 observed, which is not a surprise (P = 0.21). The
conclusion needs N of order 10⁷ or more, where the expectation drops to 0.4.
N is a count of cells, not a density, so what is missing is the culture
volume. At 2.9 × 10⁵ non-spores per mL, roughly 25 mL or more puts N past 10⁷
and the conclusion holds. `output/null_population_size_sweep.csv` now carries
the implied volume for each N — this is one measurement away from settled and
worth settling rather than publishing conditionally.

**Choosing 1000 generations is conservative in the opposite direction from the
draft's reasoning.** The draft picked it "to ensure that any deviation from the
expected distribution was due to positive selection and not an underestimation
of generations". More generations means both more mutations and more drift, so
it makes neutral sharing *more* likely, not less — which is why the test is
hard to pass at small N.

**One variant is excluded from the counts, but not from the analysis.**
Position 3,529,981 (*epsA-slrR*) is carried by 47 of the 84 total-fraction
clones (56%), far too common to have arisen during the experiment. It is
excluded from Fig. 2a and the mutation counts, as it was from the published
figure — but it is kept in `R/06_lineage_structure.R`, because what it marks
turns out to matter more than the variant itself.

**The total fraction is two lineages, not one population.** No clone carrying
a *sinR* mutation carries *epsA-slrR*, and none of the 47 that carry
*epsA-slrR* has a *sinR* mutation. Treating the 20 *sinR* clones as
independent makes that look decisive; they are not independent, they descend
from five mutational events. At the event level the probability that all five
landed among the 44% of clones lacking the marker is 0.44⁵ = **0.017** —
suggestive, not decisive. Both numbers are reported.

**But the clade structure does not distort the counts.** Of the 34 mutations
other than the marker, 17 fall inside the *epsA-slrR* clade and 17 outside,
against 19 expected from the clade's size (binomial P = 0.49). And neither
clade mutates faster (0.28 / 0.20 / 0.41 private mutations per clone, Poisson
GLM P = 0.49). So clonal structure is *not* inflating allele frequencies, and
the earlier worry that it might does not survive checking.

Some mutations do nest — the two *comEC* clones are a strict subset of the
eight carrying *sinR* 2,552,884, and three of the four *ypeB* clones carry
*sinR* 2,552,905 — so a variant in four clones can still be one event carried
up by an expanding lineage. No variant spans the two clades, which sounds
meaningful but is not: any mutation arising after the split necessarily lands
on one side.

**Gene-level parallelism is the cleanest evidence, and it needs no
simulation.** Of the 34 mutations, 27 are coding; collapsing variants that
share a gene and a carrier set leaves **26 independent events**. Scattered
across the 3,684,498 bp of coding sequence, *sinR* (336 bp) expects 0.0024 of
them and has five. Family-wise P < 5 × 10⁻⁵ by min-P permutation over all
4,171 CDS. *ywcC* (*slrC*, 672 bp) has two genuine events, P = 0.009. None of
this depends on allele frequency, population size, or the sample being well
mixed — unlike both the frequency statement and the neutral simulation the
manuscript currently leans on. `R/07_parallelism.R`.

**Three corrections that changed the answer:**

- **Events, not variant rows.** *recX*'s two variants are both in clone 4 and
  no other — a 49 bp deletion at 925,845 and a 1 bp insertion at 925,897, 52 bp
  apart. That is one mutational event at best and an alignment artifact at the
  deletion breakpoint at worst. Counted as two it looked significant; counted
  correctly it drops out entirely. **The draft's "three genes that received
  multiple independent mutations" is two.** Clone 4 also carries a *ywcC*
  indel, so one clone accounts for three of the eight indels in this fraction —
  worth a coverage check if the reads surface.
- **Exact CDS lengths, mapped by position.** `data/NC_000964.3_cds.tsv` is
  extracted from the RefSeq GFF for the accession already in the Methods.
  Mutations are assigned to whatever CDS contains them rather than matched by
  name. Size matters — *pksN* is 16 kb, so two hits there are unremarkable.
  Three legacy names are genuine renames: *ywcC* → **slrC** (BSU_38220),
  *yllA* → **bshC**, *ysiB* → **fadB**.
- **One RefSeq symbol is wrong, and is overridden.** BSU_22920 comes back from
  the GFF as *sleC*, but *B. subtilis* has no *sleC* — its cortex-lytic
  enzymes are SleB and CwlJ. UniProt P38490 for that locus is **YpeB**, and
  *sleB* (BSU_22930) sits immediately downstream: this is the bicistronic
  *ypeB–sleB* operon. `data/gene_name_overrides.csv` resolves names by
  locus_tag and keeps *ypeB*, with the reasoning recorded. This matters —
  under "sleC" the most interesting non-*sinR* hit would have read as noise.
- **Coding mutations as the denominator.** An intergenic variant cannot land in
  a CDS, so counting the 7 intergenic ones would dilute the per-gene rate.

**On the correction itself.** Bonferroni over every CDS is reported but is not
the number quoted: the genes were nominated by the data, so the right null is a
permutation. Two versions are computed. The crude one — does *any* gene collect
this many events? — is dominated by large genes and is *harsher* on a small
gene than Bonferroni (0.13 vs 0.047 for *ywcC*), which is not what a family-wise
correction should do. The reported one is min-P, which judges each gene against
its own length, and gives *ywcC* P = 0.009.

**The two blocks are not the same kind of object**, and Fig. 2c says so.
*epsA-slrR* is one clade — 47 clones sharing a single ancestral variant, and
remarkably uniform: 26 of the 47 carry it and nothing else, with four small
sub-lineages confined to one edge. That is what a standing variant looks like,
expanded early and accumulated little since. The *sinR* block is five
independent mutational events that share only the absence of that variant, so
it is drawn one shade per allele rather than as a single group.

The third group is labelled **no clade marker**, not "private variants only":
clone 88's sole mutation is *ypeB*, shared by four clones, and clones 4 and 38
carry *ywcC* 3,922,952, shared by two. What they have in common is the absence
of a clade marker, not privacy.

**Neither clade mutates faster.** Counting only private variants, which carry
no signal of expansion, gives 0.28 mutations per clone in the *epsA-slrR*
clade, 0.20 in the *sinR* clones and 0.41 in the rest (Poisson GLM P = 0.49).
That rules out the obvious alternative to selection — that one background
simply accumulates mutations more quickly — and it is a strong negative result
that the frequency rings cannot show. The row is on Fig. 2c.

**What the panel supports saying**: two alternative routes to altered biofilm
regulation, one standing variant that expanded and recurrent de novo mutation
of *sinR*, that never co-occur. That is what would be expected if either alone
achieves the phenotype, and it is a stronger claim than the bare observation
that they are mutually exclusive.

**Naming is now consistent across the projects.** *slrC* (*ywcC*) is the
display name in `3.GrowthCurves` and `4.Biofilm` too, mapped in each project's
`00_setup.R`. No data file was changed — `treatments_original_corrected.csv`,
`comp_data.csv` and `biofil.csv` all still say *ywcC* (or *ywcC/slrR*), so the
archive checks still pass and the original labels stay on record. Convention is
**slrC (ywcC)** at first mention, since the lab's records and the 2016-era
literature use *ywcC*.

**Constants that differ from the draft.** The equilibrium population is
1.33 × 10⁶ here, computed from the counts, against 1.36 × 10⁶ in the draft;
and the mutation rate is 3.28 × 10⁻¹⁰ throughout — the figure in the draft's
results text — while its Methods say 3.2 × 10⁻¹⁰. The difference is 2.5% and
changes nothing, but one of each should be chosen.

**One endpoint spore call is 860 bp from the excluded marker.** Position
3,529,121 falls in *epsB*; the excluded variant at 3,529,981 is intergenic
between *epsA* and *slrR*. Two calls that close in a repeat-rich operon, one of
them already treated as suspect, deserve an alignment check.
`output/variants_near_excluded_marker.csv`.

**Ten clones carry no mutation at all**: 23, 26, 32, 43, 46, 52, 70, 74, 87
and 96. M23 and M26, the two clones representing the spore group in Figure 3,
are two of ten equally qualified candidates — nothing in the sequencing
distinguishes them from the other eight. Whoever picked them had a reason, and
it is not recorded in any file recovered here. The list is in
`output/clones_without_mutations.txt`.

**The 17-generation result reproduces exactly.** For one mutant lineage to
reach the observed 9.5% of a population of 1.33 × 10⁶ cells by division alone
takes log2(0.095 × 1.33 × 10⁶) = 17 generations. That population size comes
from the `populationDynamics` project's equilibrium total.

**Not used.** `113_All.variants.txt`, `changecolors.sh`, the ticks files and
`SampleFiles/` in the Circos folder are inputs to the original Circos run,
which `circlize` does not read. `B_sub_NoTreat.conf` and `B_sub_LT_Heat.conf`
hold the original ring spec and were used as reference for radii and order.

## Source

Recreated from `SporeMut/OneDrive_3_7-15-2026/SporeMut_Data_MB/Circos_Plot_Figure_2/`
(published figure and ring inputs), `OneDrive_3_7-15-2026/Bacillus_Mutations/`
(variant matrices), and `SporeMut/try/figure1/` (the coalescent null model this
one is built on).
