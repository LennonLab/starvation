# Supplementary figures

Built by `R/03_supplementary.R`. Each is the figure its own project draws;
this project only collects them and gives them supplementary numbering.

| | figure | source | file |
| --- | --- | --- | --- |
| **S1** | Allele-frequency spectra against the neutral expectation | `2.Mutations/R/04_figures.R` | `figureS1_allele_frequency_spectra.pdf` |
| **S2** | Clone genotypes and the lineage structure behind them | `2.Mutations/R/06_lineage_structure.R` | `figureS2_lineage_structure.pdf` |
| **S3** | Biofilm by lineage group | `4.Biofilm/R/06_lineage_groups.R` | `figureS3_biofilm_by_lineage_group.pdf` |
| **S4** | Diagnostic: the 2023 wild-type biofilm re-read | `4.Biofilm/R/07_ancestor_2023.R` | `figureS4_ancestor_2023_diagnostic.pdf` |

## Notes

**S1. Allele-frequency spectra against the neutral expectation.** Observed spectra for the three sequenced fractions against the neutral coalescent. The two spore fractions are all singletons, which is what neutrality predicts; the endpoint total fraction is not.

**S2. Clone genotypes and the lineage structure behind them.** Which clone carries which mutation, ordered so the clades are visible. The apparent mutual exclusivity of sinR and epsA-slrR is a property of the lineages, not of the mutations: they descend from five events, not twenty independent observations.

**S3. Biofilm by lineage group.** The three-group comparison: sinR against slrC/epsA-slrR against spore. This never uses the ancestor, so it is unaffected by the ancestor problem described for Figure 2.

**S4. Diagnostic: the 2023 wild-type biofilm re-read.** Where the 2023 wild-type reading falls against the 2020 plate. The two runs read different wavelengths against different blanks and share no strain, so this is not a calibration and nothing in the analysis depends on it. Include only if the ancestor question is discussed.

