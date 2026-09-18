## Methods: statistics

### Biofilm assay

Biofilm production was measured as crystal-violet retention read at 550 nm in a microtitre plate, with 8 replicate wells per strain and 88 wells in total across 11 strains (the ancestor and ten evolved clones). Each reading was corrected by subtracting the mean of the plate's blank wells. The evolved clones carrying sporulation mutations were assayed as spores and the remainder as vegetative cells. Replicate wells are technical replicates of a single strain, not independent isolates.

### Strain-level estimates

Corrected optical densities span two orders of magnitude across strains, so they were modelled on the log scale. For strain *j*, well *i* was modelled as x[i,j] ~ LogNormal(µ[j], τ[j]), each strain having its own variance; the log-normal keeps estimates positive during sampling. Priors were µ[j] ~ Normal(0, precision 0.001) and τ[j] ~ Gamma(0.01, 0.01). Because every full conditional is conjugate, the posterior was sampled with a Gibbs sampler written in base R: four chains of 10,000 draws each after 2,000 burn-in iterations. Estimates are reported as posterior medians with 95% equal-tailed credible intervals. Relative biofilm is exp(µ[j]) divided by the ancestor's exp(µ), computed draw by draw so that the interval carries the uncertainty in both terms, and is plotted on a log10 axis.

### Group comparisons

To ask whether the strain groupings explain variation between strains, the log-transformed readings were fit with a two-level model: y[i,j] ~ Normal(θ[j], σ²e) at the well level and θ[j] ~ Normal(x[j]'β, σ²s) at the strain level, where x[j] is a cell-means indicator for the grouping under test, so each element of β is one group's mean on the log scale. Priors were β ~ Normal(0, 100) and inverse-Gamma(0.01, 0.01) on both variances. Candidate groupings were: a single global mean; ancestor vs. evolved; spore vs. vegetative; and mutation identity (*sinR*, *slrC* -- labelled *ywcC/slrR* in the assay file -- sporulation mutant). Mutation vs. no mutation is the same partition as ancestor vs. evolved on this strain set and was not fit separately. Groupings that apply to a subset of the strains were compared against a global-mean model fit to that same subset, and no comparison is made across subsets.

Support for a grouping is reported as the reduction in the between-strain standard deviation σs relative to the global-mean model, together with a likelihood-ratio test against that model fit by maximum likelihood (lme4) and the posterior ratio of group means with its credible interval. WAIC was also computed but is reported only in the supplementary output: with a strain-level effect in the model it scores prediction of a further well from a strain already observed, which is not the question being asked. In the ancestor vs. evolved comparison the ancestor group contains a single strain, so its group mean and that strain's own effect are the same quantity and σs is informed only by the ten evolved clones; that contrast should be read with the ancestor's lack of biological replication in mind.

### Software

Analyses were run in R version 4.6.0 (2026-04-24) with ggplot2 4.0.3, ggridges 0.5.7, lme4 2.0.1 and loo 2.9.0. The Gibbs samplers are conjugate implementations in base R. For the strain-level model the nuisance precision can be integrated out analytically, and the sampler reproduces the resulting exact marginal posterior to within 0.2% on medians and 1.1% on the bounds of the 95% credible intervals; the group-level sampler agrees with the corresponding lme4 fit. Code and data are in the biofilm project.

