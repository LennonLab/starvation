## Methods: statistics

### Biofilm assay

Biofilm production was measured as crystal-violet retention read at 550 nm in a microtitre plate, with 8 replicate wells per strain and 88 wells in total across 11 strains (the ancestor and ten evolved clones). Each reading was corrected by subtracting the mean of the plate's blank wells. The two spore-fraction clones were assayed as spores and the remainder as vegetative cells. Replicate wells are technical replicates of a single strain, not independent isolates.

### Strain-level estimates

Corrected optical densities span two orders of magnitude across strains, so they were modelled on the log scale: for strain *j*, the log reading of well *i* was Normal(µ[j], σ[j]), each strain with its own SD and no weighting, the wells being equally precise. The model was fit in brms with Stan, with µ[j] ~ Normal(0, 31.6) and log σ[j] ~ Normal(0, 5). Four chains of 11,000 iterations with 1,000 warm-up gave 40,000 draws; all R-hat values were below 1.01. Estimates are posterior medians with 95% equal-tailed credible intervals, and are what the figures show.

### Group comparisons

Whether the strain groupings explain variation between strains was tested with linear mixed models (lme4): log OD with one fixed mean per group and a random intercept for strain. Candidate groupings were a single global mean; ancestor vs. evolved; mutation vs. no mutation; spore vs. total; and mutation identity (*sinR*, *ywcC* (BSU_38220), or none). Groupings that apply to a subset of strains were compared against a global-mean model on the same subset, and no comparison is made across subsets. Groupings were compared by AICc (MuMIn) from maximum-likelihood fits, reported as Akaike weights, and tested by a Kenward-Roger F-test against the global-mean model on REML fits (pbkrtest), which gives the reported *P*. Pairwise ratios of group means with 95% confidence intervals and Kenward-Roger *P* values were computed with emmeans, Tukey-adjusted where a grouping has three levels. Throughout, the strain is the unit of replication: the eight wells of a strain are technical replicates and are not treated as independent.

The reference well is recorded as "Ancestor" in the assay file but named *B. subtilis* 168 Δ6 in the plate record, and is probably not this experiment's ancestor. Every comparison that involves it --- ancestor vs. evolved, and mutation vs. no mutation on all eleven strains --- is provisional. Comparisons among the evolved clones do not use it.

### Software

Analyses were run in R version 4.6.0 (2026-04-24) with brms 2.23.0, rstan 2.32.7, lme4 2.0.1, pbkrtest 0.5.5, MuMIn 1.48.19, emmeans 2.0.3 and ggplot2 4.0.3. For the strain-level model the posterior of µ[j] can be computed exactly on a grid under the same priors, and the brms fit reproduces it to within 1% on medians and 3% on the bounds of the 95% intervals. Code and data are in the biofilm project.

