## Methods: statistics

### Growth curves

Optical density at 600 nm was recorded every 15 min for each well on a microplate reader across 4 plate runs. Each trace was fit with a modified Gompertz model, yielding a lower asymptote, the maximum yield (A), the maximum specific growth rate (µmax) and the lag time (L), each with an asymptotic standard error. Curves that failed visual inspection of the fit diagnostics were discarded and the six lowest-RMSE curves per strain retained, giving 66 curves across 11 strains (the ancestor and ten evolved clones). The six curves of a strain are technical replicates, not independent isolates.

### Strain-level estimates

Because the fitted parameters differ in precision from curve to curve, replicates were combined by inverse-variance weighting in a Bayesian model rather than averaged. For strain *j*, the log of the fitted value of curve *i* was modelled as Normal(µ[j], σ[j] / √w[i,j]), where w[i,j] is the inverse of the squared standard error of that curve's fit, scaled to mean one within a strain, and each strain has its own σ[j]. The model was fit in brms with Stan: µ[j] ~ Normal(0, 31.6) and log σ[j] ~ Normal(0, 5), the latter effectively flat over any plausible range. Four chains of 3,500 iterations with 1,000 warm-up gave 10,000 draws; all R-hat values were below 1.01. Estimates are posterior medians with 95% equal-tailed credible intervals, and are what the figures show. Values relative to the ancestor are exp(µ[j]) divided by the ancestor's exp(µ), computed draw by draw so that the interval carries the uncertainty in both terms.

### Group comparisons

Whether the strain groupings explain variation between strains was tested with linear mixed models (lme4): the log-transformed parameter with one fixed mean per group, a random intercept for strain, and curves weighted by the inverse variance of their fit. Candidate groupings were a single global mean; ancestor vs. evolved; mutation vs. no mutation; spore vs. total; and mutation identity (*sinR*, *ywcC* (BSU_38220), or none). Groupings that apply to a subset of the strains were compared against a global-mean model fit to that same subset, and no comparison is made across subsets.

Groupings were compared two ways. As a multimodel comparison, by AICc (MuMIn) from maximum-likelihood fits, reported as Akaike weights within each set of strains. And as a test, by a Kenward-Roger F-test against the global-mean model on REML fits (pbkrtest), which gives the reported *P*. Kenward-Roger rather than a chi-square likelihood-ratio test because there are only 8 to 11 strains, and the chi-square test is anticonservative with so few groups. Pairwise ratios of group means, with 95% confidence intervals and Kenward-Roger *P* values, were computed with emmeans, Tukey-adjusted where a grouping has three levels.

Two features of the strain set limit what these comparisons can show. The ancestor is a single strain, so in the ancestor vs. evolved comparison its group mean and its own effect are the same quantity. And among the ten evolved clones, the two without a detected mutation (M23 and M26) are also the two from the spore fraction, so mutation status and cell type are the same partition and cannot be separated.

### Software

Analyses were run in R version 4.6.0 (2026-04-24) with brms 2.23.0, rstan 2.32.7, lme4 2.0.1, pbkrtest 0.5.5, MuMIn 1.48.19, emmeans 2.0.3 and ggplot2 4.0.3. The brms strain-level model reproduces a reference JAGS fit's posterior medians to within 0.25%; its 95% intervals are narrower than the reference's because the reference's Gamma(0.01, 0.01) prior on the precision was not vague at the scale the weights put it on. Code and data are in the growthCurves project.

