## Methods: statistics

### Population sampling

Four replicate populations (A–D) were sampled at 463 timepoints over 898 days. At each timepoint, total and spore abundances were measured as colony-forming units: total counts from a plated dilution series, spore counts from the same sample after treatment that leaves only spores viable. Counts were converted to CFU per mL by multiplying by the dilution factor and by ten, the plated volume being 100 µL. Sampling paused between days 554 and 588.

### Vegetative abundance

Vegetative cells were not counted directly: they are the difference between the total and the spore count, and because both counts carry plating error that difference is non-positive at 174 of 1844 replicate-timepoints. Vegetative and spore abundances were therefore treated as latent states inferred jointly. For timepoint *i*, the spore state S[i] and vegetative state V[i] were given gamma priors, and the observations modelled as Sobs[i] ~ Normal(S[i], tau_1) and Tobs[i] ~ Normal(V[i] + S[i], tau_2), so that the counted total constrains the sum and the counted spores constrain one component. Gamma shape and rate were set by moment-matching the observed abundances, separately for the transition phase (before day 10) and the equilibrium phase after it, since the two differ by orders of magnitude. Observation standard deviations were given uniform(0, 1000) priors. The model was fit per replicate in JAGS with five chains, 1,000 burn-in iterations and 5,000 samples; posterior means are used as the state estimates.

### Abundance trajectories

Each series was described with a sigmoidal curve on log-log axes, log10(N) = b + (a − b) / (1 + exp((m − log10(t)) / w)), where *b* is the level the population starts from, *a* the level it settles at, *m* the midpoint in log10 days and *w* the width of the transition. Curves were fit by maximum likelihood assuming normal error on log10(N), to the four replicates pooled and to each replicate separately. Day-zero readings were retained: log10(0) sends the curve to exactly *b*, so they inform the starting asymptote and nothing else.

### Software

Analyses were run in R version 4.6.0 (2026-04-24) with bbmle 1.0.26, ggplot2 4.0.3 and patchwork 1.3.2. The latent-state model was fit in JAGS via rjags. Code and data are in the populationDynamics project.

