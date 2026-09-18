## Populations persist under resource limitation

### 1. Energy limitation alters population dynamics

Energy limitation led to a large drop in total population size (*T*). There was a 99.70 ± 0.05% reduction over the first 10 days, a 357 ± 57-fold fall. Following this bottleneck, there was an extended period of equilibrium where total cell densities remained relatively unchanged for nearly 900 days, converging on a population size of 1.33 x 10^6^ ± 9.8 x 10^4^ cell/mL (Fig. 1, Table 2).

### 2. Energy limitation creates population heterogeneity

The population dynamics involved a large transition from mostly vegetative cell types to spores. In other words, energy limitation generated a shift in the equilibrium abundances of the two sub-populations. This is not entirely surprising, since the *Bacillus* stress response is known to let individuals form endospores under extreme energy limitation, and it is also known that not all individuals commit to sporulation synchronously. Despite endosporulation being understood as a responsive process, there appears to be some degree of bet-hedging, albeit consistent across populations. During the first 24 h, sporulation efficiency was low (0.107 ± 0.020%), but 81.5 ± 1.4% of remaining individuals were spores by day 10 (Fig. 1).

This was not due entirely to the death of vegetative cells: the absolute abundance of spores rose 12-fold over the experiment, from 8.59 x 10^4^ to 9.89 x 10^5^ CFU/mL, while non-spores fell 201-fold. Maximum-likelihood fitting of a sigmoidal function put the midpoint of this metabolic transition at 4.28 ± 0.58 d following the onset of energy limitation for the spore sub-population and 5.68 ± 0.21 d for the non-spore sub-population. By the end of the experiment, 76 ± 1% of cells were spores and 24 ± 1% were vegetative. A non-trivial fraction of a *B. subtilis* population can therefore withstand extended resource limitation as vegetative cells, whether actively dividing on resources turned over by dead cells and cannibalism, or in some non-spore form of dormancy.

---

### Where these numbers differ from the earlier draft

Every value above is generated from the analysis rather than typed in, so it moves if the data or the fits do. Four differences from the draft are worth knowing about.

- **Averages across replicates.** The original script used `mean(a, b, c, d)` to average the four replicate populations. R takes only the first argument as data and passes the rest to `trim` and `na.rm`, so those means were replicate A alone. The bottleneck is 99.70 ± 0.05% rather than 99.6%, and day-10 sporulation efficiency 81.5 ± 1.4% rather than 79.3%.
- **Replicate D.** The original read replicate C's state estimates in place of D's, so D never contributed its own spore and non-spore numbers. Fixed here; D's own estimates are used.
- **Uncertainty on the percentages.** The draft reports 99.6 +/- 0.0005% and 0.13 +/- 0.021%, mixing a percentage with the standard error of the underlying proportion. Both are on the same scale here.
- **Transition midpoint.** The draft gives 4.4 ± 2.27 d. The midpoint itself reproduces (4.28 ± 0.58 d for the spore sub-population), but the ± 2.27 could not be reproduced from these fits by any route tried -- the standard error across replicates is 0.58 d. Worth checking against whatever produced the draft before reusing it.

