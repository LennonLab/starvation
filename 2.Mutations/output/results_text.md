## Mutation accumulation under energy limitation

### 1. Mutations do not accumulate in spores

The developmental steps of creating a spore, including genome replication, could lead to mutations. Of the 52 clones from the S~10~ treatment that were whole-genome sequenced, we detected only 2 mutations (3.8%). Without active repair, mutations could also accumulate during the 1000-day experiment even with no turnover or cryptic replication, and other mechanisms such as mobile genetic elements might contribute to genomic change in inactive spores. We see little evidence for either. Of the 96 clones from the endpoint spore population, we detected only 4 mutations (3.1% of clones carried one). Every mutation in both spore fractions was a singleton.

That is exactly the number expected from the growth phase alone, and it takes no simulation to see. Spores experience no generations once the population contracts, so every mutation they carry was acquired during the 29 generations of growth beforehand. At a mutation rate of 3.28 x 10^-10^ per base per generation that is 0.0401 mutations per clone: 3.8 expected in the 96 endpoint spore clones against 4 observed, and 2.1 expected in the 52 early spore clones against 2 observed. Both land inside the 95% Poisson interval of the expectation (1-8 and 0-5 respectively). Nothing accumulated while they were dormant.

### 2. Mutations accumulate in the total fraction

After 1000 days of resource limitation the total fraction had acquired far more. Of the 84 clones sequenced, 48 (57%) carried at least one mutation, against 2 of 52 (4%) in the early spore fraction. In total 34 distinct mutations were identified, appearing 63 times across clones, 10 of them in more than one clone (Fig. 2a, Table S1). Of the clones carrying mutations, 36 had 1, 10 had 2, 1 had 3, 1 had 4.

### 3. The neutral simulation, and why it is now supplementary

We used a neutral coalescent to ask whether the variation in the total fraction needs selection to explain it. Cells were assumed to experience log2(10^9^) = 29 generations of growth before resource limitation, with a *B. subtilis* mutation rate of 3.28 x 10^-10^ per base per generation (Sung et al. 2015). Spores were given no further generations; the non-spore subpopulation was given 1000, one per day, a deliberate upper bound so that any excess could not be blamed on undercounting generations. Across 100 iterations, the neutral expectation for the two spore fractions is essentially all singletons, which is what was observed (Fig. 2b) -- the same conclusion section 1 reaches by arithmetic, and the arithmetic is the version to quote.

For the total fraction the answer depends on the size of the non-spore population, and this is the one place where the analysis does not support the earlier draft. At N = 2.9 x 10^5^ -- the equilibrium non-spore density in CFU/mL -- the neutral model expects 7.8 mutations in more than one clone, against the 10 observed (P = 0.26): sharing at these frequencies is what drift alone produces over 1000 generations. The conclusion only holds once N is around 10^7^ or larger, where the expectation falls to 0.38 shared mutations and the observed count is highly improbable (P = 1.23e-11). Whether N should be cells per mL or cells in the whole culture therefore decides the result; see `output/null_population_size_sweep.csv`.

Independent of that, for a single mutant lineage to reach the observed 9.5% of a population of 1.33 x 10^6^ cells by division alone requires at least 17 generations over the 1000 days. The non-spore cells are therefore likely vegetative, and the equilibrium in population size during extreme resource limitation is consistent with balanced cryptic growth.

### 4. The total fraction is two lineages, not one population

The variant at *epsA-slrR* is carried by 47 of the 84 sequenced clones. It is far too common to have arisen during the experiment and the figures exclude it, but what it marks matters: no clone carrying a *sinR* mutation carries it (Fig. 2c). Treating the 20 *sinR* clones as independent gives a vanishingly small probability, but they are not independent -- they descend from five mutational events. At the level of those events the question is whether all five landed among the 44% of clones lacking the clade marker, which has probability 0.017. Suggestive; not decisive.

The clade structure does not, however, distort the mutation counts. Of the 34 mutations other than the clade marker, 17 fall inside the *epsA-slrR* clade and 17 outside it, against 19 expected from the clade's size (binomial P = 0.49). Mutations are spread between the two clades as evenly as chance predicts, so clonal structure is not inflating allele frequencies. Nor does either clade mutate faster: counting only private variants gives 0.28 mutations per clone in the *epsA-slrR* clade, 0.20 in the *sinR* clones and 0.41 in the rest (Poisson GLM P = 0.49).

The two clades are different kinds of object, and Fig. 2c draws them that way. *epsA-slrR* is one clade: 47 clones sharing a single ancestral variant, and remarkably uniform -- 26 of them carry it and nothing else, with four small sub-lineages (*ybdN*, *pksN*, *citZ-ytwI*, *yutK*) confined to one edge. That is what a standing variant looks like: expanded early, accumulated little since. The *sinR* block is the opposite -- five separate events, each with a small set of descendants, with *comEC* nested inside one and *ypeB* inside another. No variant spans the two, which sounds meaningful but is not: any mutation arising after the split necessarily lands on one side.

### 5. Genes hit more often than chance allows

The strongest evidence of selection here comes from neither allele frequency nor the neutral simulation. Of the 34 mutations in the total fraction, 27 fall in a coding sequence; collapsing variants that share a gene and a carrier set leaves 26 independent events. Scattered at random across the 3,684,498 bp of coding sequence, a gene the size of *sinR* (336 bp) expects 0.0024 of them. It has five, in five alleles carried by non-overlapping sets of clones. A min-P permutation over all 4,171 coding sequences puts the family-wise probability below 5 x 10^-5^.

*ywcC* -- *slrC* in the current annotation -- is hit twice, at 3,922,748 in one clone and 3,922,952 in two others, giving a family-wise P = 0.009. *recX* does **not** belong on this list, though the draft names it: both of its variants, a 49 bp deletion at 925,845 and a 1 bp insertion at 925,897, are carried by clone 4 and no other clone. Fifty-two bases apart in a single clone, they are one mutational event at best and an alignment artifact around the deletion breakpoint at worst. Two genes show parallelism, not three.

This argument needs no coalescent, no assumption about population size and no assumption that the sample is well mixed. It is the line the manuscript should lead with; the frequency statement and the simulation comparison both rest on assumptions the data strain.

### 6. What were the mutations?

Two genes took more than one independent mutation in the total fraction: *sinR* (5 independent events) and *slrC* (2 independent events). Beyond repeated hits in a single gene, 10 of the 34 mutations in the total fraction were present in more than one clone, and one mutation in *sinR* was found in 8 (9.5%) of the sequenced clones. The remaining 24 were singletons.

---

### Where these numbers differ from the earlier draft

- **Section 2 describes a different dataset.** The draft's "65 clones, 21 (32%) had acquired a mutation, 16 had one, two had three, three had five" matches `EvolvedBacillus.filtered.xlsx` (65 clones; 22 = 34% with a mutation; 17, 2 and 3 clones with one, three and five). That file shares **no positions and no genes** with `ControlBacillus.filtered.xlsx`, which is the 84-clone matrix behind Fig. 2a and behind the draft's own section 4 (*sinR* x5, *recX* x2, *ywcC* x2, sinR in 8 = 9.5%). The two sections describe different sequencing sets. This recreation uses the 84-clone matrix throughout, for consistency with the figure.

- **"There were also 34 singleton mutations"** should be 24. There are 34 mutations in total, 10 of them in more than one clone.

- **The neutral-model conclusion is conditional**, as set out above. The draft reports P < 2.2 x 10^-16^ without stating N, and the model assumes a well-mixed population that the co-occurrence structure rules out.

- **The draft names three genes with multiple independent mutations; it is two.** *recX*'s two variants are both in clone 4 alone, 52 bp apart, one of them a 49 bp deletion -- one event, or an artifact at the deletion breakpoint. Clone 4 also carries a *ywcC* indel, so a single clone accounts for three of the eight indels reported in this fraction; worth a coverage check if the reads can be found.

- **The omission of the other high-frequency variants is defensible but unstated.** The draft's criterion appears to be genes with multiple independent mutations, which only *sinR* and *ywcC* meet once *recX* is corrected. That is the right criterion given the lineage structure -- but *ypeB* (a germination protein, missense, in four clones) and *comEC* (a frameshift in a competence gene) are worth a sentence in a paper about the boundary between sporulating and vegetative cells.

- **One variant is excluded.** The variant at position 3,529,981 is carried by 47 of the 84 total-fraction clones (56%). It is far too common to have arisen during the experiment and was excluded from the published figure; it is excluded here too, which is what leaves 34 mutations.

- **Two constants differ from the draft.** The equilibrium population is 1.33 x 10^6^ here, computed from the counts, against 1.36 x 10^6^ in the draft; and the mutation rate is 3.28 x 10^-10^ throughout, which is the figure in the draft's own results text, while its Methods say 3.2 x 10^-10^. The 2.5% difference changes nothing, but one should be chosen.

- **The endpoint spore calls were 10 kb out, and are fixed.** They used to be read from `S1000_endpoint_spore.circos.txt`, whose rows are ±10 kb Circos highlight windows rather than point positions, so every position was 10,000 bp low. The correct values come from the surviving summary of the 96-clone GSF1925 run (see `data/S1000_endpoint_spore.tsv`): *yetA* 777,008; *ylmG* 1,611,379; *kdgA* 2,323,251; and *levB*–*aspP* 3,539,121. The offset is verified against the two tracks whose variants also survive in a matrix — all 34 total-fraction rows and both S10 rows match at +10,000 and none at +0 — and is asserted at every run. This retires an earlier caveat: the fourth call appeared to sit 860 bp from the excluded *epsA*–*slrR* marker and inside *epsB*, which looked like a second eps-operon hit worth an alignment check. At its true position it is 9,140 bp away and in no gene, so the proximity was an artefact of the offset and there is nothing to check.

- **Ten clones carry no mutation at all**: 23, 26, 32, 43, 46, 52, 70, 74, 87, 96. M23 and M26, which represent the spore group in Figure 3, are two of ten equally qualified candidates; nothing in the sequencing distinguishes them from the other eight. Whoever chose them had a reason, and it is not recorded in any of the files recovered here.

