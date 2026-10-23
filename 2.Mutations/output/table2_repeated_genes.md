## Table 2. Genes hit by more than one independent mutation

Endpoint total fraction. Variants in the same gene carried by the same clone and no other are counted once, since they cannot be shown to be independent -- this is what removes *recX*, whose two variants are both in clone 4 alone and always together in its reads. Conversely, the two *ywcC* changes at 3,922,952 are different alleles in different clones and count as two events. Expected is the number of events a gene of that length would collect if all events were scattered across the coding genome at random. The family-wise column is a min-P permutation over all 4,171 CDS, which is the right correction when the genes were nominated by the data. Gene names and lengths are from the RefSeq annotation of NC_000964.3.

| Gene | Locus tag | Length (bp) | Independent events | Expected | P (Poisson) | P (family-wise) |
| :--- | :--- | ---: | ---: | ---: | :--- | :--- |
| sinR | BSU_24610 | 336 | 5 | 0.0025 | 7.53e-16 | <5e-05 |
| ywcC | BSU_38220 | 672 | 3 | 0.0049 | 1.98e-08 | <5e-05 |

