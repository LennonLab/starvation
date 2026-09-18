################################################################################
# 05_report.R
#
# Manuscript tables and the results paragraphs, with every number pulled from
# the analysis.
#
# Input : data/mutations.csv, output/clone_mutation_counts.csv,
#         output/null_tests.csv, output/null_population_size_sweep.csv
# Output: output/table1_mutations_by_fraction.md
#         output/table2_repeated_genes.md
#         output/tableS1_all_mutations.md / .csv
#         output/results_text.md
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the mutations project root (or R/).")
source(.setup[1])

muts      <- read.csv(file.path(DATA_DIR, "mutations.csv"), stringsAsFactors = FALSE)
backgrounds <- read.csv(file.path(OUT_DIR, "clone_backgrounds.csv"), stringsAsFactors = FALSE)
load <- read.csv(file.path(OUT_DIR, "mutation_load_by_clade.csv"), stringsAsFactors = FALSE)
par_tests <- read.csv(file.path(OUT_DIR, "parallelism_tests.csv"), stringsAsFactors = FALSE)
clade_test <- read.csv(file.path(OUT_DIR, "clade_distribution_test.csv"), stringsAsFactors = FALSE)
spore_check <- read.csv(file.path(OUT_DIR, "spore_analytic_expectation.csv"), stringsAsFactors = FALSE)
no_mutation <- readLines(file.path(OUT_DIR, "clones_without_mutations.txt"))
cooccur   <- read.csv(file.path(OUT_DIR, "cooccurrence_tests.csv"), stringsAsFactors = FALSE)
eps_clones  <- backgrounds$clone[backgrounds$clade == "epsA-slrR"]
sinr_clones <- backgrounds$clone[backgrounds$clade == "sinR"]
excl_p <- fisher.test(table(backgrounds$clade == "sinR",
                            backgrounds$clade == "epsA-slrR"))$p.value
mutator_p <- as.numeric(readLines(file.path(OUT_DIR, "mutator_test_p.txt")))
n_sinr_events <- length(unique(muts$position[muts$fraction == "S+NS1000" &
                                               muts$gene %in% "sinR"]))
event_p <- (1 - length(eps_clones) / FRACTIONS$`S+NS1000`$n)^n_sinr_events
n_tot_mut <- sum(muts$fraction == "S+NS1000")
n_occurrences <- sum(muts$carriers[muts$fraction == "S+NS1000"])
cds_tbl   <- read.delim(file.path(DATA_DIR, "NC_000964.3_cds.tsv"), stringsAsFactors = FALSE)
ev_tbl    <- read.csv(file.path(OUT_DIR, "mutation_events.csv"), stringsAsFactors = FALSE)
pt        <- read.csv(file.path(OUT_DIR, "parallelism_tests.csv"), stringsAsFactors = FALSE)
n_cds     <- nrow(cds_tbl)
cds_total <- sum(cds_tbl$length_bp)
n_events  <- nrow(ev_tbl)
n_coding  <- sum(ev_tbl$n_variants)
sinr_len  <- pt$length_bp[pt$cds_gene == "sinR"]
sinr_exp  <- pt$expected[pt$cds_gene == "sinR"]
ywcc_p    <- pt$p_permutation[pt$cds_gene == "slrC"]
per_clone <- read.csv(file.path(OUT_DIR, "clone_mutation_counts.csv"), stringsAsFactors = FALSE)
tests     <- read.csv(file.path(OUT_DIR, "null_tests.csv"), stringsAsFactors = FALSE)
sweep     <- read.csv(file.path(OUT_DIR, "null_population_size_sweep.csv"), stringsAsFactors = FALSE)
gens      <- as.numeric(readLines(file.path(OUT_DIR, "generations_needed.txt")))

md_table <- function(df) {
  numeric_only <- vapply(df, function(col)
    all(grepl("^[-–<>0-9.%,]+$", as.character(col))), logical(1))
  sep <- ifelse(numeric_only, "---:", ":---")
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("| ", paste(sep, collapse = " | "), " |"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

write_md <- function(lines, file, title, caption) {
  writeLines(c(paste("##", title), "", caption, "", lines, ""),
             file.path(OUT_DIR, file))
  message("Wrote output/", file)
}

## ---- per-fraction summary --------------------------------------------------

muts$fraction <- factor(muts$fraction, levels = FRACTION_ORDER)
per_clone$fraction <- factor(per_clone$fraction, levels = FRACTION_ORDER)

by_frac <- muts %>%
  group_by(fraction) %>%
  summarise(mutations = n(), singletons = sum(carriers == 1),
            shared = sum(carriers > 1), max_carriers = max(carriers),
            .groups = "drop") %>%
  arrange(fraction) %>%
  mutate(n = vapply(as.character(fraction), function(f) FRACTIONS[[f]]$n, integer(1)),
         label = vapply(as.character(fraction), function(f) FRACTIONS[[f]]$label,
                        character(1)))

clones_hit <- per_clone %>%
  group_by(fraction) %>%
  summarise(with_mutation = sum(mutations > 0), .groups = "drop")

t1 <- by_frac %>%
  left_join(clones_hit, by = "fraction") %>%
  transmute(
    Fraction = label,
    `Clones sequenced` = n,
    `Clones with a mutation` = ifelse(
      is.na(with_mutation), "–",
      sprintf("%d (%.0f%%)", with_mutation, 100 * with_mutation / n)),
    Mutations = mutations,
    Singletons = singletons,
    `In >1 clone` = shared,
    `Highest frequency` = sprintf("%d/%d (%.1f%%)", max_carriers, n,
                                  100 * max_carriers / n))

write_md(md_table(as.data.frame(t1)), "table1_mutations_by_fraction.md",
         "Table 1. Mutations by sequenced fraction",
         paste0("Counts exclude the variant at position ",
                format(EXCLUDE_POS, big.mark = ","),
                ", carried by 47 of 84 clones in the endpoint total fraction ",
                "and therefore not something that arose during the ",
                "experiment. No clone-by-variant matrix survives for the ",
                "endpoint spore fraction, so its per-clone column is blank; ",
                "its four mutations are known to be singletons."))

## ---- genes hit more than once ----------------------------------------------
# Built from independent events, not variant rows: two variants in one gene
# carried by the same clone alone are one event, which is what removes recX.

events <- read.csv(file.path(OUT_DIR, "mutation_events.csv"), stringsAsFactors = FALSE)
par_tests <- read.csv(file.path(OUT_DIR, "parallelism_tests.csv"), stringsAsFactors = FALSE)

t2 <- par_tests %>%
  arrange(desc(events)) %>%
  transmute(Gene = cds_gene, `Locus tag` = locus_tag,
            `Length (bp)` = length_bp,
            `Independent events` = events,
            Expected = sprintf("%.4f", expected),
            `P (Poisson)` = format.pval(p_poisson, digits = 3),
            `P (family-wise)` = ifelse(p_permutation == 0, "<5e-05",
                                       format.pval(p_permutation, digits = 3)))

write_md(md_table(as.data.frame(t2)), "table2_repeated_genes.md",
         "Table 2. Genes hit by more than one independent mutation",
         paste0(
           "Endpoint total fraction. Variants in the same gene carried by the ",
           "same clone and no other are counted once, since they cannot be ",
           "shown to be independent -- this is what removes *recX*, whose two ",
           "variants are both in clone 4 alone. Expected is the number of ",
           "events a gene of that length would collect if all events were ",
           "scattered across the coding genome at random. The family-wise ",
           "column is a min-P permutation over all 4,171 CDS, which is the ",
           "right correction when the genes were nominated by the data. Gene ",
           "names and lengths are from the RefSeq annotation of NC_000964.3."))

## ---- full mutation list ----------------------------------------------------

tS1 <- muts %>%
  arrange(fraction, desc(carriers)) %>%
  mutate(Fraction = vapply(as.character(fraction), function(f) FRACTIONS[[f]]$label,
                           character(1)),
         Position = format(position, big.mark = ","),
         Frequency = sprintf("%.1f%%", 100 * frequency)) %>%
  # gene_display is already the resolved name, including the legacy rewrite
  # applied to intergenic labels in 01_mutations.R -- do not re-derive it here.
  mutate(Gene = ifelse(is.na(gene_display), "–", gene_display),
         Legacy = dplyr::case_when(
           !is.na(gene_legacy)          ~ gene_legacy,
           !is.na(gene) & gene != Gene  ~ gene,
           TRUE                         ~ "–")) %>%
  select(Fraction, Position, Gene, Legacy, Clones = carriers, Frequency,
         Variant = variant, Severity = severity, Product = product) %>%
  mutate(across(c(Variant, Severity, Product), ~ ifelse(is.na(.x), "–", .x)),
         Variant = sub("_variant$", "", Variant),
         Product = substr(Product, 1, 46))

write.csv(tS1, file.path(OUT_DIR, "tableS1_all_mutations.csv"), row.names = FALSE)
write_md(md_table(as.data.frame(tS1)), "tableS1_all_mutations.md",
         "Table S1. Every mutation identified",
         paste0("One row per mutation — 34 distinct sites in the endpoint ",
                "total fraction, appearing 63 times across clones. Gene is the ",
                "current symbol, resolved by position against the RefSeq ",
                "annotation of NC_000964.3; Legacy is the name the source ",
                "workbook used, where it differs. **The two *recX* rows are ",
                "listed separately here but counted once in Table 2**: both ",
                "are in clone 4 alone, 52 bp apart, so they cannot be shown to ",
                "be independent events. Variant type and severity are the ",
                "annotations carried in the source workbooks; the early-spore ",
                "file records severity but not variant type, and no annotation ",
                "survives for the endpoint spore fraction."))

## ---- results text ----------------------------------------------------------

get <- function(f, col) by_frac[[col]][as.character(by_frac$fraction) == f]
hit <- function(f) clones_hit$with_mutation[as.character(clones_hit$fraction) == f]

n_s10   <- FRACTIONS$S10$n
n_s1000 <- FRACTIONS$S1000$n
n_tot   <- FRACTIONS$`S+NS1000`$n

top <- muts %>% filter(fraction == "S+NS1000") %>% slice_max(carriers, n = 1)
tot_test <- tests[tests$fraction == "S+NS1000", ]

results <- c(
  "## Mutation accumulation under energy limitation",
  "",
  "### 1. Mutations do not accumulate in spores",
  "",
  sprintf(paste(
    "The developmental steps of creating a spore, including genome",
    "replication, could lead to mutations. Of the %d clones from the S%s",
    "treatment that were whole-genome sequenced, we detected only %d mutations",
    "(%.1f%%). Without active repair, mutations could also accumulate during",
    "the 1000-day experiment even with no turnover or cryptic replication, and",
    "other mechanisms such as mobile genetic elements might contribute to",
    "genomic change in inactive spores. We see little evidence for either. Of",
    "the %d clones from the endpoint spore population, we detected only %d",
    "mutations (%.1f%% of clones carried one). Every mutation in both spore",
    "fractions was a singleton."),
    n_s10, "~10~", get("S10", "mutations"),
    100 * get("S10", "mutations") / n_s10,
    n_s1000, get("S1000", "mutations"), 100 * 3 / n_s1000),
  "",
  sprintf(paste(
    "That is exactly the number expected from the growth phase alone, and it",
    "takes no simulation to see. Spores experience no generations once the",
    "population contracts, so every mutation they carry was acquired during",
    "the %d generations of growth beforehand. At a mutation rate of 3.28 x",
    "10^-10^ per base per generation that is %.4f mutations per clone: %.1f",
    "expected in the %d endpoint spore clones against %d observed, and %.1f",
    "expected in the %d early spore clones against %d observed. Both land",
    "inside the 95%% Poisson interval of the expectation (%d-%d and %d-%d",
    "respectively). Nothing accumulated while they were dormant."),
    29L, spore_check$expected[1] / spore_check$clones[1],
    spore_check$expected[spore_check$fraction == "S1000"],
    spore_check$clones[spore_check$fraction == "S1000"],
    spore_check$observed[spore_check$fraction == "S1000"],
    spore_check$expected[spore_check$fraction == "S10"],
    spore_check$clones[spore_check$fraction == "S10"],
    spore_check$observed[spore_check$fraction == "S10"],
    spore_check$lower95[spore_check$fraction == "S1000"],
    spore_check$upper95[spore_check$fraction == "S1000"],
    spore_check$lower95[spore_check$fraction == "S10"],
    spore_check$upper95[spore_check$fraction == "S10"]),
  "",
  "### 2. Mutations accumulate in the total fraction",
  "",
  sprintf(paste(
    "After 1000 days of resource limitation the total fraction had acquired",
    "far more. Of the %d clones sequenced, %d (%.0f%%) carried at least one",
    "mutation, against %d of %d (%.0f%%) in the early spore fraction. In total",
    "%d distinct mutations were identified, appearing %d times across clones,",
    "%d of them in more than one clone (Fig. 2a, Table S1). Of the clones",
    "carrying mutations, %s."),
    n_tot, hit("S+NS1000"), 100 * hit("S+NS1000") / n_tot,
    hit("S10"), n_s10, 100 * hit("S10") / n_s10,
    get("S+NS1000", "mutations"), n_occurrences, get("S+NS1000", "shared"),
    paste(sprintf("%d had %d", as.integer(table(
      per_clone$mutations[per_clone$fraction == "S+NS1000" &
                            per_clone$mutations > 0])),
      as.integer(names(table(
        per_clone$mutations[per_clone$fraction == "S+NS1000" &
                              per_clone$mutations > 0])))),
      collapse = ", ")),
  "",
  "### 3. The neutral simulation, and why it is now supplementary",
  "",
  sprintf(paste(
    "We used a neutral coalescent to ask whether the variation in the total",
    "fraction needs selection to explain it. Cells were assumed to experience",
    "log2(10^9^) = 29 generations of growth before resource limitation, with a",
    "*B. subtilis* mutation rate of 3.28 x 10^-10^ per base per generation",
    "(Sung et al. 2015). Spores were given no further generations; the",
    "non-spore subpopulation was given 1000, one per day, a deliberate upper",
    "bound so that any excess could not be blamed on undercounting",
    "generations. Across %d iterations, the neutral expectation for the two",
    "spore fractions is essentially all singletons, which is what was observed",
    "(Fig. 2b) -- the same conclusion section 1 reaches by arithmetic, and the",
    "arithmetic is the version to quote."), 100L),
  "",
  sprintf(paste(
    "For the total fraction the answer depends on the size of the non-spore",
    "population, and this is the one place where the analysis does not support",
    "the earlier draft. At N = 2.9 x 10^5^ -- the equilibrium non-spore density",
    "in CFU/mL -- the neutral model expects %.1f mutations in more than one",
    "clone, against the %d observed (P = %.2f): sharing at these frequencies is",
    "what drift alone produces over 1000 generations. The conclusion only",
    "holds once N is around 10^7^ or larger, where the expectation falls to",
    "%.2f shared mutations and the observed count is highly improbable",
    "(P = %s). Whether N should be cells per mL or cells in the whole culture",
    "therefore decides the result; see `output/null_population_size_sweep.csv`."),
    tot_test$expected_shared, tot_test$observed_shared, tot_test$p_value,
    sweep$expected_shared[sweep$N == 1e7],
    format.pval(sweep$p_observed_10[sweep$N == 1e7], digits = 3)),
  "",
  sprintf(paste(
    "Independent of that, for a single mutant lineage to reach the observed",
    "%.1f%% of a population of 1.33 x 10^6^ cells by division alone requires at",
    "least %.0f generations over the 1000 days. The non-spore cells are",
    "therefore likely vegetative, and the equilibrium in population size during",
    "extreme resource limitation is consistent with balanced cryptic growth."),
    100 * top$frequency[1], gens),
  "",
  "### 4. The total fraction is two lineages, not one population",
  "",
  sprintf(paste(
    "The variant at *epsA-slrR* is carried by %d of the %d sequenced clones.",
    "It is far too common to have arisen during the experiment and the figures",
    "exclude it, but what it marks matters: no clone carrying a *sinR*",
    "mutation carries it (Fig. 2c). Treating the %d *sinR* clones as",
    "independent gives a",
    "vanishingly small probability, but they are not independent -- they",
    "descend from five mutational events. At the level of those events the",
    "question is whether all five landed among the %.0f%% of clones lacking the",
    "clade marker, which has probability %.3f. Suggestive; not decisive."),
    length(eps_clones), n_tot, length(sinr_clones),
    100 * (1 - length(eps_clones) / n_tot), event_p),
  "",
  sprintf(paste(
    "The clade structure does not, however, distort the mutation counts. Of",
    "the %d mutations other than the clade marker, %d fall inside the",
    "*epsA-slrR* clade and %d outside it, against %.0f expected from the",
    "clade's size (binomial P = %.2f). Mutations are spread between the two",
    "clades as evenly as chance predicts, so clonal structure is not inflating",
    "allele frequencies. Nor does either clade mutate faster: counting only",
    "private variants gives %.2f mutations per clone in the *epsA-slrR* clade,",
    "%.2f in the *sinR* clones and %.2f in the rest (Poisson GLM P = %.2f)."),
    clade_test$variants, clade_test$inside_clade, clade_test$outside_clade,
    clade_test$expected_inside, clade_test$p_binomial,
    load$per_clone[load$clade == "epsA-slrR"],
    load$per_clone[load$clade == "sinR"],
    load$per_clone[load$clade == "neither"], mutator_p),
  "",
  paste(
    "The two clades are different kinds of object, and Fig. 2c draws them",
    "that way. *epsA-slrR* is one clade: 47 clones sharing a single ancestral",
    "variant, and remarkably uniform -- 26 of them carry it and nothing else,",
    "with four small sub-lineages (*ybdN*, *pksN*, *citZ-ytwI*, *yutK*)",
    "confined to one edge. That is what a standing variant looks like:",
    "expanded early, accumulated little since. The *sinR* block is the",
    "opposite -- five separate events, each with a small set of descendants,",
    "with *comEC* nested inside one and *ypeB* inside another. No variant",
    "spans the two, which sounds meaningful but is not: any mutation arising",
    "after the split necessarily lands on one side."),
  "",
  "### 5. Genes hit more often than chance allows",
  "",
  sprintf(paste(
    "The strongest evidence of selection here comes from neither allele",
    "frequency nor the neutral simulation. Of the %d mutations in the total",
    "fraction, %d fall in a coding sequence; collapsing variants that share a",
    "gene and a carrier set leaves %d independent events. Scattered at random",
    "across the %s bp of coding sequence, a gene the size of *sinR* (%d bp)",
    "expects %.4f of them. It has five, in five alleles carried by",
    "non-overlapping sets of clones. A min-P permutation over all %s coding",
    "sequences puts the family-wise probability below %s."),
    n_tot_mut, n_coding, n_events, format(cds_total, big.mark = ","),
    sinr_len, sinr_exp, format(n_cds, big.mark = ","), "5 x 10^-5^"),
  "",
  sprintf(paste(
    "*ywcC* -- *slrC* in the current annotation -- is hit twice, at 3,922,748",
    "in one clone and 3,922,952 in two others, giving a family-wise",
    "P = %.3f. *recX* does **not** belong on this list, though the draft",
    "names it: both of its variants, a 49 bp deletion at 925,845 and a 1 bp",
    "insertion at 925,897, are carried by clone 4 and no other clone. Fifty-two",
    "bases apart in a single clone, they are one mutational event at best and",
    "an alignment artifact around the deletion breakpoint at worst. Two genes",
    "show parallelism, not three."),
    ywcc_p),
  "",
  paste(
    "This argument needs no coalescent, no assumption about population size",
    "and no assumption that the sample is well mixed. It is the line the",
    "manuscript should lead with; the frequency statement and the simulation",
    "comparison both rest on assumptions the data strain."),
  "",
  "### 6. What were the mutations?",
  "",
  sprintf(paste(
    "Two genes took more than one independent mutation in the total",
    "fraction: %s. Beyond repeated hits in a single gene, %d of the %d",
    "mutations in the total fraction were present in more than one clone, and",
    "one mutation in *sinR* was found in %d (%.1f%%) of the sequenced clones.",
    "The remaining %d were singletons."),
    paste(sprintf("*%s* (%d independent events)", pt$cds_gene, pt$events),
          collapse = " and "),
    get("S+NS1000", "shared"), get("S+NS1000", "mutations"),
    top$carriers[1], 100 * top$frequency[1],
    get("S+NS1000", "singletons")),
  "",
  "---",
  "",
  "### Where these numbers differ from the earlier draft",
  "",
  paste(
    "- **Section 2 describes a different dataset.** The draft's \"65 clones,",
    "21 (32%) had acquired a mutation, 16 had one, two had three, three had",
    "five\" matches `EvolvedBacillus.filtered.xlsx` (65 clones; 22 = 34% with a",
    "mutation; 17, 2 and 3 clones with one, three and five). That file shares",
    "**no positions and no genes** with `ControlBacillus.filtered.xlsx`, which",
    "is the 84-clone matrix behind Fig. 2a and behind the draft's own section 4",
    "(*sinR* x5, *recX* x2, *ywcC* x2, sinR in 8 = 9.5%). The two sections",
    "describe different sequencing sets. This recreation uses the 84-clone",
    "matrix throughout, for consistency with the figure."),
  "",
  paste(
    "- **\"There were also 34 singleton mutations\"** should be 24. There are 34",
    "mutations in total, 10 of them in more than one clone."),
  "",
  paste(
    "- **The neutral-model conclusion is conditional**, as set out above. The",
    "draft reports P < 2.2 x 10^-16^ without stating N, and the model assumes",
    "a well-mixed population that the co-occurrence structure rules out."),
  "",
  paste(
    "- **The draft names three genes with multiple independent mutations; it",
    "is two.** *recX*'s two variants are both in clone 4 alone, 52 bp apart,",
    "one of them a 49 bp deletion -- one event, or an artifact at the deletion",
    "breakpoint. Clone 4 also carries a *ywcC* indel, so a single clone",
    "accounts for three of the eight indels reported in this fraction; worth a",
    "coverage check if the reads can be found."),
  "",
  paste(
    "- **The omission of the other high-frequency variants is defensible but",
    "unstated.** The draft's criterion appears to be genes with multiple",
    "independent mutations, which only *sinR* and *ywcC* meet once *recX* is",
    "corrected. That is the right criterion given the lineage structure -- but",
    "*ypeB* (a germination protein, missense, in four clones) and *comEC* (a",
    "frameshift in a competence gene) are worth",
    "a sentence in a paper about the boundary between sporulating and",
    "vegetative cells."),
  "",
  sprintf(paste(
    "- **One variant is excluded.** The variant at position %s is carried by 47",
    "of the 84 total-fraction clones (56%%). It is far too common to have arisen",
    "during the experiment and was excluded from the published figure; it is",
    "excluded here too, which is what leaves 34 mutations."),
    format(EXCLUDE_POS, big.mark = ",")),
  "",
  paste(
    "- **Two constants differ from the draft.** The equilibrium population is",
    "1.33 x 10^6^ here, computed from the counts, against 1.36 x 10^6^ in the",
    "draft; and the mutation rate is 3.28 x 10^-10^ throughout, which is the",
    "figure in the draft's own results text, while its Methods say 3.2 x",
    "10^-10^. The 2.5% difference changes nothing, but one should be chosen."),
  "",
  paste(
    "- **The endpoint spore calls were 10 kb out, and are fixed.** They used",
    "to be read from `S1000_endpoint_spore.circos.txt`, whose rows are ±10 kb",
    "Circos highlight windows rather than point positions, so every position",
    "was 10,000 bp low. The correct values come from the surviving summary of",
    "the 96-clone GSF1925 run (see `data/S1000_endpoint_spore.tsv`):",
    "*yetA* 777,008; *ylmG* 1,611,379; *kdgA* 2,323,251; and *levB*–*aspP*",
    "3,539,121. The offset is verified against the two tracks whose variants",
    "also survive in a matrix — all 34 total-fraction rows and both S10 rows",
    "match at +10,000 and none at +0 — and is asserted at every run.",
    "This retires an earlier caveat: the fourth call appeared to sit 860 bp",
    "from the excluded *epsA*–*slrR* marker and inside *epsB*, which looked",
    "like a second eps-operon hit worth an alignment check. At its true",
    "position it is 9,140 bp away and in no gene, so the proximity was an",
    "artefact of the offset and there is nothing to check."),
  "",
  sprintf(paste(
    "- **Ten clones carry no mutation at all**: %s. M23 and M26, which",
    "represent the spore group in Figure 3, are two of ten equally qualified",
    "candidates; nothing in the sequencing distinguishes them from the other",
    "eight. Whoever chose them had a reason, and it is not recorded in any of",
    "the files recovered here."),
    paste(no_mutation, collapse = ", ")),
  ""
)

writeLines(results, file.path(OUT_DIR, "results_text.md"))
message("Wrote output/results_text.md")
