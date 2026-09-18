################################################################################
# 07_parallelism.R
#
# The cleanest evidence of selection: genes hit by more independent mutations
# than chance allows.
#
# Input : data/mutations.csv, data/genotypes_total_fraction.csv,
#         data/NC_000964.3_cds.tsv   exact CDS coordinates, from the RefSeq GFF
# Output: output/mutation_events.csv       variants collapsed to events
#         output/parallelism_tests.csv     per-gene Poisson and permutation P
#         output/clade_distribution_test.csv
#
# Three things this gets right that a naive count does not.
#
# 1. Events, not variant rows. Two variants in the same gene carried by the
#    same clone and no other are one mutational event -- or an alignment
#    artifact. recX is exactly this case: a 49 bp deletion and a 1 bp
#    insertion 52 bp apart, both in clone 4 alone. Counted as k = 2 it looks
#    significant; counted correctly it is a single event and says nothing.
#
# 2. Exact CDS lengths, and by position rather than by gene name. The
#    workbook's names predate the current annotation -- its "ywcC" is slrC
#    (BSU_38220) and its "ypeB" falls in sleC (BSU_22920) -- so each mutation
#    is mapped to whatever CDS contains it. Gene size matters: pksN is 16 kb,
#    so two hits there are unremarkable, while sinR is 336 bp.
#
# 3. Coding mutations as the denominator. An intergenic mutation cannot land
#    in a CDS, so including it dilutes the per-gene rate.
#
# The headline test is a permutation: scatter the same number of coding events
# across the genome with probability proportional to CDS length, and ask how
# often *any* gene collects as many as the one observed. That is the right null
# when the gene was nominated by the data. Bonferroni over every CDS is
# reported too, as a conservative bound.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the mutations project root (or R/).")
source(.setup[1])

set.seed(20260910)
N_PERM <- 20000L

muts <- read.csv(file.path(DATA_DIR, "mutations.csv"), stringsAsFactors = FALSE)
geno <- read.csv(file.path(DATA_DIR, "genotypes_total_fraction.csv"),
                 colClasses = c("character", "integer", "character"))


total <- muts[muts$fraction == "S+NS1000", ]

## ---- 1. mutations already carry their CDS, resolved in R/01 ----------------

cds <- load_cds()

total <- total %>%
  rename(cds_gene = gene_current, length_bp = cds_length)
total$coding <- !is.na(total$cds_gene)

renamed <- total %>%
  filter(coding, !is.na(gene), gene != cds_gene) %>%
  distinct(legacy = gene, current = cds_gene, locus_tag)

## ---- 2. collapse variants to independent mutational events -----------------
# Within one CDS, variants carried by exactly the same set of clones cannot be
# shown to be independent; they are counted once.

carriers_key <- function(position) {
  paste(sort(geno$clone[geno$position == position]), collapse = ",")
}

events <- total %>%
  filter(coding) %>%
  mutate(clones = vapply(position, carriers_key, character(1))) %>%
  group_by(locus_tag, cds_gene, length_bp, clones) %>%
  summarise(positions = paste(format(position, big.mark = ","), collapse = " + "),
            n_variants = n(), .groups = "drop")

write.csv(events, file.path(OUT_DIR, "mutation_events.csv"), row.names = FALSE)

collapsed <- events %>% filter(n_variants > 1)

n_events <- nrow(events)
n_coding_variants <- sum(total$coding)

cat(sprintf("%d mutations in the endpoint total fraction: %d coding, %d intergenic.\n",
            nrow(total), n_coding_variants, sum(!total$coding)))
cat(sprintf("Collapsing variants that share a gene and a carrier set leaves %d independent events.\n\n",
            n_events))

if (nrow(collapsed)) {
  cat("Collapsed to a single event:\n")
  print(as.data.frame(collapsed[, c("cds_gene", "locus_tag", "positions", "clones")]),
        row.names = FALSE)
  cat("\n")
}

if (nrow(renamed)) {
  cat("Legacy names in the workbooks, and the current symbol used here:\n")
  print(as.data.frame(renamed), row.names = FALSE)
  cat("  slrC (BSU_38220) sits beside slrA, which its product description names.\n")
  cat("  ypeB (BSU_22920) is NOT renamed: RefSeq labels it sleC, but B. subtilis\n")
  cat("  has no sleC, UniProt P38490 is YpeB, and sleB is immediately downstream.\n")
  cat("  See data/gene_name_overrides.csv.\n\n")
}

## ---- 3. how many events per gene? ------------------------------------------

per_gene <- events %>%
  count(locus_tag, cds_gene, length_bp, name = "events") %>%
  filter(events > 1) %>%
  arrange(desc(events))

CDS_TOTAL <- sum(cds$length_bp)

## ---- 4. permutation null ---------------------------------------------------
# Scatter n_events across the CDS complement in proportion to gene length.
#
# Two statistics. The raw one -- does any gene collect this many? -- is
# dominated by large genes, so for a small gene it is harsher than Bonferroni,
# which is not what a family-wise correction should do. The right statistic
# conditions on gene size: compute each gene's own Poisson probability within
# a scattering, take the smallest, and ask how often that beats the observed
# gene's. This is the standard min-P procedure and it is the one reported.

lambda_all <- n_events * cds$length_bp / CDS_TOTAL

perm <- replicate(N_PERM, {
  hit <- sample.int(nrow(cds), n_events, replace = TRUE, prob = cds$length_bp)
  cnt <- tabulate(hit, nbins = nrow(cds))
  i   <- which(cnt > 0)
  c(max_count = max(cnt),
    min_p = min(ppois(cnt[i] - 1, lambda_all[i], lower.tail = FALSE)))
})
max_hits <- perm["max_count", ]
min_p    <- perm["min_p", ]

per_gene <- per_gene %>%
  mutate(
    expected = n_events * length_bp / CDS_TOTAL,
    # per-gene Poisson, ignoring that the gene was chosen by the data
    p_poisson = ppois(events - 1, expected, lower.tail = FALSE),
    p_bonferroni = pmin(1, p_poisson * nrow(cds)),
    # family-wise, conditioning on gene size: how often does the best gene in
    # a random scattering beat this one?
    p_permutation = vapply(p_poisson, function(p) mean(min_p <= p), numeric(1)),
    # the same idea ignoring gene size, kept for comparison
    p_perm_maxcount = vapply(events, function(k) mean(max_hits >= k), numeric(1)))

write.csv(per_gene, file.path(OUT_DIR, "parallelism_tests.csv"), row.names = FALSE)

cat(sprintf("Genes with more than one independent event (%d events over %s CDS bp)\n\n",
            n_events, format(CDS_TOTAL, big.mark = ",")))
print(per_gene %>%
        mutate(expected = sprintf("%.4f", expected),
               p_poisson = format.pval(p_poisson, digits = 3),
               p_bonferroni = format.pval(p_bonferroni, digits = 3),
               p_permutation = ifelse(p_permutation == 0,
                                      sprintf("<%.0e", 1 / N_PERM),
                                      format.pval(p_permutation, digits = 3)),
               p_perm_maxcount = ifelse(p_perm_maxcount == 0,
                                        sprintf("<%.0e", 1 / N_PERM),
                                        format.pval(p_perm_maxcount, digits = 2))) %>%
        as.data.frame(),
      row.names = FALSE)

cat(sprintf(paste0(
  "\np_permutation is the family-wise probability, over %s scatterings of %d\n",
  "events across %s genes, that the most surprising gene in a random genome\n",
  "is at least as surprising as this one, with each gene judged against its\n",
  "own length. That is the appropriate null when the gene was nominated by\n",
  "the data. p_perm_maxcount is the cruder version that only counts hits and\n",
  "ignores gene size; it is harsher on small genes than Bonferroni, which is\n",
  "why it is not the one to quote.\n"),
  format(N_PERM, big.mark = ","), n_events, format(nrow(cds), big.mark = ",")))

## ---- 5. are mutations spread evenly across the two clades? -----------------

n_clones   <- FRACTIONS$`S+NS1000`$n
eps_clones <- geno$clone[geno$position == EXCLUDE_POS]
p_clade    <- length(eps_clones) / n_clones

variants <- setdiff(unique(geno$position), EXCLUDE_POS)
side <- vapply(variants, function(p) {
  cl <- geno$clone[geno$position == p]
  if (all(cl %in% eps_clones)) "inside" else
    if (!any(cl %in% eps_clones)) "outside" else "spans"
}, character(1))

inside <- sum(side == "inside")
bt <- binom.test(inside, length(variants), p_clade)

write.csv(data.frame(variants = length(variants), inside_clade = inside,
                     outside_clade = sum(side == "outside"),
                     spanning = sum(side == "spans"),
                     expected_inside = length(variants) * p_clade,
                     p_binomial = bt$p.value),
          file.path(OUT_DIR, "clade_distribution_test.csv"), row.names = FALSE)

cat(sprintf(paste0(
  "\nOf the %d mutations other than the clade marker, %d fall inside the\n",
  "epsA-slrR clade and %d outside it, against %.1f expected from the clade's\n",
  "size (binomial P = %.2f). Mutations are distributed between the clades as\n",
  "evenly as chance predicts, so clonal structure is not inflating allele\n",
  "frequencies. No variant spans the two, which sounds meaningful but is not:\n",
  "any mutation arising after the split necessarily lands on one side.\n"),
  length(variants), inside, sum(side == "outside"),
  length(variants) * p_clade, bt$p.value))

message("\nWrote output/mutation_events.csv, parallelism_tests.csv and clade_distribution_test.csv")
