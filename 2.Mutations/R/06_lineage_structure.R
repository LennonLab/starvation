################################################################################
# 06_lineage_structure.R
#
# The total fraction is not one well-mixed population.
#
# Input : data/genotypes_total_fraction.csv, data/mutations.csv
# Output: output/clone_backgrounds.csv     each clone's genetic background
#         output/cooccurrence_tests.csv    mutual exclusivity and nesting
#         output/fig2c_lineages.pdf / .png
#
# The variant at epsA-slrR is carried by 47 of the 84 sequenced clones. It is
# too common to have arisen during the experiment, so the figures exclude it --
# but excluding it hides what it marks. No clone carrying a sinR mutation
# carries it, and vice versa: the fraction splits into two backgrounds that
# share no derived variants.
#
# That matters for how the other variants are read. A mutation in four clones
# can be one event carried up by an expanding lineage rather than four
# independent hits, and the neutral model in R/03_null_model.R assumes a single
# well-mixed population, which this violates. Clonal structure raises allele
# frequencies with no selection at all.
#
# sinR is the exception: distinct alleles at distinct positions in
# non-overlapping sets of clones is parallelism, and clonal structure does not
# produce it.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the mutations project root (or R/).")
source(.setup[1])

geno <- read.csv(file.path(DATA_DIR, "genotypes_total_fraction.csv"),
                 colClasses = c("character", "integer", "character"))
muts <- read.csv(file.path(DATA_DIR, "mutations.csv"), stringsAsFactors = FALSE)

# Every sequenced clone, including the ones carrying no variant at all -- the
# long genotype table only has rows where a clone carries something.
ALL_CLONES <- read.csv(file.path(DATA_DIR, "clones_total_fraction.csv"),
                       colClasses = "character")$clone
N_CLONES   <- FRACTIONS$`S+NS1000`$n
stopifnot(length(ALL_CLONES) == N_CLONES,
          all(geno$clone %in% ALL_CLONES))

carriers_of <- function(position) geno$clone[geno$position == position]
carriers_of_gene <- function(g) unique(geno$clone[geno$gene %in% g])

## ---- the two backgrounds ---------------------------------------------------

eps_clones  <- carriers_of(EXCLUDE_POS)
sinr_clones <- carriers_of_gene("sinR")

# The two blocks are not the same kind of object, and the labels say so. The
# epsA-slrR block is one clade: 47 clones sharing a single ancestral variant.
# The sinR block is five independent mutational events that share only the
# absence of that variant, so each clone is named by the allele it carries.
# What is left over splits again, into clones with private variants only and
# clones carrying nothing at all.
sinr_allele_of <- function(clone) {
  p <- geno$position[geno$clone == clone & geno$gene == "sinR"]
  if (!length(p)) NA_character_ else sprintf("sinR %s", format(p[1], big.mark = ","))
}

backgrounds <- data.frame(
  clone = ALL_CLONES,
  background = vapply(ALL_CLONES, function(cl) {
    if (cl %in% eps_clones)  return("epsA-slrR clade")
    if (cl %in% sinr_clones) return(sinr_allele_of(cl))
    if (cl %in% geno$clone)  return("no clade marker")
    "no mutations"
  }, character(1)),
  stringsAsFactors = FALSE)

backgrounds$clade <- ifelse(backgrounds$background == "epsA-slrR clade",
                            "epsA-slrR",
                     ifelse(grepl("^sinR", backgrounds$background), "sinR",
                            backgrounds$background))

# A three-way version for the mutation-load comparison. The finer split above
# separates clones by whether they carry a private variant, which is the very
# thing being counted -- comparing on that would be circular.
backgrounds$clade3 <- ifelse(backgrounds$clade %in% c("epsA-slrR", "sinR"),
                             backgrounds$clade, "neither")

write.csv(backgrounds, file.path(OUT_DIR, "clone_backgrounds.csv"),
          row.names = FALSE)

tab <- table(sinR = ALL_CLONES %in% sinr_clones,
             epsA = ALL_CLONES %in% eps_clones)
excl <- fisher.test(tab)

# The clone-level test above overstates the case: the 20 sinR clones are not
# 20 independent observations, they descend from five mutational events. At
# the event level the question is whether five mutations all landed in the
# clones that lack epsA-slrR, which is (1 - p)^5 for a clade covering p of
# the sample -- suggestive, not decisive.
n_sinr_events <- length(unique(geno$position[geno$gene == "sinR"]))
p_clade  <- length(eps_clones) / N_CLONES
event_p  <- (1 - p_clade)^n_sinr_events

cat("The two backgrounds do not co-occur\n\n")
cat(sprintf("  epsA-slrR carriers : %d of %d\n", length(eps_clones), N_CLONES))
cat(sprintf("  sinR carriers      : %d of %d\n", length(sinr_clones), N_CLONES))
cat(sprintf("  carrying both      : %d\n", length(intersect(eps_clones, sinr_clones))))
cat(sprintf("  P treating clones as independent : %s -- overstates it\n",
            format.pval(excl$p.value, digits = 3)))
cat(sprintf("  P at the level of the %d sinR mutational events : %.3f = %.2f^%d\n\n",
            n_sinr_events, event_p, 1 - p_clade, n_sinr_events))

## ---- where does each mutation sit? -----------------------------------------

shared <- muts %>%
  filter(fraction == "S+NS1000", carriers > 1) %>%
  arrange(desc(carriers))

# For each shared mutation: which background its carriers belong to, and
# whether its carriers are contained within some other mutation's carriers.
nesting <- lapply(seq_len(nrow(shared)), function(i) {
  p  <- shared$position[i]
  cl <- carriers_of(p)

  bg <- backgrounds$background[match(cl, backgrounds$clone)]
  bg_label <- if (all(bg == bg[1])) bg[1] else "mixed"

  # Any other variant whose carriers contain all of this one's.
  others <- setdiff(unique(geno$position), p)
  within <- vapply(others, function(q) {
    qc <- carriers_of(q)
    length(qc) > length(cl) && all(cl %in% qc)
  }, logical(1))

  parents <- others[within]
  parent_lab <- if (!length(parents)) NA_character_ else
    paste(sprintf("%s@%s", geno$gene[match(parents, geno$position)],
                  format(parents, big.mark = ",")), collapse = "; ")

  data.frame(
    gene = ifelse(is.na(shared$gene[i]), "-", shared$gene[i]),
    position = p, carriers = shared$carriers[i],
    clones = paste(cl, collapse = ","),
    background = bg_label,
    nested_within = parent_lab,
    stringsAsFactors = FALSE)
})
nesting <- bind_rows(nesting)

write.csv(nesting, file.path(OUT_DIR, "cooccurrence_tests.csv"), row.names = FALSE)

cat("Shared mutations, by background and containment\n\n")
print(nesting %>%
        mutate(nested_within = ifelse(is.na(nested_within), "-", nested_within),
               position = format(position, big.mark = ",")),
      row.names = FALSE)

## ---- is one clade simply mutating faster? ----------------------------------
# The obvious alternative to selection is that one background carries more
# mutations because it mutates more. Counting only private (singleton)
# variants, which carry no signal of expansion, settles it.

singleton_pos <- geno %>%
  filter(position != EXCLUDE_POS) %>%
  count(position) %>%
  filter(n == 1) %>%
  pull(position)

load_df <- backgrounds %>%
  mutate(singletons = vapply(clone, function(cl)
    sum(geno$position[geno$clone == cl] %in% singleton_pos), integer(1)))

load_summary <- load_df %>%
  group_by(clade = clade3) %>%
  summarise(clones = n(), per_clone = mean(singletons),
            total = sum(singletons), .groups = "drop")

mutator_p <- anova(glm(singletons ~ clade3, poisson, load_df),
                   test = "Chisq")$`Pr(>Chi)`[2]

write.csv(load_summary, file.path(OUT_DIR, "mutation_load_by_clade.csv"),
          row.names = FALSE)
writeLines(format(mutator_p, digits = 4),
           file.path(OUT_DIR, "mutator_test_p.txt"))

cat("\nPrivate (singleton) mutations per clone, by background\n\n")
print(load_summary %>% mutate(per_clone = round(per_clone, 3)) %>% as.data.frame(),
      row.names = FALSE)
cat(sprintf("\nPoisson GLM, background effect: P = %.2f -- indistinguishable,\n",
            mutator_p))
cat("so the difference between clades is not a difference in mutation rate.\n")

## ---- the sinR alleles ------------------------------------------------------

sinr_alleles <- geno %>%
  filter(gene == "sinR") %>%
  group_by(position) %>%
  summarise(clones = paste(sort(as.integer(clone)), collapse = ","),
            n = n(), .groups = "drop") %>%
  arrange(position)

overlaps <- combn(nrow(sinr_alleles), 2, function(ij) {
  a <- carriers_of(sinr_alleles$position[ij[1]])
  b <- carriers_of(sinr_alleles$position[ij[2]])
  length(intersect(a, b))
})

cat(sprintf(
  "\nsinR carries %d alleles at %d positions across %d clones; the largest\n",
  sum(sinr_alleles$n), nrow(sinr_alleles), length(sinr_clones)))
cat(sprintf("overlap between any two allele sets is %d clone(s).\n",
            max(overlaps)))
print(as.data.frame(sinr_alleles), row.names = FALSE)

## ---- the clones carrying nothing -------------------------------------------
# Worth recording because Figure 3 picks two clones to represent the spore
# group, and on this evidence they are two of ten equally qualified ones.

no_mutation <- sort(as.integer(setdiff(ALL_CLONES, geno$clone)))
writeLines(as.character(no_mutation),
           file.path(OUT_DIR, "clones_without_mutations.txt"))

cat(sprintf("\n%d clones carry no variant at all: %s\n",
            length(no_mutation), paste(no_mutation, collapse = ", ")))
cat("Nothing in the sequencing distinguishes any one of them from the others.\n")

## ---- figure ----------------------------------------------------------------
# Clone-by-mutation presence, clones grouped by background. This is the
# lineage structure the frequency rings cannot show.

# Show every variant carried by more than one clone, plus every variant in a
# gene that was hit more than once -- which brings in the two singleton sinR
# alleles, and with them the full picture of five independent hits.
repeat_genes <- muts %>%
  filter(fraction == "S+NS1000", !is.na(gene)) %>%
  count(gene) %>% filter(n > 1) %>% pull(gene)

show_pos <- sort(unique(c(
  EXCLUDE_POS,
  shared$position,
  geno$position[geno$gene %in% repeat_genes])))

locus_label <- function(p) sprintf("%s (%s)",
                                   geno$gene[match(p, geno$position)],
                                   format(p, big.mark = ","))

present <- geno %>% filter(position %in% show_pos)

# Order clones by genotype within each block, so clones sharing a variant sit
# together and the nesting is visible; the ordering key is the presence
# pattern read as a string.
pattern <- vapply(ALL_CLONES, function(cl) {
  paste(as.integer(show_pos %in% present$position[present$clone == cl]),
        collapse = "")
}, character(1))

BLOCK_ORDER <- c("epsA-slrR clade",
                 sort(unique(backgrounds$background[
                   grepl("^sinR", backgrounds$background)])),
                 "no clade marker", "no mutations")

clone_order <- backgrounds %>%
  mutate(pattern = pattern[match(clone, ALL_CLONES)],
         block = factor(background, levels = BLOCK_ORDER)) %>%
  arrange(block, desc(pattern), as.integer(clone)) %>%
  pull(clone)

# One colour family per kind of block: blue for the single clade, a magma ramp
# across the five independent sinR alleles, greys for the remainder.
sinr_blocks <- BLOCK_ORDER[grepl("^sinR", BLOCK_ORDER)]
BG_COLOURS <- c(
  setNames("#3B6EA5", "epsA-slrR clade"),
  setNames(viridisLite::magma(length(sinr_blocks), begin = 0.35, end = 0.8),
           sinr_blocks),
  setNames(c("grey62", "grey86"), c("no clade marker", "no mutations")))

BG_ROW  <- "Background"
SNG_ROW <- "Private mutations"
row_levels <- rev(c(BG_ROW, SNG_ROW, locus_label(show_pos)))

tiles <- bind_rows(
  backgrounds %>% transmute(clone, row = BG_ROW, fill = background),
  present %>% transmute(clone, row = locus_label(position), fill = "mutation")
) %>%
  mutate(clone = factor(clone, levels = clone_order),
         row   = factor(row, levels = row_levels),
         fill  = factor(fill, levels = c(names(BG_COLOURS), "mutation")))

counts <- load_df %>%
  filter(singletons > 0) %>%
  transmute(clone = factor(clone, levels = clone_order),
            row = factor(SNG_ROW, levels = row_levels),
            label = as.character(singletons))

fig <- ggplot(tiles, aes(clone, row)) +
  geom_tile(aes(fill = fill), width = 0.85, height = 0.78) +
  geom_text(data = counts, aes(label = label), size = 2.4, colour = "grey20") +
  scale_fill_manual(values = c(BG_COLOURS, mutation = "grey15"),
                    breaks = BLOCK_ORDER, name = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = sprintf("Clone (n = %d), grouped by background then genotype", N_CLONES),
       y = NULL,
       caption = paste0(
         "The two blocks are different kinds of thing: epsA-slrR is one clade ",
         "of 47 clones sharing an ancestral variant, while the sinR block is ",
         "five\nindependent mutations that share only the absence of it, ",
         "shown here one shade per allele. Private mutations per clone are ",
         sprintf("%.2f, %.2f and %.2f\nin the epsA-slrR, sinR and remaining ",
                 load_summary$per_clone[load_summary$clade == "epsA-slrR"],
                 load_summary$per_clone[load_summary$clade == "sinR"],
                 mean(load_df$singletons[!load_df$clade %in% c("epsA-slrR", "sinR")])),
         sprintf("clones (Poisson GLM P = %.2f), so neither clade mutates faster.",
                 mutator_p))) +
  mytheme +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 8),
        legend.position = "top",
        legend.text = element_text(size = 8),
        plot.caption = element_text(size = 7, colour = "grey35", hjust = 0))

ggsave(file.path(OUT_DIR, "fig2c_lineages.pdf"), fig, width = 9, height = 5.6)
ggsave(file.path(OUT_DIR, "fig2c_lineages.png"), fig, width = 9, height = 5.6,
       dpi = 200)

message("\nWrote output/clone_backgrounds.csv, cooccurrence_tests.csv and fig2c_lineages")
