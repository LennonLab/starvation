################################################################################
# 01_mutations.R
#
# Read the variant matrices and assemble one tidy table of mutations.
#
# Input : data/NewControlBacillus.goodcoverage.xlsx   early spore fraction, S10
#         data/ControlBacillus.filtered.xlsx          endpoint total, S+NS1000
#         data/S1000_endpoint_spore.tsv               endpoint spore, S1000
#         data/EvolvedBacillus.filtered.xlsx          a separate 65-clone set
# Output: data/mutations.csv          one row per mutation, with carrier count
#         output/clone_mutation_counts.csv   mutations per clone, per fraction
#
# Each matrix holds one row per variant and two blocks of one column per clone:
# the called allele, then a 0/1 indicator of whether that clone carries the
# variant. The trailing `Total` column is the row sum of the indicator block,
# which is the number of clones carrying that mutation -- the quantity the
# rings are coloured by. The column layouts differ between files, so the
# positions are named here rather than reused.
#
# The indicator block says whether a clone differs from the reference at a
# site, not how. At 3,922,952 in ywcC the two carriers hold different alleles:
# clone 4 has an extra T in an 8-T run and clone 38 one T fewer, which the
# reads confirm (128 of 130 and 105 of 105 spanning reads; tools/site_check.py).
# Those are two mutations, each in one clone, not one mutation in two clones.
# So the allele block is read as well and a site whose carriers disagree is
# split into one row per allele. Positions then no longer identify a mutation
# on their own; `variant_id` does, and the downstream scripts key on it.
#
# No 96-clone matrix survives for the endpoint spore fraction (the GSF1925 run,
# 96 BS_Heatkilled libraries), so its four mutations are read from the summary
# table that does survive -- see the header of data/S1000_endpoint_spore.tsv.
# All four are singletons, so no carrier counts are lost.
#
# They are NOT read from S1000_endpoint_spore.circos.txt. That file is a Circos
# highlight track whose rows are +/-10 kb windows, so its start column sits
# 10,000 bp below the true position. Reading the start as the position put the
# levB-aspP call at 3,529,121, inside epsB and 860 bp from the excluded eps
# marker; the true position is 3,539,121, which is 9,140 bp away and in no
# gene. The offset is verified in check_circos_offset() below.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the mutations project root (or R/).")
source(.setup[1])

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. install.packages(\"readxl\")")
}

#' Read one variant matrix.
#'
#' @param file      the workbook
#' @param pos_col   column holding the genome position
#' @param gene_col  column holding the gene name, or NA
#' @param n_clones  how many clones the matrix covers
#' @param first_geno first column of the 0/1 carrier block
#' @param severity_col,type_col,product_col annotation columns, where present
#' @param ref_col   column holding the reference allele
#'
#' The called-allele block is the n_clones columns just before the carrier
#' block, in the same clone order; that is checked rather than assumed.
read_matrix <- function(file, pos_col, gene_col, n_clones, first_geno,
                        severity_col = NA, type_col = NA, product_col = NA,
                        ref_col) {
  d <- suppressMessages(
    readxl::read_excel(file.path(DATA_DIR, file), sheet = 1,
                       .name_repair = "minimal"))

  pos <- suppressWarnings(as.numeric(d[[pos_col]]))
  keep <- !is.na(pos)                      # trailing notes rows carry no POS

  geno <- as.matrix(d[keep, seq(first_geno, first_geno + n_clones - 1)])
  mode(geno) <- "numeric"
  colnames(geno) <- names(d)[seq(first_geno, first_geno + n_clones - 1)]

  total <- suppressWarnings(as.numeric(d[[ncol(d)]]))[keep]
  if (!isTRUE(all.equal(unname(rowSums(geno)), total))) {
    stop(file, ": the carrier block does not sum to the Total column; ",
         "check the column positions.")
  }

  col_or_na <- function(i) {
    if (is.na(i)) rep(NA_character_, sum(keep)) else as.character(d[[i]][keep])
  }

  allele_cols <- seq(first_geno - n_clones, first_geno - 1)
  if (!identical(names(d)[allele_cols], colnames(geno))) {
    stop(file, ": the allele block does not line up with the carrier block.")
  }
  alleles <- as.matrix(d[keep, allele_cols])

  variants <- data.frame(
    position = as.integer(pos[keep]),
    gene     = col_or_na(gene_col),
    carriers = as.integer(total),
    variant  = col_or_na(type_col),
    severity = col_or_na(severity_col),
    product  = col_or_na(product_col),
    ref      = col_or_na(ref_col),
    stringsAsFactors = FALSE)

  # One row per allele carried. Almost every site has a single allele and
  # passes through unchanged.
  rows <- lapply(seq_len(nrow(variants)), function(k) {
    car <- which(geno[k, ] == 1)
    alt <- alleles[k, car]
    lapply(unique(alt), function(a) {
      g <- geno[k, ] * 0
      g[car[alt == a]] <- 1
      list(v = transform(variants[k, ], alt = a, carriers = sum(alt == a)), g = g)
    })
  })
  rows <- unlist(rows, recursive = FALSE)

  variants <- bind_rows(lapply(rows, `[[`, "v"))
  geno <- do.call(rbind, lapply(rows, `[[`, "g"))
  variants$change <- allele_change(variants$ref, variants$alt)
  split <- duplicated(variants$position) |
    duplicated(variants$position, fromLast = TRUE)
  variants$variant_id <- ifelse(split,
                                paste0(variants$position, variants$change),
                                as.character(variants$position))

  list(variants = variants, genotypes = geno)
}

#' A short description of the change: C>T, +T, -T, -49 bp.
allele_change <- function(ref, alt) {
  out <- rep(NA_character_, length(ref))
  ok <- !is.na(ref) & !is.na(alt)
  r <- ref[ok]; a <- alt[ok]
  nr <- nchar(r); na <- nchar(a)
  short <- function(sign, s) ifelse(nchar(s) <= 3, paste0(sign, s),
                                    sprintf("%s%d bp", sign, nchar(s)))
  out[ok] <- ifelse(nr == 1 & na == 1, paste0(r, ">", a),
             ifelse(na > nr & startsWith(a, r), short("+", substring(a, nr + 1)),
             ifelse(na < nr & startsWith(r, a), short("-", substring(r, na + 1)),
                    "complex")))
  out
}

## ---- the three fractions of the starvation experiment ----------------------

# The two workbooks are laid out differently -- the early-spore file has no
# variant-type column and puts severity and product in unnamed columns -- so
# the positions are given per file rather than reused.
s10 <- read_matrix("NewControlBacillus.goodcoverage.xlsx",
                   pos_col = 5, gene_col = 3, n_clones = 52, first_geno = 59,
                   severity_col = 2, product_col = 4, ref_col = 6)

total <- read_matrix("ControlBacillus.filtered.xlsx",
                     pos_col = 2, gene_col = 7, n_clones = 84, first_geno = 94,
                     severity_col = 6, type_col = 5, product_col = 8,
                     ref_col = 9)

s1000 <- read.delim(file.path(DATA_DIR, "S1000_endpoint_spore.tsv"),
                    comment.char = "#", stringsAsFactors = FALSE)
stopifnot(nrow(s1000) == 4L, !anyNA(s1000$position))

#' Confirm the Circos tracks are +/-10 kb windows, not point positions.
#'
#' The two tracks whose variants also survive in a matrix are checked against
#' it. If the convention ever changes this fails loudly rather than silently
#' shifting four positions by 10 kb.
check_circos_offset <- function(track, truth, label) {
  f <- file.path(DATA_DIR, track)
  if (!file.exists(f)) return(invisible(NULL))
  # the tracks have no trailing newline; that warning is not informative here
  st <- suppressWarnings(read.table(f, sep = "\t"))$V2
  if (!all((st + CIRCOS_HALF_WINDOW) %in% truth)) {
    stop(sprintf("%s: Circos starts do not match %s positions at +%d bp.",
                 track, label, CIRCOS_HALF_WINDOW))
  }
  invisible(NULL)
}

check_circos_offset("LT_NoTreat.mutations.txt", total$variants$position, "total-fraction")
check_circos_offset("ST_Heat.mutations.txt",    s10$variants$position,   "S10")

## ---- assemble --------------------------------------------------------------

muts <- bind_rows(
  data.frame(fraction = "S10", s10$variants, stringsAsFactors = FALSE),
  data.frame(fraction = "S1000",
             position = as.integer(s1000$position), gene = s1000$gene_legacy,
             carriers = 1L, variant = s1000$variant, severity = NA_character_,
             product = s1000$product, stringsAsFactors = FALSE),
  data.frame(fraction = "S+NS1000", total$variants, stringsAsFactors = FALSE)
) %>%
  mutate(fraction = factor(fraction, levels = FRACTION_ORDER),
         variant_id = ifelse(is.na(variant_id), as.character(position),
                             variant_id))

split_sites <- muts %>% filter(variant_id != as.character(position))
if (nrow(split_sites)) {
  message("Sites whose carriers hold different alleles, split into one row per allele:")
  print(as.data.frame(split_sites[, c("fraction", "position", "gene", "ref",
                                      "alt", "carriers")]), row.names = FALSE)
}

excluded <- muts %>% filter(position == EXCLUDE_POS)
muts     <- muts %>% filter(position != EXCLUDE_POS)

muts <- muts %>%
  mutate(n_clones = vapply(as.character(fraction),
                           function(f) FRACTIONS[[f]]$n, integer(1)),
         frequency = carriers / n_clones,
         gene = ifelse(is.na(gene) | !nzchar(gene), NA_character_, gene)) %>%
  arrange(fraction, position)

# Resolve names against the annotation once, so the tables, the figures and
# the parallelism test cannot disagree. `gene` keeps whatever the workbook
# said; `gene_current` is the name to display.
muts <- bind_cols(muts, resolve_genes(muts$position))
muts$gene_display <- ifelse(is.na(muts$gene_current), muts$gene,
                            muts$gene_current)

# Intergenic calls have no CDS to resolve against, so gene_display falls back
# to whatever the workbook wrote -- which is the legacy name, e.g. levB-yveA.
# Rewrite those from the same override table the CDS names come from, so a
# renamed gene reads the same whether it was hit or merely flanks the site.
.ovr <- read.csv(file.path(DATA_DIR, "gene_name_overrides.csv"),
                 comment.char = "#", stringsAsFactors = FALSE)
.ovr <- .ovr[nzchar(.ovr$legacy) & .ovr$legacy != .ovr$preferred, ]
if (nrow(.ovr)) {
  muts$gene_display[is.na(muts$gene_current)] <-
    vapply(muts$gene_display[is.na(muts$gene_current)], function(lbl) {
      if (is.na(lbl)) return(NA_character_)
      parts <- strsplit(lbl, "-", fixed = TRUE)[[1]]
      hit <- match(parts, .ovr$legacy)
      parts[!is.na(hit)] <- .ovr$preferred[hit[!is.na(hit)]]
      paste(parts, collapse = "-")
    }, character(1), USE.NAMES = FALSE)
}

write.csv(muts, file.path(DATA_DIR, "mutations.csv"), row.names = FALSE)

## ---- mutations per clone ---------------------------------------------------

per_clone <- bind_rows(
  data.frame(fraction = "S10", clone = colnames(s10$genotypes),
             mutations = as.integer(colSums(s10$genotypes)),
             stringsAsFactors = FALSE),
  data.frame(fraction = "S+NS1000", clone = colnames(total$genotypes),
             mutations = as.integer(colSums(
               total$genotypes[total$variants$position != EXCLUDE_POS, ,
                               drop = FALSE])),
             stringsAsFactors = FALSE))

write.csv(per_clone, file.path(OUT_DIR, "clone_mutation_counts.csv"),
          row.names = FALSE)

## ---- clone-by-variant matrix ------------------------------------------------
# Long form, and keeping the high-frequency variant, because the lineage
# analysis in R/06_lineage_structure.R needs it: it is the marker that splits
# the total fraction into two backgrounds.

geno_long <- as.data.frame(total$genotypes) %>%
  mutate(position   = total$variants$position,
         gene       = total$variants$gene,
         variant_id = total$variants$variant_id) %>%
  pivot_longer(-c(position, gene, variant_id), names_to = "clone",
               values_to = "carries") %>%
  mutate(clone = sub("^Sample_control_", "", clone)) %>%
  filter(carries == 1) %>%
  select(clone, position, gene, variant_id)

write.csv(geno_long, file.path(DATA_DIR, "genotypes_total_fraction.csv"),
          row.names = FALSE)

# The long form only lists clone-variant pairs that exist, so clones carrying
# nothing drop out of it. Write the full roster separately -- it is the
# denominator for every frequency, and the lineage analysis needs all 84.
write.csv(data.frame(clone = sub("^Sample_control_", "",
                                 colnames(total$genotypes)),
                     stringsAsFactors = FALSE),
          file.path(DATA_DIR, "clones_total_fraction.csv"), row.names = FALSE)

## ---- report ----------------------------------------------------------------

if (nrow(excluded)) {
  message(sprintf(
    "Excluded the variant at %s, carried by %d of %d clones (%.0f%%) -- too ",
    format(excluded$position[1], big.mark = ","), excluded$carriers[1],
    FRACTIONS[["S+NS1000"]]$n,
    100 * excluded$carriers[1] / FRACTIONS[["S+NS1000"]]$n),
    "common to have arisen during the experiment.")
}

## ---- anything suspiciously close to the excluded marker? -------------------
# The excluded variant sits in the eps operon, which is repeat-rich. A second
# call within a couple of kb of it is worth an eye before it is trusted.

NEAR_WINDOW <- 2000L
near <- muts %>%
  filter(abs(position - EXCLUDE_POS) <= NEAR_WINDOW) %>%
  select(fraction, position, gene_current, locus_tag, carriers)

if (nrow(near)) {
  message(sprintf(
    "\nNote: %d retained variant(s) lie within %s bp of the excluded marker at %s:",
    nrow(near), format(NEAR_WINDOW, big.mark = ","),
    format(EXCLUDE_POS, big.mark = ",")))
  print(as.data.frame(near), row.names = FALSE)
  message(paste0(
    "  The marker is intergenic between epsA and slrR; this region is ",
    "repeat-rich.\n  Both calls sit in the eps operon. Probably coincidence, ",
    "but worth checking\n  the alignments before either is relied on."))
  write.csv(near, file.path(OUT_DIR, "variants_near_excluded_marker.csv"),
            row.names = FALSE)
}

cat("\nMutations by fraction\n\n")
summ <- muts %>%
  group_by(fraction) %>%
  summarise(mutations = n(),
            singletons = sum(carriers == 1),
            shared = sum(carriers > 1),
            max_carriers = max(carriers),
            .groups = "drop") %>%
  mutate(clones = vapply(as.character(fraction),
                         function(f) FRACTIONS[[f]]$n, integer(1)),
         max_freq = sprintf("%.1f%%", 100 * max_carriers / clones))
print(as.data.frame(summ), row.names = FALSE)

cat("\nClones carrying at least one mutation\n\n")
pc <- per_clone %>%
  group_by(fraction) %>%
  summarise(clones = n(), with_mutation = sum(mutations > 0), .groups = "drop") %>%
  mutate(percent = sprintf("%.0f%%", 100 * with_mutation / clones))
print(as.data.frame(pc), row.names = FALSE)

message("\nWrote data/mutations.csv and output/clone_mutation_counts.csv")
