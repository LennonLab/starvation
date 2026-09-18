################################################################################
# 06_lineage_groups.R
#
# Biofilm formation split by which biofilm regulator the clone carries.
#
# Input : data/biofil.csv              the analysis dataset (plate 1)
#         data/Biofilm_06_16_20.txt    the raw two-plate reading
# Output: output/table4_lineage_groups.md
#         output/fig_biofilm_lineage_groups.pdf / .png
#
# An early manuscript draft splits the likely-vegetative clones by regulator --
# "clones containing a mutation in sinR did not also contain mutations in ywcC
# or upstream of slrR" -- and compares:
#
#   sinR              m17, m19, m21, m41, m54
#   slrC/epsA-slrR    m4, m13, m79        (the draft writes ywcC/epsA-slrR)
#   spore             m23, m26
#
# "Upstream of slrR" is the intergenic epsA-slrR variant at 3,529,981 -- the
# one carried by 47 of 84 sequenced clones and excluded from Figure 2. So the
# "/slrR" in the biofil.csv label was never the slrR gene; it names that
# marker, and the draft's two subgroups are the two clades the sequencing
# finds. The mutual exclusivity was already known.
#
# The draft's own numbers are from an earlier round of the assay and are not
# reproduced here.
#
# On scope: data/biofil.csv is plate 1 of a two-plate reading. Plate 2 holds
# S1, S6, S11, S22, S51 and S95 in alternating columns with paired media
# controls -- a different layout, and not part of this analysis. It is read
# below only as a sensitivity check, because it is the one thing that changes
# the sinR-versus-spore result.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the biofilm project root (or R/).")
source(.setup[1])

RAW <- file.path(DATA_DIR, "Biofilm_06_16_20.txt")

## ---- primary: the analysis dataset -----------------------------------------

biofilm <- read.csv(file.path(DATA_DIR, "biofil.csv"), stringsAsFactors = FALSE) %>%
  transmute(plate = Plate, well = Replicate, clone = clones,
            OD_corrected = OD550_Corrected)

GROUPS <- list(
  sinR             = c("m17", "m19", "m21", "m41", "m54"),
  `slrC/epsA-slrR` = c("m4", "m13", "m79"),
  spore            = c("m23", "m26"))

assign_groups <- function(d) {
  d$group <- NA_character_
  for (g in names(GROUPS)) d$group[d$clone %in% GROUPS[[g]]] <- g
  d[!is.na(d$group), ]
}

grouped <- assign_groups(biofilm)

summary_tbl <- grouped %>%
  group_by(group) %>%
  summarise(strains = n_distinct(clone), wells = n(),
            mean = mean(OD_corrected), sem = sem(OD_corrected), .groups = "drop") %>%
  arrange(mean)

pw <- pairwise.t.test(grouped$OD_corrected, grouped$group,
                      p.adjust.method = "holm")

cat("\nBiofilm by lineage group (blank-corrected OD550, data/biofil.csv)\n\n")
print(summary_tbl %>%
        mutate(mean = sprintf("%.3f", mean), sem = sprintf("%.3f", sem)) %>%
        as.data.frame(), row.names = FALSE)
cat("\nPairwise t-tests, pooled SD, Holm-adjusted\n\n")
print(signif(pw$p.value, 3))

## ---- sensitivity: what plate 2 would do ------------------------------------
# The raw file stacks two plates, each with its own header row. Plate 2's six
# spore isolates are not part of the analysis above; this is only to record
# what including them would change.

lines <- readLines(RAW)
header_rows <- grep("\tReplicate\t", lines)

read_plate <- function(hdr) {
  cols <- strsplit(lines[hdr], "\t")[[1]]
  body <- do.call(rbind, strsplit(lines[(hdr + 1):(hdr + 8)], "\t"))
  colnames(body) <- cols
  d <- as.data.frame(body, stringsAsFactors = FALSE)
  names(d)[1:2] <- c("plate", "well")
  pivot_longer(d, -c(plate, well), names_to = "strain", values_to = "OD") %>%
    mutate(OD = as.numeric(OD))
}

plates <- bind_rows(lapply(header_rows, read_plate))

# Blank correction is the median, not the mean: plate 1 has two contaminated
# blank wells (0.262 and 0.179 against a baseline near 0.064). The median
# gives 0.0645, which is what data/biofil.csv carries.
blanks <- plates %>%
  filter(grepl("^Blank", strain)) %>%
  group_by(plate) %>%
  summarise(blank = median(OD), .groups = "drop")

all_strains <- plates %>%
  filter(!grepl("^Blank", strain)) %>%
  left_join(blanks, by = "plate") %>%
  mutate(clone = tolower(strain), OD_corrected = OD - blank)

check <- all_strains %>%
  filter(plate == "plate1") %>%
  inner_join(biofilm, by = c("clone", "well"), suffix = c("", "_ref"))

message(sprintf("Plate 1 of the raw file reproduces biofil.csv: %s (%d wells, max diff %.3g)",
                isTRUE(all.equal(check$OD_corrected, check$OD_corrected_ref,
                                 tolerance = 1e-8)),
                nrow(check), max(abs(check$OD_corrected - check$OD_corrected_ref))))

GROUPS$spore <- c(GROUPS$spore, "s1", "s6", "s11", "s22", "s51", "s95")
with_plate2 <- assign_groups(all_strains)
pw2 <- pairwise.t.test(with_plate2$OD_corrected, with_plate2$group,
                       p.adjust.method = "holm")

cat(sprintf(paste0(
  "\nSensitivity: adding plate 2's six S isolates to the spore group takes it\n",
  "from %d strains to %d, and sinR vs spore from P = %.2g to P = %.2g. The\n",
  "other two comparisons stay significant either way. Plate 2 is a different\n",
  "layout and is not part of the analysis; this is recorded only so the\n",
  "dependence is visible.\n"),
  summary_tbl$strains[summary_tbl$group == "spore"],
  length(GROUPS$spore),
  pw$p.value["spore", "sinR"], pw2$p.value["spore", "sinR"]))

contrasts <- as.data.frame(as.table(pw$p.value)) %>%
  filter(!is.na(Freq)) %>%
  transmute(Comparison = paste(Var1, "vs", Var2), p = Freq)

write_md <- function(df, file, title, caption) {
  sep <- ifelse(vapply(df, is.numeric, logical(1)), "---:", ":---")
  writeLines(c(paste("##", title), "", caption, "",
               paste0("| ", paste(names(df), collapse = " | "), " |"),
               paste0("| ", paste(sep, collapse = " | "), " |"),
               apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")),
               ""),
             file.path(OUT_DIR, file))
  message("Wrote output/", file)
}

write_md(as.data.frame(summary_tbl %>%
           transmute(Group = group, Strains = strains, Wells = wells,
                     `Mean OD550` = sprintf("%.3f ± %.3f", mean, sem))),
         "table4_lineage_groups.md",
         "Table 4. Biofilm formation by lineage group",
         paste0(
           "Blank-corrected OD550 from `data/biofil.csv`. Groups follow the ",
           "split an early draft describes: clones carrying a *sinR* mutation ",
           "never carry one in *slrC* (*ywcC*) or upstream of *slrR*, so the ",
           "likely-vegetative clones divide in two. Holm-adjusted pairwise ",
           "t-tests on pooled SD: ",
           paste(sprintf("%s P = %s", contrasts$Comparison,
                         format.pval(contrasts$p, digits = 2)),
                 collapse = "; "),
           ". The *sinR*-versus-spore comparison is the weak one, and it is ",
           "weak because the two spore clones differ 16-fold from each other."))

## ---- figure ---## ---- figure ----------------------------------------------------------------

grouped$group <- factor(grouped$group, levels = summary_tbl$group)

fig <- ggplot(grouped, aes(group, OD_corrected)) +
  geom_jitter(width = 0.14, height = 0, colour = "grey72", size = 2, alpha = 0.75) +
  stat_summary(fun = mean, geom = "point", size = 3.6, shape = 21,
               colour = "black", fill = "white", stroke = 0.9) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.16,
               linewidth = 0.6) +
  scale_y_log10(sec.axis = dup_axis()) +
  labs(x = NULL, y = expression("Biofilm (OD"[550] * ", blank corrected)")) +
  mytheme

ggsave(file.path(OUT_DIR, "fig_biofilm_lineage_groups.pdf"), fig,
       width = 5.6, height = 4.6)
ggsave(file.path(OUT_DIR, "fig_biofilm_lineage_groups.png"), fig,
       width = 5.6, height = 4.6, dpi = 200)
message("Wrote output/fig_biofilm_lineage_groups.pdf and .png")
