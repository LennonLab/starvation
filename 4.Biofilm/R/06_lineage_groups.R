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
#   ywcC/epsA-slrR    m4, m13, m79        (BSU_38220; RefSeq now calls it slrC)
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
  `ywcC/epsA-slrR` = c("m4", "m13", "m79"),
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

#' Pairwise group contrasts with the strain as the unit of replication.
#'
#' This used to be pairwise.t.test() on individual wells, which counts the
#' eight technical replicate wells of a strain as eight independent
#' observations. That put sinR against ywcC/epsA-slrR at P = 1.6e-12; the same
#' contrast with strain as a random effect is P ~ 0.02. The mixed model is the
#' one consistent with R/04_group_models.R, and it is the one reported.
#'
#' log(OD) ~ 0 + group + (1 | strain), REML, Kenward-Roger df, Tukey-adjusted.
strain_contrasts <- function(d) {
  dd <- data.frame(y = log(d$OD_corrected), grp = factor(d$group),
                   strain = factor(d$clone))
  m  <- suppressMessages(lme4::lmer(y ~ 0 + grp + (1 | strain), dd, REML = TRUE))
  em <- suppressMessages(emmeans::emmeans(m, ~ grp, lmer.df = "kenward-roger"))
  pw <- as.data.frame(summary(emmeans::contrast(em, "pairwise", adjust = "tukey"),
                              infer = c(TRUE, TRUE)))
  # the log was taken by hand, so emmeans has nothing to back-transform
  gm <- as.data.frame(summary(em))
  list(model = m,
       pairs = data.frame(contrast = gsub(" - ", " vs ", pw$contrast),
                          ratio = exp(pw$estimate), lower = exp(pw$lower.CL),
                          upper = exp(pw$upper.CL), df = pw$df, p = pw$p.value),
       means = data.frame(group = gm$grp, mean = exp(gm$emmean),
                          lower = exp(gm$lower.CL), upper = exp(gm$upper.CL)))
}

sc <- strain_contrasts(grouped)
pair_p <- function(res, a, b) {
  i <- which(res$pairs$contrast %in% c(paste(a, "vs", b), paste(b, "vs", a)))
  res$pairs$p[i]
}

cat("\nBiofilm by lineage group (blank-corrected OD550, data/biofil.csv)\n\n")
print(summary_tbl %>%
        mutate(mean = sprintf("%.3f", mean), sem = sprintf("%.3f", sem)) %>%
        as.data.frame(), row.names = FALSE)
cat("\nPairwise contrasts, strain as random effect, Kenward-Roger df, Tukey\n\n")
print(transform(sc$pairs, ratio = round(ratio, 2), lower = round(lower, 2),
                upper = round(upper, 2), df = round(df, 1), p = signif(p, 2)),
      row.names = FALSE)

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
sc2 <- strain_contrasts(with_plate2)

sens <- data.frame(
  contrast = c("sinR vs spore", "spore vs ywcC/epsA-slrR", "sinR vs ywcC/epsA-slrR"),
  plate1_only = c(pair_p(sc, "sinR", "spore"),
                  pair_p(sc, "spore", "ywcC/epsA-slrR"),
                  pair_p(sc, "sinR", "ywcC/epsA-slrR")),
  with_plate2 = c(pair_p(sc2, "sinR", "spore"),
                  pair_p(sc2, "spore", "ywcC/epsA-slrR"),
                  pair_p(sc2, "sinR", "ywcC/epsA-slrR")))
write.csv(sens, file.path(OUT_DIR, "lineage_groups_plate2_sensitivity.csv"),
          row.names = FALSE)

cat(sprintf(paste0(
  "\nSensitivity: adding plate 2's six S isolates takes the spore group from\n",
  "%d strains to %d, and sinR vs spore from P = %.2g to P = %.2g (strain as the\n",
  "unit, as above). CAUTION: no strain was read on both plates, so a plate\n",
  "effect cannot be separated from the six S isolates -- any difference between\n",
  "the plates lands entirely on the spore group. Recorded so the dependence is\n",
  "visible, not as a result.\n"),
  summary_tbl$strains[summary_tbl$group == "spore"], length(GROUPS$spore),
  sens$plate1_only[1], sens$with_plate2[1]))

contrasts <- transmute(sc$pairs, Comparison = contrast, p = p)

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
           "never carry one in *ywcC* (BSU_38220) or upstream of *slrR*, so the ",
           "likely-vegetative clones divide in two. Means \u00b1 SEM are over ",
           "wells and are descriptive. Tests treat the strain as the unit of ",
           "replication -- a linear mixed model of log OD with strain as a random ",
           "effect, Kenward\u2013Roger degrees of freedom, Tukey-adjusted: ",
           paste(sprintf("%s P = %s", contrasts$Comparison,
                         format.pval(contrasts$p, digits = 2)),
                 collapse = "; "),
           ". The *sinR*-versus-spore comparison is the weak one: the spore ",
           "group is two strains, and they differ 16-fold from each other."))

## ---- figure ---## ---- figure ----------------------------------------------------------------

grouped$group <- factor(grouped$group, levels = summary_tbl$group)

# Wells in light grey are the raw data. The larger grey points are strain means
# -- the unit the test uses. Black is the mixed model's group mean and 95% CI,
# which replaces an error bar that used to be the SE of pooled wells.
strain_means <- grouped %>%
  group_by(group, clone) %>%
  summarise(OD = exp(mean(log(OD_corrected))), .groups = "drop")
grp_est <- sc$means
grp_est$group <- factor(grp_est$group, levels = summary_tbl$group)
strain_means$group <- factor(strain_means$group, levels = summary_tbl$group)

fig <- ggplot(grouped, aes(group, OD_corrected)) +
  geom_jitter(width = 0.14, height = 0, colour = "grey85", size = 1.6, alpha = 0.8) +
  geom_point(data = strain_means, aes(group, OD), colour = "grey45", size = 2.6,
             position = position_nudge(x = 0.22)) +
  geom_errorbar(data = grp_est, aes(x = group, ymin = lower, ymax = upper),
                inherit.aes = FALSE, width = 0.14, linewidth = 0.7) +
  geom_point(data = grp_est, aes(group, mean), inherit.aes = FALSE, size = 3.6,
             shape = 21, colour = "black", fill = "white", stroke = 0.9) +
  scale_y_log10(sec.axis = dup_axis()) +
  labs(x = NULL, y = expression("Biofilm (OD"[550] * ", blank corrected)")) +
  mytheme

ggsave(file.path(OUT_DIR, "fig_biofilm_lineage_groups.pdf"), fig,
       width = 5.6, height = 4.6)
ggsave(file.path(OUT_DIR, "fig_biofilm_lineage_groups.png"), fig,
       width = 5.6, height = 4.6, dpi = 200)
message("Wrote output/fig_biofilm_lineage_groups.pdf and .png")
