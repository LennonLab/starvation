################################################################################
# 03_supplementary.R
#
# Supplementary figures, collected into output/ so that everything going to the
# manuscript sits in one folder.
#
# Most of these are single panels that their own project already draws well.
# Nothing is redrawn: each is pulled out of the project script that owns it,
# through load_project(), and re-saved under a supplementary name and a
# consistent size. Add an entry to SUPPLEMENTARY below to add a figure.
#
# Output: output/figureS*.pdf / .png
#         output/supplementary_index.md
################################################################################

.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the 5.MsFigures project root (or R/).")
source(.setup[1])

## ---- what goes in the supplement --------------------------------------------

# `object` is the name the project script gives its finished plot.
SUPPLEMENTARY <- list(
  list(id      = "S1",
       stem    = "figureS1_allele_frequency_spectra",
       project = "mut", script = "04_figures.R", object = "fig",
       width   = 5.6, height = 6.6,
       title   = "Allele-frequency spectra against the neutral expectation",
       note    = paste(
         "Observed spectra for the three sequenced fractions against the",
         "neutral coalescent. The two spore fractions are all singletons,",
         "which is what neutrality predicts; the endpoint total fraction is",
         "not.")),

  list(id      = "S2",
       stem    = "figureS2_lineage_structure",
       project = "mut", script = "06_lineage_structure.R", object = "fig",
       width   = 9.2, height = 5.8,
       title   = "Clone genotypes and the lineage structure behind them",
       note    = paste(
         "Which clone carries which mutation, ordered so the clades are",
         "visible. The apparent mutual exclusivity of sinR and epsA-slrR is a",
         "property of the lineages, not of the mutations: they descend from",
         "five events, not twenty independent observations.")),

  list(id      = "S3",
       stem    = "figureS3_biofilm_by_lineage_group",
       project = "biofilm", script = "06_lineage_groups.R", object = "fig",
       width   = 5.8, height = 5.2,
       title   = "Biofilm by lineage group",
       note    = paste(
         "The three-group comparison: sinR against slrC/epsA-slrR against",
         "spore. This never uses the ancestor, so it is unaffected by the",
         "ancestor problem described for Figure 2.")),

  list(id      = "S4",
       stem    = "figureS4_ancestor_2023_diagnostic",
       project = "biofilm", script = "07_ancestor_2023.R", object = "p",
       width   = 7.6, height = 6,
       title   = "Diagnostic: the 2023 wild-type biofilm re-read",
       note    = paste(
         "Where the 2023 wild-type reading falls against the 2020 plate. The",
         "two runs read different wavelengths against different blanks and",
         "share no strain, so this is not a calibration and nothing in the",
         "analysis depends on it. Include only if the ancestor question is",
         "discussed."))
)

## ---- build ------------------------------------------------------------------

built <- list()

for (fg in SUPPLEMENTARY) {
  env <- load_project(fg$project, fg$script)

  if (!exists(fg$object, envir = env, inherits = FALSE)) {
    warning(sprintf("%s: '%s' not found in %s/R/%s -- skipped.",
                    fg$id, fg$object, PROJECTS[[fg$project]], fg$script),
            call. = FALSE)
    next
  }

  p <- get(fg$object, envir = env)
  if (!inherits(p, "ggplot")) {
    warning(sprintf("%s: '%s' is not a ggplot -- skipped.", fg$id, fg$object),
            call. = FALSE)
    next
  }

  save_figure(p, fg$stem, width = fg$width, height = fg$height)
  built[[length(built) + 1]] <- fg
  cat(sprintf("  %-3s %s\n", fg$id, fg$stem))
}

## ---- index ------------------------------------------------------------------

index <- c(
  "# Supplementary figures",
  "",
  "Built by `R/03_supplementary.R`. Each is the figure its own project draws;",
  "this project only collects them and gives them supplementary numbering.",
  "",
  "| | figure | source | file |",
  "| --- | --- | --- | --- |")

for (fg in built) {
  index <- c(index, sprintf("| **%s** | %s | `%s/R/%s` | `%s.pdf` |",
                            fg$id, fg$title, PROJECTS[[fg$project]],
                            fg$script, fg$stem))
}

index <- c(index, "", "## Notes", "")
for (fg in built) {
  index <- c(index, sprintf("**%s. %s.** %s", fg$id, fg$title, fg$note), "")
}

writeLines(index, file.path(OUT_DIR, "supplementary_index.md"))

cat(sprintf("Wrote %d supplementary figures and supplementary_index.md\n",
            length(built)))
