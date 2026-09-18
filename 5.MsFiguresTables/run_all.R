################################################################################
# run_all.R -- build both manuscript figures.
#
# The four analysis projects must have been run first: the panels are built
# from their saved posteriors and tables, not recomputed here.
################################################################################

scripts <- c("R/01_figure1.R", "R/02_figure2.R", "R/03_supplementary.R",
             "R/04_tables.R")

for (s in scripts) {
  cat("\n=====", s, "=====\n")
  source(s, echo = FALSE)
}

cat("\nDone. Figures are in output/.\n")
