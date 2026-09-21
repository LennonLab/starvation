################################################################################
# run_all.R -- build both manuscript figures.
#
# The four analysis projects must have been run first: the panels are built
# from their saved posteriors and tables, not recomputed here.
################################################################################

# 05 before 04: 05 writes the chosen Figure 2, and 04 stages what exists for
# upload, so the order decides which Figure 2 reaches Overleaf.
scripts <- c("R/01_figure1.R", "R/02_figure2.R", "R/03_supplementary.R",
             "R/05_figure2_variants.R", "R/04_tables.R")

for (s in scripts) {
  cat("\n=====", s, "=====\n")
  source(s, echo = FALSE)
}

cat("\nDone. Figures are in output/.\n")
