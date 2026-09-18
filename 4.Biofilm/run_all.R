################################################################################
# run_all.R -- posteriors, figures, group models, and manuscript tables
#
#   Rscript run_all.R
#
# Takes a few seconds.
################################################################################

for (script in c("R/01_bayes_biofilm.R",   # assay -> posterior draws
                 "R/02_figures.R",         # posterior draws -> figures
                 "R/03_validate.R",        # check against the exact posterior
                 "R/04_group_models.R",    # do the strain groupings explain it?
                 "R/06_lineage_groups.R",  # the early draft's three-group comparison
                 "R/07_ancestor_2023.R",   # diagnostic: the 2023 wild-type re-read
                 "R/05_report.R")) {       # manuscript tables and methods text
  message("\n===== ", script, " =====")
  source(script, echo = FALSE)
}

message("\nDone. Figures, tables and methods text are in output/.")
