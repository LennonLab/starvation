################################################################################
# run_all.R -- rebuild the analysis input, the figures, and the group models
#
#   Rscript run_all.R
#
# Takes about half a minute. R/90_gompertz_fits.R (raw plates -> curve fits) is
# the slow step upstream and is deliberately not included; see README.md.
################################################################################

for (script in c("R/91_build_comp_data.R",  # archived fits -> analysis input
                 "R/01_bayes_fitness.R",    # posterior draws for umax, A, L
                 "R/02_figures.R",          # ridgeline and mean +/- SD figures
                 "R/03_validate.R",         # check against the original JAGS run
                 "R/04_group_models.R",     # do the strain groupings explain it?
                 "R/05_report.R")) {        # manuscript tables and methods text
  message("\n===== ", script, " =====")
  source(script, echo = FALSE)
}

message("\nDone. Figures, tables and methods text are in output/.")
