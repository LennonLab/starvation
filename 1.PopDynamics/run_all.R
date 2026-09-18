################################################################################
# run_all.R -- abundances, curve fits, the figure, and the manuscript tables
#
#   Rscript run_all.R
#
# Takes a few seconds. R/90_latent_states_jags.R, which infers the spore and
# vegetative states, is the slow step upstream and needs JAGS; its archived
# output is in data/latent_states/. See README.md.
################################################################################

for (script in c("R/01_abundances.R",       # counts -> CFU/mL, + latent states
                 "R/02_sigmoidal_fits.R",   # abundance series -> curve fits
                 "R/03_figures.R",          # the two-panel figure
                 "R/04_report.R",          # descriptors, tables, methods text
                 "R/05_results_text.R")) {  # the results paragraphs
  message("\n===== ", script, " =====")
  source(script, echo = FALSE)
}

message("\nDone. Figures, tables and text are in output/.")
