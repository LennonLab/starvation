################################################################################
# run_all.R -- mutation tables, the two Figure 2 panels, and the results text
#
#   Rscript run_all.R
#
# Takes a few seconds.
################################################################################

# 07 runs before 02: the circos marks the genes 07 finds, so its table has to
# exist first. Nothing else depends on the order.
for (script in c("R/01_mutations.R",    # variant matrices -> tidy mutation table
                 "R/07_parallelism.R",  # genes hit more often than chance
                 "R/02_circos.R",       # Fig 2a: mutations around the genome
                 "R/03_null_model.R",   # neutral coalescent expectation
                 "R/04_figures.R",      # Fig 2b: allele-frequency spectra
                 "R/06_lineage_structure.R",  # co-occurrence and clone backgrounds
                 "R/05_report.R")) {    # tables and results text
  message("\n===== ", script, " =====")
  source(script, echo = FALSE)
}

message("\nDone. Figures, tables and text are in output/.")
