################################################################################
# 01_figure1.R
#
# Figure 1. What the population did over 1000 days, and what changed in its
# genomes.
#
#   A  spore and non-spore abundance, with the counted total over them
#                                             1.PopDynamics/R/03_figures.R
#   B  mutations around the genome            2.Mutations/R/02_circos.R
#
# Output: output/figure1_dynamics_mutations.pdf / .png
#         output/figure1_caption_notes.md
#
# The panels are the projects' own objects, not copies -- see load_project()
# in R/00_setup.R.
#
# The parallelism result used to be a third panel here. It is two numbers, and
# two numbers do not need a panel: the circos names the genes with their event
# counts in bold, and the figures below are written into the caption notes.
################################################################################

.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the manuscript-figures project root (or R/).")
source(.setup[1])

## ---- A: population dynamics -------------------------------------------------

pop_env <- load_project("pop", "03_figures.R")
require_objects(pop_env, "cells_panel", "CELL_COLOURS", "curves",
                what = "1.PopDynamics")

# The merged figure gives this panel less width than the standalone version, so
# the two greys are pushed apart before the panel is built. CELL_COLOURS is
# read inside cells_panel(), so assigning it here is enough.
# The total joins the same legend rather than getting its own floating label.
# cells_panel() reads CELL_COLOURS when it is called and hands its values and
# names straight to scale_colour_manual, so adding a third entry here adds a
# third key. The `cells` data frame was built when the script was sourced and
# still has two levels, so no point is drawn for it -- only the key.
TOTAL_LABEL <- "Total"
pop_env$CELL_COLOURS <- c(CELL_GREYS, setNames(TOTAL_COLOUR, TOTAL_LABEL))

total_curve <- pop_env$curves[pop_env$curves$series == "total", ]
stopifnot(nrow(total_curve) > 0)

panel_a <- pop_env$cells_panel() +
  geom_line(data = total_curve, aes(time, fitted, colour = TOTAL_LABEL),
            inherit.aes = FALSE, linewidth = 1) +
  # Per-key glyphs: points for the two cell types, a line for the total.
  # override.aes takes one value per key, so each is given a vector.
  guides(colour = guide_legend(override.aes = list(
    shape     = c(16, 16, NA),
    size      = c(3.4, 3.4, 0),
    linewidth = c(0, 0, 1),
    alpha     = 1))) +
  theme(legend.position.inside = c(0.78, 0.88),
        legend.text = element_text(size = 9.5),
        legend.key.spacing.y = unit(1, "pt"),
        plot.margin = margin(4, 8, 4, 4))

## ---- B: circos --------------------------------------------------------------

mut_env <- load_project("mut", "02_circos.R")
require_objects(mut_env, "draw_circos", "parallel_genes", what = "2.Mutations")

# circlize draws with base graphics. cowplot::as_grob replays the drawing onto
# a grid device, so the panel stays vector rather than being rasterised in.
panel_b <- cowplot::as_grob(function() {
  op <- par(mar = c(0, 0, 0, 0)); on.exit(par(op), add = TRUE)
  mut_env$draw_circos()
})

## ---- assemble ---------------------------------------------------------------

figure1 <- (panel_a | patchwork::wrap_elements(full = panel_b)) +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  TAG_THEME

save_figure(figure1, "figure1_dynamics_mutations", width = 12.6, height = 5.8)

## ---- caption notes ----------------------------------------------------------

# The parallelism numbers, written out so the caption can quote them without
# anyone reading them off a figure that no longer shows them.
N_PERM <- local({
  src <- readLines(file.path(project_dir("mut"), "R", "07_parallelism.R"),
                   warn = FALSE)
  hit <- regmatches(src, regexpr("N_PERM\\s*<-\\s*[0-9]+", src))
  if (!length(hit)) stop("Could not read N_PERM from 2.Mutations/R/07_parallelism.R")
  as.numeric(sub(".*<-\\s*", "", hit[1]))
})

pg <- mut_env$parallel_genes
pg <- pg[order(-pg$events), ]

fmt_p <- function(p) ifelse(p <= 0, sprintf("P < %s", format(1 / N_PERM)),
                            sprintf("P = %s", signif(p, 2)))

sentences <- sprintf(
  "*%s* took %d independent events against %.4f expected for a gene of its length (%.0f-fold; %s)",
  pg$cds_gene, pg$events, pg$expected, pg$events / pg$expected,
  fmt_p(pg$p_permutation))

notes <- c(
  "# Figure 1 - notes for the caption",
  "",
  "**A.** Spore and non-spore abundance over 1000 days. Points are replicate",
  "counts; the black dashed curves are the fitted sigmoids. Spore and non-spore",
  "are estimated by the latent-state model, since the difference between the",
  "two counts is not usable directly. The dark red curve is the total (S + V),",
  "which is counted rather than modelled.",
  "",
  "**B.** Mutations around the *Bacillus subtilis* genome, one ring per",
  "sequenced fraction: endpoint total (outer), endpoint spore, early spore",
  "(inner). Each mutation is a tile whose colour and height give the number of",
  "clones carrying it. Named genes are those hit more than once or carried by",
  "more than one clone.",
  "",
  "Genes in **bold** carry more independent mutational events than the",
  "genome-wide rate allows, with the event count in parentheses:",
  paste0(paste(sentences, collapse = "; "), "."),
  "",
  sprintf(paste(
    "P values are min-P permutation family-wise values over %s scatterings of",
    "the events across all coding sequence, each gene judged against its own",
    "length - the appropriate correction when the gene was nominated by the",
    "data rather than chosen in advance. These are the two genes phenotyped",
    "in Figure 2."), format(N_PERM, big.mark = ",")),
  "")

writeLines(notes, file.path(OUT_DIR, "figure1_caption_notes.md"))

cat("Wrote output/figure1_dynamics_mutations.pdf / .png and figure1_caption_notes.md\n")
