################################################################################
# 02_circos.R
#
# Figure 2A: mutations around the Bacillus subtilis genome, one ring per
# sequenced fraction.
#
# Input : data/mutations.csv
# Output: output/fig2a_circos.pdf / .png
#
# Rings run outside in, matching the published figure: endpoint total, endpoint
# spore, early spore. Each mutation is a tile coloured by how many clones carry
# it, and -- unlike the Circos original, which showed presence only -- drawn
# with a height proportional to that count, so a mutation in eight clones reads
# differently from a singleton at a glance.
#
# Tiles are widened to a fixed span so that a single base is visible on a 4.2
# Mb circle; the span is cosmetic and does not mean the mutation covers it.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the mutations project root (or R/).")
source(.setup[1])

if (!requireNamespace("circlize", quietly = TRUE)) {
  stop("Package 'circlize' is required. install.packages(\"circlize\")")
}

muts <- read.csv(file.path(DATA_DIR, "mutations.csv"), stringsAsFactors = FALSE)
muts$fraction <- factor(muts$fraction, levels = FRACTION_ORDER)

# Rings are drawn outermost first, so the fractions go in reverse order.
RING_ORDER <- rev(FRACTION_ORDER)

TILE_SPAN  <- 20000    # as in the published track: visible at genome scale
MAX_CARRIERS <- max(muts$carriers)

# Genes to label: anything that took more than one independent hit, plus
# anything carried by more than one clone. Between them these are every locus
# the figure shows in a colour other than singleton blue, so each coloured
# tile can be read off by name.
# Labels use gene_display -- the resolved current symbol -- so the figure reads
# slrC, not the legacy ywcC the workbook carries. See data/gene_name_overrides.csv.
label_genes <- muts %>%
  filter(!is.na(gene_display)) %>%
  group_by(gene_display) %>%
  summarise(hits = n(), max_carriers = max(carriers), .groups = "drop") %>%
  filter(hits > 1 | max_carriers > 1) %>%
  pull(gene_display)

# Genes taking more independent events than the genome-wide rate allows are
# named on the figure with their event count and set in bold, so the parallelism
# result can be read off the plot instead of needing a panel of its own.
# R/07_parallelism.R writes this, and run_all.R runs it first.
par_file <- file.path(OUT_DIR, "parallelism_tests.csv")
if (!file.exists(par_file)) {
  stop("output/parallelism_tests.csv not found -- run R/07_parallelism.R first.")
}
parallel_genes <- read.csv(par_file, stringsAsFactors = FALSE)

LABEL_CEX     <- 0.78   # gene names on the outer ring
RING_KEY_CEX  <- 0.70   # ring key in the centre
LABEL_CEX_PAR <- 0.86   # the parallel ones, a little larger again

draw_circos <- function() {
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 3,
                       cell.padding = c(0, 0, 0, 0),
                       track.margin = c(0.006, 0.006),
                       canvas.xlim = c(-1.12, 1.12),
                       canvas.ylim = c(-1.12, 1.12))

  circlize::circos.initialize(factors = GENOME, xlim = c(0, GENOME_BP))

  ## outermost: labels for genes hit more than once
  hits <- muts %>%
    filter(gene_display %in% label_genes) %>%
    group_by(gene_display) %>%
    slice_max(carriers, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gene = gene_display, position) %>%
    arrange(position)

  # Parallel genes carry their event count and are set bold italic (font 4)
  # against plain italic (font 3) for the rest.
  i_par <- match(hits$gene, parallel_genes$cds_gene)
  hits$label <- ifelse(is.na(i_par), hits$gene,
                       sprintf("%s (%d)", hits$gene,
                               parallel_genes$events[i_par]))
  hits$font <- ifelse(is.na(i_par), 3, 4)
  hits$cex  <- ifelse(is.na(i_par), LABEL_CEX, LABEL_CEX_PAR)
  hits$col  <- ifelse(is.na(i_par), "grey25", "black")

  circlize::circos.track(
    ylim = c(0, 1), track.height = 0.10, bg.border = NA,
    panel.fun = function(x, y) {
      if (!nrow(hits)) return(invisible(NULL))
      circlize::circos.text(hits$position, rep(0.45, nrow(hits)), hits$label,
                            cex = hits$cex, font = hits$font, col = hits$col,
                            facing = "clockwise", niceFacing = TRUE,
                            adj = c(0, 0.5))
    })

  ## genome axis
  circlize::circos.track(
    ylim = c(0, 1), track.height = 0.04, bg.border = NA,
    panel.fun = function(x, y) {
      circlize::circos.axis(
        h = "bottom",
        major.at = seq(0, 4e6, by = 1e6),
        labels = paste0(0:4, " Mb"),
        labels.cex = 0.7, major.tick.length = 0.5,
        labels.niceFacing = TRUE)
    })

  ## one ring per fraction, outermost first
  for (i in seq_along(RING_ORDER)) {
    frac <- RING_ORDER[i]
    d <- muts[muts$fraction == frac, ]

    circlize::circos.track(
      ylim = c(0, 1), track.height = 0.11,
      bg.col = RING_TINTS[i], bg.border = "grey70",
      panel.fun = function(x, y) {
        if (!nrow(d)) return(invisible(NULL))
        # Height on a square-root scale. On a linear scale a singleton would
        # be an eighth of the ring and invisible; sqrt keeps the ordering
        # while leaving the rarest mutations legible.
        h <- sqrt(d$carriers / MAX_CARRIERS)
        circlize::circos.rect(
          xleft   = pmax(d$position - TILE_SPAN / 2, 0),
          xright  = pmin(d$position + TILE_SPAN / 2, GENOME_BP),
          ybottom = 0, ytop = h,
          col = freq_colour(d$carriers, MAX_CARRIERS), border = NA)
      })
  }

  ## centre: ring key, each label carrying its ring's own tint as a swatch so
  ## there is no doubt which circle is which
  ys <- c(0.36, 0.27, 0.18)
  for (i in seq_along(RING_ORDER)) {
    lab <- FRACTIONS[[RING_ORDER[i]]]$label
    w   <- strwidth(lab, cex = RING_KEY_CEX)
    rect(-w / 2 - 0.075, ys[i] - 0.026, -w / 2 - 0.025, ys[i] + 0.026,
         col = RING_TINTS[i], border = "grey55", lwd = 0.6)
    text(0, ys[i], lab, cex = RING_KEY_CEX, col = "grey25")
  }
  # What the bold labels mean belongs in the caption, not in the middle of the
  # circle -- two lines of explanation here pushed the rest of the key into the
  # allele-frequency legend. 5.MsFigures writes it into figure1_caption_notes.md.
  text(0, 0.03, expression(italic("Bacillus subtilis")), cex = 0.92)
  text(0, -0.07, sprintf("%.1f Mb", GENOME_BP / 1e6), cex = 0.82, col = "grey30")

  ## allele-frequency legend, inside the circle below the centre
  n_stop <- 60
  counts <- seq(1, MAX_CARRIERS, length.out = n_stop)
  x0 <- -0.05; x1 <- 0.05; y0 <- -0.38; y1 <- -0.21
  ys_leg <- seq(y0, y1, length.out = n_stop + 1)
  for (i in seq_len(n_stop)) {
    rect(x0, ys_leg[i], x1, ys_leg[i + 1],
         col = freq_colour(counts[i], MAX_CARRIERS), border = NA)
  }
  rect(x0, y0, x1, y1, border = "grey40", lwd = 0.6)
  text(x1 + 0.02, y0, "Singleton", cex = 0.62, adj = c(0, 0.5), col = "grey25")
  text(x1 + 0.02, y1,
       sprintf("%d clones (%.1f%%)", MAX_CARRIERS,
               100 * MAX_CARRIERS / FRACTIONS[["S+NS1000"]]$n),
       cex = 0.62, adj = c(0, 0.5), col = "grey25")
  text(0, y1 + 0.055, "Allele frequency", cex = 0.68, col = "grey25")

  circlize::circos.clear()
}

for (dev in c("pdf", "png")) {
  f <- file.path(OUT_DIR, paste0("fig2a_circos.", dev))
  if (dev == "pdf") pdf(f, width = 7, height = 7)
  else png(f, width = 7, height = 7, units = "in", res = 200)
  par(mar = c(1, 1, 1, 1))
  draw_circos()
  dev.off()
}

message("Wrote output/fig2a_circos.pdf and .png")
