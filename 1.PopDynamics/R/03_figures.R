################################################################################
# 03_figures.R
#
# Output: fig_spore_vegetative.pdf  spore and non-spore abundance -- the panel
#                                   the manuscript uses
#         fig_total.pdf             total abundance (S + V), counted directly
#         fig_pop_dynamics.pdf      both stacked, the archived two-panel layout
#
# Input : data/pop_dynamics.csv, output/sigmoidal_curves.csv
#
# Total is counted; spore and non-spore abundances come from the latent-state
# model, since the difference between the two counts is not usable directly
# (see R/01_abundances.R).
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the populationDynamics project root (or R/).")
source(.setup[1])

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required. install.packages(\"patchwork\")")
}

pop    <- read.csv(file.path(DATA_DIR, "pop_dynamics.csv"))
curves <- read.csv(file.path(OUT_DIR, "sigmoidal_curves.csv"))

# Day zero cannot be drawn on a log time axis; the first reading is at 0.5 d.
pop <- pop[pop$time > 0, ]

XLIM <- c(0.4, 1400)
YLIM <- c(1e4, 1e9)
YLAB <- expression("Bacteria (CFU ml"^-1 * ")")

log_breaks <- function(range) 10^seq(log10(range[1]), log10(range[2]))

base_panel <- function() {
  list(
    # Plain numbers on the time axis, powers of ten on the abundance axis --
    # the way the published panel reads.
    scale_x_log10(limits = XLIM, breaks = c(1, 10, 100, 1000),
                  labels = c("1", "10", "100", "1000"),
                  sec.axis = dup_axis()),
    scale_y_log10(limits = YLIM, breaks = 10^(4:9),
                  labels = scales::label_log(), sec.axis = dup_axis()),
    annotation_logticks(sides = "", outside = FALSE),
    mytheme
  )
}

## ---- total ----------------------------------------------------------------

total_panel <- function(with_x_axis = TRUE) {
  p <- ggplot(pop, aes(time, total)) +
    geom_point(colour = "grey78", size = 2.2, alpha = 0.85) +
    geom_line(data = curves[curves$series == "total", ],
              aes(time, fitted), linetype = "dashed", linewidth = 0.9,
              colour = "black") +
    annotate("text", x = XLIM[2] * 0.95, y = 5e6, label = "Total (S + V)",
             hjust = 1, size = 4.2, colour = "grey20") +
    base_panel() +
    labs(x = "Time (d)", y = YLAB)

  if (!with_x_axis) {
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  }
  p
}

## ---- spore and non-spore ---------------------------------------------------

# Legend wording follows the manuscript draft: the figure says "Non-spore",
# the text calls the same cells vegetative.
CELL_COLOURS <- c(Spore = "grey35", `Non-spore` = "grey72")

cells <- pop %>%
  select(time, rep, Spore = spore_est, `Non-spore` = veg_est) %>%
  pivot_longer(c(Spore, `Non-spore`), names_to = "cell_type",
               values_to = "abundance") %>%
  mutate(cell_type = factor(cell_type, levels = names(CELL_COLOURS)))

cell_curves <- curves %>%
  filter(series %in% c("spore", "veg")) %>%
  mutate(cell_type = ifelse(series == "spore", "Spore", "Non-spore"))

cells_panel <- function() {
  ggplot(cells, aes(time, abundance, colour = cell_type)) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_line(data = cell_curves, aes(time, fitted, group = cell_type),
              linetype = "dashed", linewidth = 0.9, colour = "black",
              inherit.aes = FALSE) +
    scale_colour_manual(values = CELL_COLOURS, breaks = names(CELL_COLOURS)) +
    guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
    base_panel() +
    labs(x = "Time (d)", y = YLAB) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.83, 0.88),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key = element_blank())
}

## ---- write ------------------------------------------------------------------

save_fig <- function(plot, stem, width, height) {
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 200)
}

save_fig(cells_panel(), "fig_spore_vegetative", 6.4, 5)
save_fig(total_panel(), "fig_total", 6.4, 5)

# Both stacked, sharing one time axis -- the archived layout.
stacked <- patchwork::wrap_plots(
  total_panel(with_x_axis = FALSE) + labs(y = NULL),
  cells_panel() + labs(y = NULL),
  ncol = 1) &
  theme(plot.margin = margin(4, 6, 4, 6))

stacked <- patchwork::wrap_elements(stacked) +
  labs(tag = YLAB) +
  theme(plot.tag = element_text(size = 15, angle = 90),
        plot.tag.position = "left")

save_fig(stacked, "fig_pop_dynamics", 6, 8)

message("Wrote output/fig_spore_vegetative, fig_total and fig_pop_dynamics (.pdf/.png)")
