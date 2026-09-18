################################################################################
# 04_figures.R
#
# Figure 2B: the allele-frequency spectrum of each fraction, observed against
# the neutral expectation.
#
# Input : data/mutations.csv, output/null_spectrum.csv
# Output: output/fig2b_spectrum.pdf / .png
#
# One panel per sequenced fraction. Points are the observed number of
# polymorphisms at each unfolded allele frequency; the ribbon and line are the
# neutral coalescent expectation with its 95% interval across iterations.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the mutations project root (or R/).")
source(.setup[1])

muts <- read.csv(file.path(DATA_DIR, "mutations.csv"), stringsAsFactors = FALSE)
null <- read.csv(file.path(OUT_DIR, "null_spectrum.csv"), stringsAsFactors = FALSE)

panel_labels <- vapply(FRACTION_ORDER, function(f) FRACTIONS[[f]]$label,
                       character(1))

# Observed spectrum: how many mutations sit at each allele count.
observed <- muts %>%
  count(fraction, carriers, name = "polymorphisms") %>%
  mutate(frequency = carriers / vapply(fraction, function(f) FRACTIONS[[f]]$n,
                                       integer(1)))

# Only the low-frequency end carries anything; showing the whole 0-1 range
# would be almost all empty axis.
XMAX <- 0.25

to_panel <- function(d) {
  d$panel <- factor(panel_labels[d$fraction], levels = panel_labels)
  d
}

observed <- to_panel(observed)
null     <- to_panel(null[null$frequency <= XMAX, ])

fig <- ggplot() +
  geom_ribbon(data = null,
              aes(frequency, ymin = lower, ymax = upper),
              fill = "grey80", alpha = 0.7) +
  geom_line(data = null, aes(frequency, mean),
            colour = "grey45", linewidth = 0.6) +
  geom_point(data = observed[observed$frequency <= XMAX, ],
             aes(frequency, polymorphisms),
             size = 2.6, shape = 21, fill = "white", colour = "black",
             stroke = 0.8) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_x_continuous(limits = c(0, XMAX), expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Unfolded allele frequency", y = "Number of polymorphisms") +
  mytheme

ggsave(file.path(OUT_DIR, "fig2b_spectrum.pdf"), fig, width = 5.2, height = 6.4)
ggsave(file.path(OUT_DIR, "fig2b_spectrum.png"), fig, width = 5.2, height = 6.4,
       dpi = 200)

message("Wrote output/fig2b_spectrum.pdf and .png")
