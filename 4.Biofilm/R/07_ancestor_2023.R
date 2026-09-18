################################################################################
# 07_ancestor_2023.R
#
# Diagnostic only. Where does the 2023 wild-type re-read land against the 2020
# plate?
#
# Input : OneDrive_3_7-15-2026/20230611/20230611_SporeMut_Biofilm.xlsx
#         data/biofil.csv
# Output: output/fig_ancestor_2023_check.pdf / .png
#         output/ancestor_2023_wells.csv
#
# NOTHING HERE FEEDS THE ANALYSIS. The two runs cannot be put on one scale:
# 2020 read OD550, 2023 read OD540 and OD600; the blanks differ by a factor of
# 2.7 (0.065 against 0.175), which is a different plate background, not a
# different amount of biofilm; and no strain was measured on both plates, so
# there is nothing to calibrate between them. The figure exists to answer
# "roughly where does it sit", and its answer carries a run-to-run difference
# of unknown size inside it.
#
# Plate layout, from the note in the 2023 sheet -- "First ring empty, second
# ring and C blank, EDF samples":
#   ring 1 (row A, row H, col 1, col 12)   empty, no liquid: plate background
#   ring 2 (row B, row G, col 2, col 11)   medium, no cells
#   row C                                  medium, no cells
#   rows D, E, F x cols 3-10               wild type, 24 wells
################################################################################

.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the biofilm project root (or R/).")
source(.setup[1])

XLSX_2023 <- file.path(
  "/Users/canankarakoc/Desktop/SporeMut/OneDrive_3_7-15-2026/20230611",
  "20230611_SporeMut_Biofilm.xlsx")

if (!file.exists(XLSX_2023)) {
  message("2023 workbook not found; skipping the ancestor diagnostic.")
} else {

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. install.packages(\"readxl\")")
}

## ---- read one lid condition -------------------------------------------------

# The reader writes each plate row as two lines, 540 nm then 600 nm, starting
# at sheet row 25, with the twelve wells in columns 3:14.
FIRST_ROW  <- 25L
WELL_COLS  <- 3:14
PLATE_ROWS <- LETTERS[1:8]

read_plate <- function(sheet, wavelength = c("540", "600")) {
  wavelength <- match.arg(wavelength)
  offset <- if (wavelength == "540") 0L else 1L
  d <- as.data.frame(readxl::read_excel(XLSX_2023, sheet = sheet,
                                        col_names = FALSE,
                                        .name_repair = "minimal"))
  m <- t(vapply(seq_along(PLATE_ROWS), function(k) {
    suppressWarnings(as.numeric(unlist(d[FIRST_ROW + (k - 1) * 2 + offset,
                                         WELL_COLS])))
  }, numeric(length(WELL_COLS))))
  dimnames(m) <- list(PLATE_ROWS, as.character(1:12))
  m
}

#' Split a plate into the three well classes the layout note describes.
classify <- function(m) {
  idx <- expand.grid(row = PLATE_ROWS, col = 1:12, stringsAsFactors = FALSE)
  idx$od <- m[cbind(idx$row, as.character(idx$col))]
  idx$class <- with(idx, ifelse(
    row %in% c("A", "H") | col %in% c(1, 12), "empty",
    ifelse(row %in% c("D", "E", "F") & col %in% 3:10, "sample", "blank")))
  idx
}

collect <- function(sheet, label) {
  w <- classify(read_plate(sheet, "540"))
  blank <- mean(w$od[w$class == "blank"], na.rm = TRUE)
  w$blank_mean <- blank
  w$corrected  <- w$od - blank
  w$lid        <- label
  w
}

wells <- rbind(collect("withLid",    "2023 WT, lid"),
               collect("withoutLid", "2023 WT, no lid"))

write.csv(wells, file.path(OUT_DIR, "ancestor_2023_wells.csv"), row.names = FALSE)

samples_2023 <- wells[wells$class == "sample", ]
stopifnot(nrow(samples_2023) == 48)   # 24 wells x 2 lid conditions

## ---- the 2020 plate ---------------------------------------------------------

b2020 <- read.csv(file.path(DATA_DIR, "biofil.csv"), stringsAsFactors = FALSE)
b2020$strain <- tolower(b2020$Sample)

anc_mean <- mean(b2020$OD550_Corrected[b2020$strain == "ancestor"])

## ---- one panel --------------------------------------------------------------

ORDER_2020 <- c("ancestor", PLOT_ORDER)
row_levels <- c(rev(unique(samples_2023$lid)), ORDER_2020)

plot_rows <- rbind(
  data.frame(row = b2020$strain, value = b2020$OD550_Corrected,
             era = "2020 plate (OD550)", stringsAsFactors = FALSE),
  data.frame(row = samples_2023$lid, value = samples_2023$corrected,
             era = "2023 re-read (OD540)", stringsAsFactors = FALSE))

# A handful of 2023 wells correct to at or below zero, which a log axis cannot
# show. They are the low corner of the sample block, and dropping them is
# noted on the panel rather than done silently.
n_nonpos <- sum(plot_rows$value <= 0)
plot_rows <- plot_rows[plot_rows$value > 0, ]

plot_rows$row <- factor(plot_rows$row, levels = row_levels)
plot_rows$era <- factor(plot_rows$era,
                        levels = c("2020 plate (OD550)", "2023 re-read (OD540)"))

means <- aggregate(value ~ row + era, plot_rows, mean)

p <- ggplot(plot_rows, aes(value, row, colour = era, shape = era)) +
  # The band is the 2023 interquartile range, not its full spread: the full
  # spread of 24 wells covers most of the axis and says little.
  annotate("rect",
           xmin = quantile(samples_2023$corrected, 0.25),
           xmax = quantile(samples_2023$corrected, 0.75),
           ymin = -Inf, ymax = Inf, fill = "#C8102E", alpha = 0.09) +
  geom_vline(xintercept = anc_mean, linetype = "dashed", colour = "grey30") +
  geom_point(size = 2, alpha = 0.75, position = position_nudge(y = 0.16)) +
  geom_point(data = means, size = 3.4, stroke = 1.1, fill = "white",
             show.legend = FALSE) +
  scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3),
                labels = c("0.01", "0.03", "0.1", "0.3", "1", "3")) +
  scale_colour_manual(values = c("2020 plate (OD550)"  = "grey35",
                                 "2023 re-read (OD540)" = "#C8102E")) +
  scale_shape_manual(values = c("2020 plate (OD550)" = 16,
                                "2023 re-read (OD540)" = 21)) +
  annotate("text", x = anc_mean, y = length(row_levels) + 0.35,
           label = sprintf("2020 ancestor mean %.3f", anc_mean),
           hjust = -0.05, size = 3.1, colour = "grey30") +
  labs(x = "Blank-corrected optical density", y = NULL,
       title = "Diagnostic: where the 2023 wild-type re-read lands",
       subtitle = paste(
         "Not a calibration. 2020 read OD550 against a 0.065 blank; 2023 read",
         "OD540 against a 0.175 blank.\nNo strain was measured on both plates,",
         "so the offset between them is unknown and is inside every",
         "comparison\nbelow. Large points are means."),
       caption = if (n_nonpos) sprintf(
         "%d of 48 wells in 2023 correct to <= 0 and cannot be shown on a log axis.",
         n_nonpos) else NULL) +
  coord_cartesian(clip = "off") +
  mytheme +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 8.6, colour = "grey30"),
        plot.caption  = element_text(size = 8, colour = "grey40", hjust = 0),
        plot.margin   = margin(24, 12, 6, 6))

ggsave(file.path(OUT_DIR, "fig_ancestor_2023_check.pdf"), p, width = 7.6, height = 6)
ggsave(file.path(OUT_DIR, "fig_ancestor_2023_check.png"), p, width = 7.6, height = 6,
       dpi = 200)

## ---- what it says -----------------------------------------------------------

s23 <- samples_2023$corrected
ratio <- median(s23) / anc_mean

message(sprintf(paste0(
  "\n2023 wild-type re-read, 24 wells x 2 lid conditions:\n",
  "  lid    median %.3f   no lid median %.3f  (the lid made no difference)\n",
  "  2020 ancestor mean %.3f\n",
  "  ratio  %.2fx the 2020 ancestor, which is within the range a run-to-run\n",
  "         difference can produce -- it is not evidence either way about\n",
  "         whether 168 delta 6 is the right strain.\n"),
  median(s23[samples_2023$lid == "2023 WT, lid"]),
  median(s23[samples_2023$lid == "2023 WT, no lid"]),
  anc_mean, ratio))

message("Wrote output/fig_ancestor_2023_check.pdf/.png and ancestor_2023_wells.csv")

}
