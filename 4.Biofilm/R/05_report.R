################################################################################
# 05_report.R
#
# Turn the model output into things that go into a manuscript.
#
# Input : output/model_comparison.csv, output/model_contrasts.csv,
#         output/posterior_summary.csv, data/biofil.csv
# Output: output/table1_group_models.md   grouping vs. between-strain SD and LRT
#         output/table2_group_ratios.md   pairwise group ratios with CIs and P
#         output/table3_strain_estimates.md  per-strain posterior estimates
#         output/methods_stats.md         methods paragraph, numbers filled in
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the biofilm project root (or R/).")
source(.setup[1])

need <- file.path(OUT_DIR, c("model_comparison.csv", "model_contrasts.csv",
                             "posterior_summary.csv"))
if (!all(file.exists(need))) {
  stop("Missing model output -- run R/01_bayes_biofilm.R and ",
       "R/04_group_models.R first.")
}

comparison <- read.csv(file.path(OUT_DIR, "model_comparison.csv"))
contrasts  <- read.csv(file.path(OUT_DIR, "model_contrasts.csv"))
psummary   <- read.csv(file.path(OUT_DIR, "posterior_summary.csv"),
                       check.names = FALSE)
assay      <- read.csv(file.path(DATA_DIR, "biofil.csv"))

# How each parameter is named in print.
PARAM_LABEL <- c(biofilm = "Biofilm (OD550)")
PARAM_ORDER <- names(PARAM_LABEL)

fmt_p <- function(x) ifelse(is.na(x), "–",
                            ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))

# Minimal markdown table writer. Columns holding nothing but numbers are
# right-aligned; anything with words or interval brackets reads better left.
md_table <- function(df) {
  numeric_only <- vapply(df, function(col)
    all(grepl("^[-\u2013<>0-9.%]+$", as.character(col))), logical(1))
  align <- ifelse(numeric_only, "r", "l")
  sep <- ifelse(align == "l", ":---", "---:")
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("| ", paste(sep, collapse = " | "), " |"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

write_md <- function(lines, file, title, caption) {
  writeLines(c(paste("##", title), "", caption, "", lines, ""),
             file.path(OUT_DIR, file))
  message("Wrote output/", file)
}

## ---- Table 1: does the grouping explain between-strain variation? ----------
# w is the Akaike weight within each set of strains; P is the Kenward-Roger
# F-test against the global mean on the same strains. Between-strain SDs and
# the chi-square LRT are in output/model_comparison.csv.

fmt_w <- function(x) ifelse(is.na(x), "\u2013", sprintf("%.2f", x))

rows <- unique(comparison[, c("scope", "model", "grouping", "n_groups")])
rows <- rows[order(match(rows$scope, unique(comparison$scope)),
                   match(rows$model, c("1", "2", "3", "4", "5a", "5"))), ]

t1 <- data.frame(Strains = rows$scope, Model = rows$model,
                 Grouping = rows$grouping, Groups = rows$n_groups,
                 check.names = FALSE)

for (p in PARAM_ORDER) {
  x <- comparison[comparison$parameter == p, ]
  i <- match(paste(rows$scope, rows$model), paste(x$scope, x$model))
  t1[["w"]] <- fmt_w(ifelse(x$singular[i], NA, x$akaike_weight[i]))
  t1[["P"]] <- fmt_p(ifelse(x$singular[i], NA, x$p_kr[i]))
}

write_md(
  md_table(t1),
  "table1_group_models.md",
  "Table 1. Do the strain groupings explain variation in biofilm between strains?",
  paste0(
    "Linear mixed models of log blank-corrected OD550, one per grouping, with ",
    "strain as a random effect. *w* is the Akaike weight: the relative support ",
    "for each grouping among those fit to the same strains, summing to one ",
    "within each set. *P* is a Kenward\u2013Roger F-test against the ",
    "global-mean model on the same strains. Sets of strains are not comparable ",
    "with one another. Model 3 (mutation vs. no mutation) on all eleven strains ",
    "compares the reference well, m23 and m26 --- none of which carries a ",
    "mutation --- with the eight *sinR* and *ywcC* clones; among the ten evolved ",
    "clones it is the same partition as model 4. Models 2 and 3 on all eleven ",
    "strains include the reference well, recorded as *B. subtilis* 168 ",
    "\u03946 and probably not this experiment's ancestor; models 4, 5 and ",
    "5a do not use it.")
)

## ---- Table 2: pairwise group contrasts -------------------------------------

key <- unique(contrasts[, c("scope", "model", "contrast")])
key <- key[order(match(key$scope, unique(contrasts$scope)),
                 match(key$model, c("2", "3", "4", "5a", "5"))), ]

t2 <- data.frame(Strains = key$scope, Model = key$model,
                 Contrast = gsub(" / ", " : ", key$contrast),
                 check.names = FALSE)

for (p in PARAM_ORDER) {
  x <- contrasts[contrasts$parameter == p, ]
  i <- match(paste(key$scope, key$model, key$contrast),
             paste(x$scope, x$model, x$contrast))
  t2[["Ratio"]] <- sprintf("%.2f (%.2f\u2013%.2f)", x$ratio[i], x$lower[i], x$upper[i])
  t2[["P"]]     <- fmt_p(x$p[i])
}

write_md(
  md_table(t2),
  "table2_group_ratios.md",
  "Table 2. Ratios of group means for biofilm",
  paste0(
    "Ratio of the first group's mean to the second's on the measured scale, ",
    "with a 95% confidence interval, from the same mixed models as Table 1. ",
    "*P* uses Kenward\u2013Roger degrees of freedom and is Tukey-adjusted ",
    "within model 5a, which has three groups.")
)

## ---- Table 3: per-strain estimates -----------------------------------------

abs_est <- psummary[psummary$scale == "absolute", ]
rel_est <- psummary[psummary$scale == "relative", ]
meta    <- strain_meta(STRIP_ORDER)

t3 <- data.frame(
  Strain   = short_label(STRIP_ORDER),
  Origin   = meta$origin,
  Cells    = ifelse(is.na(meta$cell), "–", meta$cell),
  # mutation_full already carries current names; see MUTATION_LABEL in 00_setup
  Mutation = c(none = "none", sinR = "sinR", ywcC = "ywcC")[meta$mutation_full],
  check.names = FALSE
)

# Decimals chosen per parameter so the columns line up: the rates and yields
# sit below 1, the lag times above.
PARAM_DIGITS <- c(biofilm = 4)

for (p in PARAM_ORDER) {
  x <- abs_est[abs_est$parameter == p, ]
  i <- match(STRIP_ORDER, x$clone)
  f <- paste0("%.", PARAM_DIGITS[[p]], "f")
  t3[[PARAM_LABEL[[p]]]] <- sprintf(paste0(f, " [", f, ", ", f, "]"),
                                    x$`50%`[i], x$`2.5%`[i], x$`97.5%`[i])
}

r <- rel_est[rel_est$parameter == "biofilm", ]
i <- match(STRIP_ORDER, r$clone)
t3[["Relative biofilm"]] <- ifelse(
  is.na(i), "1 (reference)",
  sprintf("%.2f [%.2f, %.2f]", r$`50%`[i], r$`2.5%`[i], r$`97.5%`[i]))

write_md(
  md_table(t3),
  "table3_strain_estimates.md",
  "Table 3. Posterior estimates for each strain",
  paste0(
    "Posterior median with a 95% equal-tailed credible interval, from eight ",
    "replicate wells per strain. Relative biofilm is each strain's level ",
    "divided by the ancestor's, computed draw by draw.")
)

## ---- Methods paragraph -----------------------------------------------------

n_wells   <- nrow(assay)
n_strains <- length(MODEL_ORDER)
n_reps    <- N_REPS
pkg <- function(p) as.character(utils::packageVersion(p))

methods <- c(
  "### Biofilm assay",
  "",
  sprintf(paste(
    "Biofilm production was measured as crystal-violet retention read at 550 nm",
    "in a microtitre plate, with %d replicate wells per strain and %d wells in",
    "total across %d strains (the ancestor and ten evolved clones). Each",
    "reading was corrected by subtracting the mean of the plate's blank wells.",
    "The two spore-fraction clones were assayed as spores and the",
    "remainder as vegetative cells. Replicate wells are technical",
    "replicates of a single strain, not independent isolates."),
    n_reps, n_wells, n_strains),
  "",
  "### Strain-level estimates",
  "",
  paste(
    "Corrected optical densities span two orders of magnitude across strains,",
    "so they were modelled on the log scale: for strain *j*, the log reading",
    "of well *i* was Normal(µ[j], σ[j]), each strain with its own SD and no",
    "weighting, the wells being equally precise. The model was fit in brms",
    "with Stan, with µ[j] ~ Normal(0, 31.6) and log σ[j] ~ Normal(0, 5). Four",
    "chains of 11,000 iterations with 1,000 warm-up gave 40,000 draws; all",
    "R-hat values were below 1.01. Estimates are posterior medians with 95%",
    "equal-tailed credible intervals, and are what the figures show."),
  "",
  "### Group comparisons",
  "",
  paste(
    "Whether the strain groupings explain variation between strains was tested",
    "with linear mixed models (lme4): log OD with one fixed mean per group and",
    "a random intercept for strain. Candidate groupings were a single global",
    "mean; ancestor vs. evolved; mutation vs. no mutation; spore vs. total; and",
    "mutation identity (*sinR*, *ywcC* (BSU_38220), or none). Groupings that",
    "apply to a subset of strains were compared against a global-mean model on",
    "the same subset, and no comparison is made across subsets. Groupings were",
    "compared by AICc (MuMIn) from maximum-likelihood fits, reported as Akaike",
    "weights, and tested by a Kenward-Roger F-test against the global-mean",
    "model on REML fits (pbkrtest), which gives the reported *P*. Pairwise",
    "ratios of group means with 95% confidence intervals and Kenward-Roger *P*",
    "values were computed with emmeans, Tukey-adjusted where a grouping has",
    "three levels. Throughout, the strain is the unit of replication: the",
    "eight wells of a strain are technical replicates and are not treated as",
    "independent."),
  "",
  paste(
    "The reference well is recorded as \"Ancestor\" in the assay file but",
    "named *B. subtilis* 168 Δ6 in the plate record, and is probably not this",
    "experiment's ancestor. Every comparison that involves it --- ancestor vs.",
    "evolved, and mutation vs. no mutation on all eleven strains --- is",
    "provisional. Comparisons among the evolved clones do not use it."),
  "",
  "### Software",
  "",
  sprintf(paste(
    "Analyses were run in %s with brms %s, rstan %s, lme4 %s, pbkrtest %s,",
    "MuMIn %s, emmeans %s and ggplot2 %s. For the strain-level model the",
    "posterior of µ[j] can be computed exactly on a grid under the same priors,",
    "and the brms fit reproduces it to within 1%% on medians and 3%% on the",
    "bounds of the 95%% intervals. Code and data are in the biofilm project."),
    R.version.string, pkg("brms"), pkg("rstan"), pkg("lme4"), pkg("pbkrtest"),
    pkg("MuMIn"), pkg("emmeans"), pkg("ggplot2")),
  ""
)

writeLines(c("## Methods: statistics", "", methods),
           file.path(OUT_DIR, "methods_stats.md"))
message("Wrote output/methods_stats.md")
