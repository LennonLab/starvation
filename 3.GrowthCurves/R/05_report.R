################################################################################
# 05_report.R
#
# Turn the model output into things that go into a manuscript.
#
# Input : output/model_comparison.csv, output/model_contrasts.csv,
#         output/posterior_summary.csv, data/comp_data_annotated.csv
# Output: output/table1_group_models.md   grouping vs. between-strain SD and LRT
#         output/table2_group_ratios.md   pairwise group ratios with CIs and P
#         output/table3_strain_estimates.md  per-strain posterior estimates
#         output/methods_stats.md         methods paragraph, numbers filled in
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

need <- file.path(OUT_DIR, c("model_comparison.csv", "model_contrasts.csv",
                             "posterior_summary.csv"))
if (!all(file.exists(need))) {
  stop("Missing model output -- run R/01_bayes_fitness.R and ",
       "R/04_group_models.R first.")
}

comparison <- read.csv(file.path(OUT_DIR, "model_comparison.csv"))
contrasts  <- read.csv(file.path(OUT_DIR, "model_contrasts.csv"))
psummary   <- read.csv(file.path(OUT_DIR, "posterior_summary.csv"),
                       check.names = FALSE)
fits       <- read.csv(file.path(DATA_DIR, "comp_data_annotated.csv"))

# How each parameter is named in print.
PARAM_LABEL <- c(umax = "µmax", A = "Yield (A)", L = "Lag (L)")
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
# Two numbers per parameter. w is the Akaike weight -- the multimodel quantity,
# summing to one within a set of strains. P is the Kenward-Roger F-test against
# the global mean on the same strains. Between-strain SDs and the chi-square
# LRT are in output/model_comparison.csv.

fmt_w <- function(x) ifelse(is.na(x), "\u2013", sprintf("%.2f", x))

rows <- unique(comparison[, c("scope", "model", "grouping", "n_groups")])
rows <- rows[order(match(rows$scope, unique(comparison$scope)),
                   match(rows$model, c("1", "2", "3", "4", "5a", "5"))), ]

t1 <- data.frame(
  Strains  = rows$scope,
  Model    = rows$model,
  Grouping = rows$grouping,
  Groups   = rows$n_groups,
  check.names = FALSE
)

for (p in PARAM_ORDER) {
  x <- comparison[comparison$parameter == p, ]
  i <- match(paste(rows$scope, rows$model), paste(x$scope, x$model))
  # A singular fit has no between-strain variation left to explain. Its test
  # and its Akaike weight are both reported as missing: with the strain effect
  # at zero, the model treats every curve as an independent observation, so
  # neither measures what the table is about.
  t1[[paste0(PARAM_LABEL[[p]], " w")]] <- fmt_w(ifelse(x$singular[i], NA,
                                                       x$akaike_weight[i]))
  pv <- ifelse(x$singular[i], NA, x$p_kr[i])
  t1[[paste0(PARAM_LABEL[[p]], " P")]] <- fmt_p(pv)
}

lag_singular <- all(comparison$singular[comparison$parameter == "L"])

write_md(
  md_table(t1),
  "table1_group_models.md",
  "Table 1. Do the strain groupings explain variation between strains?",
  paste0(
    "Linear mixed models of each log-transformed parameter, one per grouping, ",
    "with strain as a random effect and curves weighted by the inverse variance ",
    "of their fit. *w* is the Akaike weight: the relative support for each ",
    "grouping among those fit to the same strains, summing to one within each ",
    "set. *P* is a Kenward\u2013Roger F-test against the global-mean model on ",
    "the same strains, which corrects for the small number of strains. Sets of ",
    "strains are not comparable with one another. Model 3 (mutation vs. no ",
    "mutation) on all eleven strains compares the ancestor, M23 and M26 --- none ",
    "of which carries a mutation --- with the eight *sinR* and ",
    "*ywcC* clones. Among the ten evolved clones it is the same ",
    "partition as model 4, because the two clones without a mutation are the ",
    "two spore-fraction clones, so mutation and cell type cannot be separated ",
    "on these strains.",
    if (lag_singular) paste0(
      " For lag time the between-strain variance is estimated at zero in every ",
      "model: there is no strain-to-strain variation for a grouping to explain, ",
      "so neither a weight nor a test is reported.") else "")
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
  t2[[paste(PARAM_LABEL[[p]], "ratio")]] <- sprintf("%.2f (%.2f\u2013%.2f)",
                                                    x$ratio[i], x$lower[i], x$upper[i])
  t2[[paste(PARAM_LABEL[[p]], "P")]] <- fmt_p(x$p[i])
}

write_md(
  md_table(t2),
  "table2_group_ratios.md",
  "Table 2. Ratios of group means",
  paste0(
    "Ratio of the first group's mean to the second's on the measured scale, ",
    "with a 95% confidence interval, from the same mixed models as Table 1. ",
    "*P* uses Kenward\u2013Roger degrees of freedom and is Tukey-adjusted ",
    "within model 5a, which has three groups. A ratio of 1 means no difference.")
)

## ---- Table 3: per-strain estimates -----------------------------------------

abs_est <- psummary[psummary$scale == "absolute", ]
rel_est <- psummary[psummary$scale == "relative", ]
meta    <- strain_meta(STRIP_ORDER)

t3 <- data.frame(
  Strain   = short_label(STRIP_ORDER),
  Origin   = meta$origin,
  Cells    = ifelse(is.na(meta$cell), "–", meta$cell),
  # mutation_full already carries current names; see GENE_SYNONYMS in 00_setup
  Mutation = c(none = "none", sinR = "sinR", ywcC = "ywcC")[meta$mutation_full],
  check.names = FALSE
)

# Decimals chosen per parameter so the columns line up: the rates and yields
# sit below 1, the lag times above.
PARAM_DIGITS <- c(umax = 3, A = 3, L = 2)

for (p in PARAM_ORDER) {
  x <- abs_est[abs_est$parameter == p, ]
  i <- match(STRIP_ORDER, x$clone)
  f <- paste0("%.", PARAM_DIGITS[[p]], "f")
  t3[[PARAM_LABEL[[p]]]] <- sprintf(paste0(f, " [", f, ", ", f, "]"),
                                    x$`50%`[i], x$`2.5%`[i], x$`97.5%`[i])
}

r <- rel_est[rel_est$parameter == "umax", ]
i <- match(STRIP_ORDER, r$clone)
t3[["Relative µmax"]] <- ifelse(
  is.na(i), "1 (reference)",
  sprintf("%.2f [%.2f, %.2f]", r$`50%`[i], r$`2.5%`[i], r$`97.5%`[i]))

write_md(
  md_table(t3),
  "table3_strain_estimates.md",
  "Table 3. Posterior estimates for each strain",
  paste0(
    "Posterior median with a 95% equal-tailed credible interval, from six ",
    "growth curves per strain weighted by the precision of each curve fit. ",
    "Relative µmax is each strain's growth rate divided by the ",
    "ancestor's, computed draw by draw.")
)

## ---- Methods paragraph -----------------------------------------------------

analysed  <- fits[fits$clones %in% MODEL_ORDER, ]
n_curves  <- nrow(analysed)
n_strains <- length(MODEL_ORDER)
n_runs    <- length(unique(analysed$run))
pkg <- function(p) as.character(utils::packageVersion(p))

methods <- c(
  "### Growth curves",
  "",
  sprintf(paste(
    "Optical density at 600 nm was recorded every 15 min for each well on a",
    "microplate reader across %d plate runs. Each trace was fit with a modified",
    "Gompertz model, yielding a lower asymptote, the maximum yield (A), the",
    "maximum specific growth rate (µmax) and the lag time (L), each with",
    "an asymptotic standard error. Curves that failed visual inspection of the",
    "fit diagnostics were discarded and the six lowest-RMSE curves per strain",
    "retained, giving %d curves across %d strains (the ancestor and ten evolved",
    "clones). The six curves of a strain are technical replicates, not",
    "independent isolates."), n_runs, n_curves, n_strains),
  "",
  "### Strain-level estimates",
  "",
  paste(
    "Because the fitted parameters differ in precision from curve to curve,",
    "replicates were combined by inverse-variance weighting in a Bayesian model",
    "rather than averaged. For strain *j*, the log of the fitted value of curve",
    "*i* was modelled as Normal(µ[j], σ[j] / √w[i,j]), where w[i,j] is the",
    "inverse of the squared standard error of that curve's fit, scaled to mean",
    "one within a strain, and each strain has its own σ[j]. The model was fit",
    "in brms with Stan: µ[j] ~ Normal(0, 31.6) and log σ[j] ~ Normal(0, 5), the",
    "latter effectively flat over any plausible range. Four chains of 3,500",
    "iterations with 1,000 warm-up gave 10,000 draws; all R-hat values were",
    "below 1.01. Estimates are posterior medians with 95% equal-tailed credible",
    "intervals, and are what the figures show. Values relative to the ancestor",
    "are exp(µ[j]) divided by the ancestor's exp(µ), computed draw by draw so",
    "that the interval carries the uncertainty in both terms."),
  "",
  "### Group comparisons",
  "",
  paste(
    "Whether the strain groupings explain variation between strains was tested",
    "with linear mixed models (lme4): the log-transformed parameter with one",
    "fixed mean per group, a random intercept for strain, and curves weighted",
    "by the inverse variance of their fit. Candidate groupings were a single",
    "global mean; ancestor vs. evolved; mutation vs. no mutation; spore vs.",
    "total; and mutation identity (*sinR*, *ywcC* (BSU_38220), or none).",
    "Groupings that apply to a subset of the strains were compared against a",
    "global-mean model fit to that same subset, and no comparison is made",
    "across subsets."),
  "",
  paste(
    "Groupings were compared two ways. As a multimodel comparison, by AICc",
    "(MuMIn) from maximum-likelihood fits, reported as Akaike weights within",
    "each set of strains. And as a test, by a Kenward-Roger F-test against the",
    "global-mean model on REML fits (pbkrtest), which gives the reported *P*.",
    "Kenward-Roger rather than a chi-square likelihood-ratio test because there",
    "are only 8 to 11 strains, and the chi-square test is anticonservative with",
    "so few groups. Pairwise ratios of group means, with 95% confidence",
    "intervals and Kenward-Roger *P* values, were computed with emmeans,",
    "Tukey-adjusted where a grouping has three levels."),
  "",
  paste(
    "Two features of the strain set limit what these comparisons can show. The",
    "ancestor is a single strain, so in the ancestor vs. evolved comparison its",
    "group mean and its own effect are the same quantity. And among the ten",
    "evolved clones, the two without a detected mutation (M23 and M26) are also",
    "the two from the spore fraction, so mutation status and cell type are the",
    "same partition and cannot be separated."),
  "",
  "### Software",
  "",
  sprintf(paste(
    "Analyses were run in %s with brms %s, rstan %s, lme4 %s, pbkrtest %s,",
    "MuMIn %s, emmeans %s and ggplot2 %s. The brms strain-level model",
    "reproduces a reference JAGS fit's posterior medians to within 0.25%%; its",
    "95%% intervals are narrower than the reference's because the reference's",
    "Gamma(0.01, 0.01) prior on the precision was not vague at the scale the",
    "weights put it on. Code and data are in the growthCurves project."),
    R.version.string, pkg("brms"), pkg("rstan"), pkg("lme4"), pkg("pbkrtest"),
    pkg("MuMIn"), pkg("emmeans"), pkg("ggplot2")),
  ""
)

writeLines(c("## Methods: statistics", "", methods),
           file.path(OUT_DIR, "methods_stats.md"))
message("Wrote output/methods_stats.md")
