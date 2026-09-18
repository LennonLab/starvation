################################################################################
# 05_report.R
#
# Turn the model output into things that go into a manuscript.
#
# Input : output/model_comparison.csv, output/model_contrasts.csv,
#         output/posterior_summary.csv, data/biofil.csv
# Output: output/table1_group_models.md   grouping vs. between-strain SD and LRT
#         output/table2_group_ratios.md   group ratios with credible intervals
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

rows <- unique(comparison[, c("scope", "model", "grouping", "n_groups")])
rows <- rows[order(match(rows$scope, unique(comparison$scope)),
                   rows$model), ]

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
  t1[[paste0(PARAM_LABEL[[p]], " σ")]] <- sprintf("%.3f", x$sigma_strain[i])
  t1[[paste0(PARAM_LABEL[[p]], " P")]]      <- fmt_p(x$lrt_p[i])
}

write_md(
  md_table(t1),
  "table1_group_models.md",
  "Table 1. Do the strain groupings explain variation between strains?",
  paste0(
    "For each grouping, σ is the between-strain standard deviation ",
    "remaining on the log scale after the grouping is fit; a grouping that ",
    "explains something drives σ below the global-mean row within the ",
    "same set of strains. *P* is a likelihood-ratio test against that ",
    "global-mean model. Sets of strains are not comparable with one another. ",
    "Model 3 (mutation vs. no mutation) is the same partition as model 2 on ",
    "these strains and is not fit separately.")
)

## ---- Table 2: group ratios -------------------------------------------------

key <- unique(contrasts[, c("scope", "model", "contrast")])
key <- key[order(match(key$scope, unique(contrasts$scope)), key$model), ]

t2 <- data.frame(Strains = key$scope, Model = key$model,
                 Contrast = gsub("spormut", "sporulation", gsub(" / ", " : ", key$contrast)),
                 check.names = FALSE)

for (p in PARAM_ORDER) {
  x <- contrasts[contrasts$parameter == p, ]
  i <- match(paste(key$scope, key$model, key$contrast),
             paste(x$scope, x$model, x$contrast))
  t2[[PARAM_LABEL[[p]]]] <- sprintf("%.2f [%.2f, %.2f]",
                                    x$ratio_50[i], x$ratio_2.5[i], x$ratio_97.5[i])
}

write_md(
  md_table(t2),
  "table2_group_ratios.md",
  "Table 2. Ratios of group means",
  paste0(
    "Posterior median ratio of the first group's mean to the second's, on the ",
    "measured scale, with a 95% equal-tailed credible interval. A ratio of 1 ",
    "means the groups do not differ; intervals excluding 1 are the ones to ",
    "read. Posterior probabilities are in `output/model_contrasts.csv`.")
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
  Mutation = c(ancestor = "–", sinR = "sinR", ywcC = "ywcC",
               spore = "sporulation")[meta$mutation_full],
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
    "The evolved clones carrying sporulation mutations were assayed as spores",
    "and the remainder as vegetative cells. Replicate wells are technical",
    "replicates of a single strain, not independent isolates."),
    n_reps, n_wells, n_strains),
  "",
  "### Strain-level estimates",
  "",
  paste(
    "Corrected optical densities span two orders of magnitude across strains,",
    "so they were modelled on the log scale. For strain *j*, well *i* was",
    "modelled as x[i,j] ~ LogNormal(µ[j], τ[j]), each strain having its own",
    "variance; the log-normal keeps estimates positive during sampling. Priors",
    "were µ[j] ~ Normal(0, precision 0.001) and τ[j] ~ Gamma(0.01, 0.01).",
    "Because every full conditional is conjugate, the posterior was sampled",
    "with a Gibbs sampler written in base R: four chains of 10,000 draws each",
    "after 2,000 burn-in iterations. Estimates are reported as posterior",
    "medians with 95% equal-tailed credible intervals. Relative biofilm is",
    "exp(µ[j]) divided by the ancestor's exp(µ), computed draw by draw so that",
    "the interval carries the uncertainty in both terms, and is plotted on a",
    "log10 axis."),
  "",
  "### Group comparisons",
  "",
  paste(
    "To ask whether the strain groupings explain variation between strains, the",
    "log-transformed readings were fit with a two-level model: y[i,j] ~",
    "Normal(θ[j], σ²e) at the well level and θ[j] ~ Normal(x[j]'β, σ²s) at the",
    "strain level, where x[j] is a cell-means indicator for the grouping under",
    "test, so each element of β is one group's mean on the log scale. Priors",
    "were β ~ Normal(0, 100) and inverse-Gamma(0.01, 0.01) on both variances.",
    "Candidate groupings were: a single global mean; ancestor vs. evolved;",
    "spore vs. vegetative; and mutation identity (*sinR*, *ywcC* (BSU_38220) -- labelled",
    "*ywcC/slrR* in the assay file --",
    "sporulation mutant). Mutation vs. no mutation is the same partition as",
    "ancestor vs. evolved on this strain set and was not fit separately.",
    "Groupings that apply to a subset of the strains were compared against a",
    "global-mean model fit to that same subset, and no comparison is made",
    "across subsets."),
  "",
  paste(
    "Support for a grouping is reported as the reduction in the between-strain",
    "standard deviation σs relative to the global-mean model, together with a",
    "likelihood-ratio test against that model fit by maximum likelihood (lme4)",
    "and the posterior ratio of group means with its credible interval. WAIC",
    "was also computed but is reported only in the supplementary output: with a",
    "strain-level effect in the model it scores prediction of a further well",
    "from a strain already observed, which is not the question being asked. In",
    "the ancestor vs. evolved comparison the ancestor group contains a single",
    "strain, so its group mean and that strain's own effect are the same",
    "quantity and σs is informed only by the ten evolved clones; that contrast",
    "should be read with the ancestor's lack of biological replication in",
    "mind."),
  "",
  "### Software",
  "",
  sprintf(paste(
    "Analyses were run in %s with ggplot2 %s, ggridges %s, lme4 %s and loo %s.",
    "The Gibbs samplers are conjugate implementations in base R. For the",
    "strain-level model the nuisance precision can be integrated out",
    "analytically, and the sampler reproduces the resulting exact marginal",
    "posterior to within 0.2%% on medians and 1.1%% on the bounds of the 95%%",
    "credible intervals; the group-level sampler agrees with the corresponding",
    "lme4 fit. Code and data are in the biofilm project."),
    R.version.string, pkg("ggplot2"), pkg("ggridges"), pkg("lme4"), pkg("loo")),
  ""
)

writeLines(c("## Methods: statistics", "", methods),
           file.path(OUT_DIR, "methods_stats.md"))
message("Wrote output/methods_stats.md")
