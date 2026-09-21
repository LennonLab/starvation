################################################################################
# 04_group_models.R
#
# Do the strain groupings explain the growth parameters?
#
# One linear mixed model per candidate grouping, fit with lme4, compared within
# each set of strains two ways:
#
#   * a Kenward-Roger F-test against the global-mean model -- the P value
#     reported. KR rather than a chi-square likelihood-ratio test because there
#     are only 8-11 strains, and the chi-square LRT is anticonservative with so
#     few groups: for yield under model 3 it gives P = 0.0015 where KR gives
#     P = 0.0058.
#   * AICc and Akaike weights -- the multimodel comparison. Weights sum to one
#     within a set of strains and read as the relative support for each
#     grouping among those considered.
#
# plus pairwise contrasts between groups, as ratios of group means with 95% CIs
# and KR-adjusted P values (emmeans).
#
# The strain-level posteriors plotted in the figures come from brms
# (R/01_bayes_fitness.R). This script is the frequentist half: the part that
# produces P values.
#
# Input : data/comp_data_annotated.csv
# Output: output/model_comparison.csv  one row per grouping x parameter
#         output/model_contrasts.csv   pairwise group ratios
#
# Model, for curve i of strain j:
#
#     log(parameter)[i,j] ~ 0 + group[j] + (1 | strain),  weights = w[i,j]
#
# w[i,j] is the inverse variance of that curve's fitted parameter, scaled to
# mean one; lme4's prior weights make the residual variance sigma2_e / w[i,j].
# Cell-means coding, so each fixed effect is one group's mean on the log scale.
# AICc needs maximum-likelihood fits and KR needs REML fits, so each model is
# fit both ways.
#
# Candidate groupings (numbering follows the project notes):
#   1   global mean
#   2   ancestor vs evolved
#   3   mutation vs no mutation
#   4   spore vs total
#   5   sinR vs ywcC (BSU_38220)
#   5a  sinR vs ywcC vs no mutation
#
# Model 3 was previously skipped on the grounds that "every evolved isolate
# carries a mutation", which would make it identical to model 2. That was
# wrong: M23 and M26 carry no mutation (see strain_meta()). On all eleven
# strains model 3 is {ancestor, M23, M26} against the eight sinR/ywcC clones --
# a different partition from model 2 -- and it is fit here.
#
# Among the ten evolved clones, model 3 coincides with model 4: the two clones
# without a mutation are exactly the two spore-fraction clones. So mutation
# status and cell type cannot be separated on this strain set, and a result
# for either is a result for both. Only the six spore-fraction isolates set
# aside in comp_data_annotated.csv (S1, S6, S95 with mutations; S11, S22, S51
# without) would break that confound.
#
# Models 4 and 5 apply to subsets of the strains and are compared against the
# global mean within the same subset. Nothing is comparable across subsets.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

for (p in c("lme4", "pbkrtest", "MuMIn", "emmeans")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. install.packages(\"%s\")", p, p))
  }
}

## ---- data ------------------------------------------------------------------

fits <- read.csv(file.path(DATA_DIR, "comp_data_annotated.csv"),
                 stringsAsFactors = FALSE)
fits <- fits[fits$clones %in% MODEL_ORDER, ]

meta <- strain_meta(MODEL_ORDER)
meta$has_mutation <- ifelse(meta$mutated, "mutation", "none")

PARAMS <- list(umax = "umax", A = "A", L = "L")

## ---- which groupings, on which strains -------------------------------------

SCOPES <- list(
  list(
    id = "all (11)",
    strains = MODEL_ORDER,
    models = list(
      list(id = "1", label = "global mean",             grouping = NULL),
      list(id = "2", label = "ancestor vs evolved",     grouping = "origin"),
      list(id = "3", label = "mutation vs none",        grouping = "has_mutation")
    )
  ),
  list(
    id = "evolved (10)",
    strains = setdiff(MODEL_ORDER, "ancestor"),
    # model 3 is the same partition as model 4 here and is not repeated
    models = list(
      list(id = "1",  label = "global mean",            grouping = NULL),
      list(id = "4",  label = "spore vs total",         grouping = "cell"),
      list(id = "5a", label = "sinR vs ywcC vs none",   grouping = "mutation_full")
    )
  ),
  list(
    id = "sinR/ywcC (8)",
    strains = meta$clone[meta$mutation_full %in% c("sinR", "ywcC")],
    models = list(
      list(id = "1", label = "global mean",             grouping = NULL),
      list(id = "5", label = "sinR vs ywcC",            grouping = "mutation_full")
    )
  )
)

## ---- fitting ---------------------------------------------------------------

#' Fit one grouping, both ways.
#'
#' REML for the between-strain SD and the KR test (REML variance components
#' are less biased, and KR requires them); ML for AICc, since likelihoods of
#' models with different fixed effects are only comparable under ML.
fit_both <- function(dd, grouping) {
  form <- if (is.null(grouping)) y ~ 1 + (1 | strain) else y ~ 0 + grp + (1 | strain)
  fit <- function(reml) suppressMessages(suppressWarnings(
    lme4::lmer(form, data = dd, weights = w, REML = reml)))
  list(reml = fit(TRUE), ml = fit(FALSE))
}

sigma_strain <- function(m) sqrt(as.numeric(lme4::VarCorr(m)$strain))

#' KR F-test of a grouping against the global mean, or NA if it cannot be run.
kr_test <- function(big, small) {
  out <- tryCatch(pbkrtest::KRmodcomp(big, small)$test["Ftest", ],
                  error = function(e) NULL)
  if (is.null(out)) return(c(F = NA, ndf = NA, ddf = NA, p = NA))
  c(F = out$stat, ndf = out$ndf, ddf = out$ddf, p = out$p.value)
}

comparison <- list()
contrasts  <- list()

for (pname in names(PARAMS)) {
  column <- PARAMS[[pname]]
  se_col <- paste0(column, ".se")

  for (scope in SCOPES) {
    d <- fits[fits$clones %in% scope$strains, ]
    w <- 1 / d[[se_col]]^2

    base <- data.frame(y = log(d[[column]]), w = w / mean(w),
                       strain = factor(d$clones, levels = scope$strains))

    fitted <- list()
    for (mod in scope$models) {
      dd <- base
      if (!is.null(mod$grouping)) {
        dd$grp <- factor(meta[[mod$grouping]][match(d$clones, meta$clone)])
      }
      fitted[[mod$id]] <- c(fit_both(dd, mod$grouping), list(dd = dd, mod = mod))
    }

    null_fit <- fitted[["1"]]
    aicc   <- vapply(fitted, function(f) MuMIn::AICc(f$ml), numeric(1))
    weight <- exp(-0.5 * (aicc - min(aicc)))
    weight <- weight / sum(weight)

    for (id in names(fitted)) {
      f   <- fitted[[id]]
      mod <- f$mod
      kr  <- if (is.null(mod$grouping)) c(F = NA, ndf = NA, ddf = NA, p = NA)
             else kr_test(f$reml, null_fit$reml)
      lrt <- if (is.null(mod$grouping)) NA
             else anova(null_fit$ml, f$ml)$`Pr(>Chisq)`[2]

      comparison[[length(comparison) + 1]] <- data.frame(
        parameter     = pname,
        scope         = scope$id,
        model         = mod$id,
        grouping      = mod$label,
        n_groups      = if (is.null(mod$grouping)) 1L else nlevels(f$dd$grp),
        sigma_strain  = sigma_strain(f$reml),
        singular      = lme4::isSingular(f$reml),
        AICc          = aicc[[id]],
        delta_AICc    = aicc[[id]] - min(aicc),
        akaike_weight = weight[[id]],
        F             = kr[["F"]],
        ndf           = kr[["ndf"]],
        ddf           = kr[["ddf"]],
        p_kr          = kr[["p"]],
        p_lrt_chisq   = lrt,
        stringsAsFactors = FALSE)

      # Pairwise contrasts on the log scale, back-transformed to ratios of
      # group means. Tukey-adjusted where there are more than two groups.
      if (!is.null(mod$grouping)) {
        em <- suppressMessages(emmeans::emmeans(f$reml, ~ grp,
                                                lmer.df = "kenward-roger"))
        pw <- as.data.frame(summary(emmeans::contrast(em, "pairwise"),
                                    infer = c(TRUE, TRUE)))
        contrasts[[length(contrasts) + 1]] <- data.frame(
          parameter = pname,
          scope     = scope$id,
          model     = mod$id,
          contrast  = gsub(" - ", " / ", pw$contrast),
          ratio     = exp(pw$estimate),
          lower     = exp(pw$lower.CL),
          upper     = exp(pw$upper.CL),
          df        = pw$df,
          p         = pw$p.value,
          adjust    = if (nlevels(f$dd$grp) > 2) "tukey" else "none",
          stringsAsFactors = FALSE)
      }
    }
  }
}

comparison <- do.call(rbind, comparison)
contrasts  <- do.call(rbind, contrasts)

write.csv(comparison, file.path(OUT_DIR, "model_comparison.csv"), row.names = FALSE)
write.csv(contrasts,  file.path(OUT_DIR, "model_contrasts.csv"),  row.names = FALSE)

## ---- print -----------------------------------------------------------------

show <- comparison
show$sigma_strain  <- round(show$sigma_strain, 3)
show$akaike_weight <- round(show$akaike_weight, 2)
show$p_kr          <- signif(show$p_kr, 2)
show$ddf           <- round(show$ddf, 1)

for (pname in names(PARAMS)) {
  cat(sprintf("\n%s\n", pname))
  print(show[show$parameter == pname,
             c("scope", "model", "grouping", "n_groups", "sigma_strain",
               "akaike_weight", "ddf", "p_kr", "singular")],
        row.names = FALSE)
}

cat("\nNotes\n")
cat("  * P is a Kenward-Roger F-test against the global mean within the same\n",
    "   strains. The chi-square LRT is kept in model_comparison.csv for\n",
    "   reference; with 8-11 strains it runs smaller than it should.\n", sep = "")
cat("  * Model 3 on all eleven strains is {ancestor, M23, M26} against the\n",
    "   eight sinR/ywcC clones. Among the ten evolved clones it is identical to\n",
    "   model 4, so mutation and cell type cannot be separated here.\n", sep = "")
cat("  * In model 2 the ancestor group is one strain; its group mean and that\n",
    "   strain's effect are the same quantity.\n", sep = "")
if (any(comparison$singular)) {
  cat("  * Singular fits (between-strain SD estimated at zero): ",
      paste(unique(with(comparison[comparison$singular, ],
                        paste(parameter, scope, "model", model))), collapse = "; "),
      ". Their F-tests are unreliable.\n", sep = "")
}

message("\nWrote output/model_comparison.csv and model_contrasts.csv")
