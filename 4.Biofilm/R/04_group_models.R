################################################################################
# 04_group_models.R
#
# Do the strain groupings explain biofilm production?
#
# One linear mixed model per candidate grouping, fit with lme4, compared within
# each set of strains two ways:
#
#   * a Kenward-Roger F-test against the global-mean model -- the P value
#     reported. KR rather than a chi-square likelihood-ratio test because there
#     are only 8-11 strains, and the chi-square LRT is anticonservative with so
#     few groups.
#   * AICc and Akaike weights -- the multimodel comparison, summing to one
#     within each set of strains.
#
# plus pairwise contrasts between groups, as ratios of group means with 95% CIs
# and KR-adjusted P values (emmeans).
#
# The strain-level posteriors plotted in the figures come from brms
# (R/01_bayes_biofilm.R). This script is the frequentist half.
#
# Input : data/biofil.csv
# Output: output/model_comparison.csv  one row per grouping
#         output/model_contrasts.csv   pairwise group ratios
#
# Model, for well i of strain j:
#
#     log(OD550 corrected)[i,j] ~ 0 + group[j] + (1 | strain)
#
# No weights: the wells of an assay plate are equally precise. Cell-means
# coding, so each fixed effect is one group's mean on the log scale. AICc needs
# maximum-likelihood fits and KR needs REML fits, so each model is fit both ways.
#
# Candidate groupings (numbering follows the project notes):
#   1   global mean
#   2   ancestor vs evolved
#   3   mutation vs no mutation
#   4   spore vs total
#   5   sinR vs ywcC (BSU_38220)
#   5a  sinR vs ywcC vs no mutation
#
# Model 3 was previously skipped on the grounds that every evolved isolate
# carries a mutation. It does not: m23 and m26 carry none. On all eleven
# strains model 3 is {ancestor, m23, m26} against the eight sinR/ywcC clones.
# Among the ten evolved clones it coincides with model 4, because the two
# clones without a mutation are exactly the two spore-fraction clones.
#
# THE REFERENCE WELL. Models 2 and 3 on all eleven strains include the well
# recorded as "Ancestor", which the plate record names B. subtilis 168 delta 6
# -- probably not this experiment's ancestor (see warn_ancestor() in
# 00_setup.R). Both models inherit that. Models 4, 5 and 5a use evolved clones
# only and do not.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the biofilm project root (or R/).")
source(.setup[1])

for (p in c("lme4", "pbkrtest", "MuMIn", "emmeans")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. install.packages(\"%s\")", p, p))
  }
}

## ---- data ------------------------------------------------------------------

assay <- read.csv(file.path(DATA_DIR, "biofil.csv"), stringsAsFactors = FALSE)
assay <- assay[assay$clones %in% MODEL_ORDER, ]
meta  <- strain_meta(MODEL_ORDER)
meta$has_mutation <- ifelse(meta$mutated, "mutation", "none")

PARAMS <- list(biofilm = "OD550_Corrected")

## ---- which groupings, on which strains -------------------------------------

SCOPES <- list(
  list(
    id = "all (11)",
    strains = MODEL_ORDER,
    models = list(
      list(id = "1", label = "global mean",           grouping = NULL),
      list(id = "2", label = "ancestor vs evolved",   grouping = "origin"),
      list(id = "3", label = "mutation vs none",      grouping = "has_mutation")
    )
  ),
  list(
    id = "evolved (10)",
    strains = setdiff(MODEL_ORDER, "ancestor"),
    # model 3 is the same partition as model 4 here and is not repeated
    models = list(
      list(id = "1",  label = "global mean",          grouping = NULL),
      list(id = "4",  label = "spore vs total",       grouping = "cell"),
      list(id = "5a", label = "sinR vs ywcC vs none", grouping = "mutation_full")
    )
  ),
  list(
    id = "sinR/ywcC (8)",
    strains = meta$clone[meta$mutation_full %in% c("sinR", "ywcC")],
    models = list(
      list(id = "1", label = "global mean",           grouping = NULL),
      list(id = "5", label = "sinR vs ywcC",          grouping = "mutation_full")
    )
  )
)

## ---- fitting ---------------------------------------------------------------

fit_both <- function(dd, grouping) {
  form <- if (is.null(grouping)) y ~ 1 + (1 | strain) else y ~ 0 + grp + (1 | strain)
  fit <- function(reml) suppressMessages(suppressWarnings(
    lme4::lmer(form, data = dd, REML = reml)))
  list(reml = fit(TRUE), ml = fit(FALSE))
}

sigma_strain <- function(m) sqrt(as.numeric(lme4::VarCorr(m)$strain))

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

  for (scope in SCOPES) {
    d <- assay[assay$clones %in% scope$strains, ]
    base <- data.frame(y = log(d[[column]]),
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
        uses_reference_well = scope$id == "all (11)" && !is.null(mod$grouping),
        stringsAsFactors = FALSE)

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
print(show[, c("scope", "model", "grouping", "n_groups", "sigma_strain",
               "akaike_weight", "ddf", "p_kr", "singular")], row.names = FALSE)

cat("\nNotes\n")
cat("  * P is a Kenward-Roger F-test against the global mean within the same\n",
    "   strains. The chi-square LRT is kept in model_comparison.csv.\n", sep = "")
cat("  * Models 2 and 3 on all eleven strains include the reference well, which\n",
    "   is provisional (168 delta 6). Models 4, 5 and 5a do not use it.\n", sep = "")
cat("  * Among the ten evolved clones model 3 is identical to model 4.\n", sep = "")
if (any(comparison$singular)) {
  cat("  * Singular fits: ",
      paste(unique(with(comparison[comparison$singular, ],
                        paste(scope, "model", model))), collapse = "; "), "\n", sep = "")
}

warn_ancestor()
message("\nWrote output/model_comparison.csv and model_contrasts.csv")
