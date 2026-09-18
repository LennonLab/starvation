################################################################################
# 04_group_models.R
#
# Do the strain groupings explain the growth parameters? Fits a hierarchical
# model per grouping and compares them by WAIC.
#
# Input : data/comp_data_annotated.csv
# Output: output/model_comparison.csv   -- one row per grouping: how much
#                                          between-strain variation it leaves,
#                                          a likelihood-ratio test, and WAIC
#         output/model_contrasts.csv    -- group means and pairwise ratios
#
# Model, for curve i of strain j:
#
#     y[i,j] = log(parameter)
#     y[i,j]   ~ Normal( theta[j], sigma2_e / w[i,j] )   curve level
#     theta[j] ~ Normal( x[j]'beta, sigma2_s )           strain level
#     beta     ~ Normal(0, 100)
#     sigma2_e, sigma2_s ~ InverseGamma(0.01, 0.01)
#
# w[i,j] is the inverse-variance weight of the curve fit, scaled to mean 1 so
# that sigma2_e stays interpretable. x[j] is a cell-means indicator for the
# grouping under test, so each element of beta is one group's mean on the log
# scale and exp(beta_a - beta_b) is a ratio of group means.
#
# Candidate groupings (the numbering follows the project notes):
#   1  global mean -- all strains share one mean
#   2  ancestor vs. evolved
#   3  mutation vs. no mutation
#   4  spore vs. total
#   5  sinR vs. ywcC (BSU_38220)
#
# On this strain set model 3 is model 2: every evolved isolate carries a
# mutation and the ancestor carries none, so the two groupings are the same
# partition. They would differ only if the unmutated spore isolates (S1, S6,
# S11, S22, S51, S95) were included; those are in
# data/comp_data_annotated.csv and are not used here.
#
# Models 4 and 5 apply to subsets of the strains, so they are compared against
# their own null within that subset. Nothing is comparable across scopes.
#
# What to read. The question is whether a grouping explains variation *between
# strains*, so the headline numbers are sigma_strain -- the between-strain SD
# left after the grouping -- and the group ratios with their credible
# intervals. A grouping that explains something shrinks sigma_strain relative
# to the global-mean model.
#
# WAIC is reported too, but with a caveat: it scores prediction of one more
# curve from a strain already in the model, and the strain effect has already
# absorbed the group difference by then, so it is close to blind to the
# question being asked. With six curves per strain it is also unreliable
# (loo flags p_waic > 0.4). Treat it as a footnote, not as the answer.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the growthCurves project root (or R/).")
source(.setup[1])

if (!requireNamespace("loo", quietly = TRUE)) {
  stop("Package 'loo' is required for WAIC. install.packages(\"loo\")")
}

set.seed(20230826)

N_CHAINS <- 4L
N_BURN   <- 2000L
N_ITER   <- 5000L    # draws kept per chain

## ---- data ------------------------------------------------------------------

fits <- read.csv(file.path(DATA_DIR, "comp_data_annotated.csv"),
                 stringsAsFactors = FALSE)
fits <- fits[fits$clones %in% MODEL_ORDER, ]

meta <- strain_meta(MODEL_ORDER)

## ---- sampler ---------------------------------------------------------------

#' Gibbs sampler for the hierarchical model above.
#'
#' Every full conditional is conjugate, so this needs no external sampler.
#'
#' @param y      log parameter, one value per curve
#' @param w      curve weight (inverse variance, scaled to mean 1)
#' @param strain factor identifying the strain of each curve
#' @param X      strain-level design matrix, one row per level of `strain`
fit_hier <- function(y, w, strain, X,
                     n_chains = N_CHAINS, n_iter = N_ITER, n_burn = N_BURN,
                     beta_prior_var = 100, a0 = 0.01, b0 = 0.01) {
  strain <- droplevels(as.factor(strain))
  J <- nlevels(strain)
  N <- length(y)
  p <- ncol(X)
  stopifnot(nrow(X) == J)

  # Per-strain sufficient statistics, recomputed only where they depend on
  # theta inside the loop.
  Sw <- tapply(w, strain, sum)
  Tw <- tapply(w * y, strain, sum)
  idx <- as.integer(strain)

  XtX_prior <- diag(1 / beta_prior_var, p)

  one_chain <- function() {
    theta    <- tapply(y, strain, mean)
    sigma2_e <- var(y)
    sigma2_s <- var(theta)
    if (!is.finite(sigma2_s) || sigma2_s <= 0) sigma2_s <- 0.01

    beta_out  <- matrix(NA_real_, n_iter, p, dimnames = list(NULL, colnames(X)))
    sig_out   <- matrix(NA_real_, n_iter, 2,
                        dimnames = list(NULL, c("sigma2_e", "sigma2_s")))
    ll_out    <- matrix(NA_real_, n_iter, N)
    beta      <- rep(0, p)

    for (it in seq_len(n_iter + n_burn)) {
      # theta | beta, variances
      mu_s  <- as.vector(X %*% beta)
      prec  <- Sw / sigma2_e + 1 / sigma2_s
      theta <- rnorm(J, (Tw / sigma2_e + mu_s / sigma2_s) / prec, 1 / sqrt(prec))

      # beta | theta, sigma2_s
      V    <- chol2inv(chol(crossprod(X) / sigma2_s + XtX_prior))
      m    <- V %*% (crossprod(X, theta) / sigma2_s)
      beta <- as.vector(m + t(chol(V)) %*% rnorm(p))

      # variances
      resid_e  <- y - theta[idx]
      sigma2_e <- 1 / rgamma(1, a0 + N / 2, b0 + sum(w * resid_e^2) / 2)
      resid_s  <- theta - as.vector(X %*% beta)
      sigma2_s <- 1 / rgamma(1, a0 + J / 2, b0 + sum(resid_s^2) / 2)

      if (it > n_burn) {
        k <- it - n_burn
        beta_out[k, ] <- beta
        sig_out[k, ]  <- c(sigma2_e, sigma2_s)
        ll_out[k, ]   <- dnorm(y, theta[idx], sqrt(sigma2_e / w), log = TRUE)
      }
    }
    list(beta = beta_out, sigma = sig_out, loglik = ll_out)
  }

  chains <- lapply(seq_len(n_chains), function(i) one_chain())
  list(
    beta   = do.call(rbind, lapply(chains, `[[`, "beta")),
    sigma  = do.call(rbind, lapply(chains, `[[`, "sigma")),
    loglik = do.call(rbind, lapply(chains, `[[`, "loglik")),
    strain_levels = levels(strain)
  )
}

## ---- what to compare -------------------------------------------------------

# Each scope is a set of strains; within it, each candidate grouping is fit and
# compared by WAIC. `grouping = NULL` means the global-mean model.
SCOPES <- list(
  list(
    id = "all (11)",
    strains = MODEL_ORDER,
    note = "ancestor + 10 evolved clones",
    models = list(
      list(id = "1", label = "global mean",          grouping = NULL),
      list(id = "2", label = "ancestor vs evolved",  grouping = "origin")
    )
  ),
  list(
    id = "evolved (10)",
    strains = setdiff(MODEL_ORDER, "ancestor"),
    note = "10 evolved clones",
    models = list(
      list(id = "1", label = "global mean",              grouping = NULL),
      list(id = "4", label = "spore vs total",           grouping = "cell"),
      list(id = "5a", label = "sinR vs ywcC vs spormut", grouping = "mutation_full")
    )
  ),
  list(
    id = "sinR/ywcC (8)",
    strains = meta$clone[meta$mutation_full %in% c("sinR", "ywcC")],
    note = "8 non-sporulation-mutant clones",
    models = list(
      list(id = "1", label = "global mean",   grouping = NULL),
      list(id = "5", label = "sinR vs ywcC",  grouping = "mutation_full")
    )
  )
)

PARAMS <- list(umax = "umax", A = "A", L = "L")

## ---- run -------------------------------------------------------------------

comparison <- list()
contrasts  <- list()

for (pname in names(PARAMS)) {
  column <- PARAMS[[pname]]
  se_col <- paste0(column, ".se")

  for (scope in SCOPES) {
    d <- fits[fits$clones %in% scope$strains, ]
    d$strain <- factor(d$clones, levels = scope$strains)

    y <- log(d[[column]])
    w <- 1 / d[[se_col]]^2
    w <- w / mean(w)                      # keeps sigma2_e interpretable

    m <- meta[match(scope$strains, meta$clone), ]

    waics <- list()
    mls   <- list()
    for (mod in scope$models) {
      # Cell-means coding: one column per group, so each element of beta is
      # that group's mean on the log scale.
      X <- if (is.null(mod$grouping)) {
        matrix(1, nrow(m), 1, dimnames = list(NULL, "(all)"))
      } else {
        g <- factor(m[[mod$grouping]])
        matrix(as.integer(outer(g, levels(g), "==")), nrow(m),
               dimnames = list(NULL, levels(g)))
      }

      f  <- fit_hier(y, w, d$strain, X)
      # loo warns that p_waic > 0.4 for some curves; that is the unreliability
      # noted in the header, not a failure. Report it once at the end.
      ww <- withCallingHandlers(loo::waic(f$loglik),
                                warning = function(cnd) invokeRestart("muffleWarning"))
      waics[[mod$id]] <- ww

      # Maximum-likelihood fit of the same model, for a likelihood-ratio test.
      # lme4's prior `weights` scale the residual variance exactly as
      # sigma2_e / w does in the Gibbs sampler.
      ml <- NULL
      if (requireNamespace("lme4", quietly = TRUE)) {
        dd <- data.frame(y = y, w = w, strain = d$strain,
                         grp = if (is.null(mod$grouping)) factor("all")
                               else factor(m[[mod$grouping]][match(d$clones, m$clone)]))
        form <- if (is.null(mod$grouping)) y ~ 1 + (1 | strain)
                else y ~ 0 + grp + (1 | strain)
        ml <- suppressMessages(suppressWarnings(
          lme4::lmer(form, data = dd, weights = w, REML = FALSE)))
      }
      mls[[mod$id]] <- ml

      comparison[[length(comparison) + 1]] <- data.frame(
        parameter = pname, scope = scope$id, model = mod$id,
        grouping = mod$label, n_groups = ncol(X),
        # between-strain SD remaining after the grouping (log scale)
        sigma_strain = median(sqrt(f$sigma[, "sigma2_s"])),
        sigma_curve  = median(sqrt(f$sigma[, "sigma2_e"])),
        aic = if (is.null(ml)) NA_real_ else AIC(ml),
        lrt_p = NA_real_,
        waic = ww$estimates["waic", "Estimate"],
        waic_se = ww$estimates["waic", "SE"],
        p_waic = ww$estimates["p_waic", "Estimate"]
      )

      # Group means and pairwise ratios, on the measured scale.
      if (ncol(X) > 1) {
        gm <- exp(f$beta)
        pairs <- utils::combn(colnames(X), 2, simplify = FALSE)
        for (pr in pairs) {
          ratio <- gm[, pr[1]] / gm[, pr[2]]
          q <- quantile(ratio, c(.025, .5, .975))
          contrasts[[length(contrasts) + 1]] <- data.frame(
            parameter = pname, scope = scope$id, model = mod$id,
            contrast = paste(pr[1], "/", pr[2]),
            ratio_2.5 = q[[1]], ratio_50 = q[[2]], ratio_97.5 = q[[3]],
            prob_lt_1 = mean(ratio < 1)
          )
        }
      }
    }

    # Expected log predictive density of each grouping relative to the
    # global-mean model. loo_compare() would report this against whichever
    # model is best rather than against the global one, so take the paired
    # pointwise difference directly -- the same quantity, fixed reference.
    base <- waics[["1"]]$pointwise[, "elpd_waic"]
    for (id in setdiff(names(waics), "1")) {
      dpt <- waics[[id]]$pointwise[, "elpd_waic"] - base
      i <- which(sapply(comparison, function(x)
        x$parameter == pname && x$scope == scope$id && x$model == id))
      comparison[[i]]$elpd_diff_vs_global <- sum(dpt)
      comparison[[i]]$se_diff             <- sqrt(length(dpt)) * sd(dpt)

      if (!is.null(mls[[id]]) && !is.null(mls[["1"]])) {
        a <- suppressMessages(anova(mls[["1"]], mls[[id]]))
        comparison[[i]]$lrt_p <- a$`Pr(>Chisq)`[2]
      }
    }
  }
}

comparison <- bind_rows(comparison)
contrasts  <- bind_rows(contrasts)
rownames(contrasts) <- NULL

comparison$elpd_diff_vs_global[comparison$model == "1"] <- 0
comparison$se_diff[comparison$model == "1"] <- 0

write.csv(comparison, file.path(OUT_DIR, "model_comparison.csv"), row.names = FALSE)
write.csv(contrasts,  file.path(OUT_DIR, "model_contrasts.csv"),  row.names = FALSE)

## ---- report ----------------------------------------------------------------

fmt_p <- function(x) ifelse(is.na(x), "-",
                            ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))

cat("\n=========== Does the grouping explain variation between strains? ==========\n")
cat("sigma_strain: between-strain SD (log scale) left after the grouping.\n")
cat("A grouping that explains something drives it below the global-mean row.\n")
cat("WAIC differences and their standard errors are in model_comparison.csv.\n")
for (pname in names(PARAMS)) {
  cat("\n--", pname, "--\n")
  x <- comparison[comparison$parameter == pname, ]
  print(data.frame(
    scope = x$scope, model = x$model, grouping = x$grouping,
    groups = x$n_groups,
    sigma_strain = round(x$sigma_strain, 3),
    LRT_p = fmt_p(x$lrt_p),
    WAIC = round(x$waic, 1)), row.names = FALSE)
}

cat("\n================ Group ratios (95% credible interval) ==================\n")
cat("Ratio of group means on the measured scale; prob_lt_1 is the posterior\n")
cat("probability that the first group is the smaller of the two.\n")
for (pname in names(PARAMS)) {
  cat("\n--", pname, "--\n")
  x <- contrasts[contrasts$parameter == pname, ]
  print(data.frame(
    scope = x$scope, contrast = x$contrast,
    ratio = sprintf("%.2f [%.2f, %.2f]", x$ratio_50, x$ratio_2.5, x$ratio_97.5),
    prob_lt_1 = round(x$prob_lt_1, 3)), row.names = FALSE)
}

cat("\nNotes\n")
cat("  * Model 3 (mutation vs. no mutation) is the same partition as model 2 on\n",
    "   these strains, so it is not fit separately. The two would differ only\n",
    "   with the unmutated spore isolates (S1, S6, S11, S22, S51, S95) in the\n",
    "   data; they are in comp_data_annotated.csv but out of scope here.\n", sep = "")
cat("  * In model 2 the ancestor group holds one strain, so its group mean and\n",
    "   that strain's own effect are the same quantity, and the between-strain\n",
    "   SD is informed only by the ten evolved clones. The ancestor's six\n",
    "   curves are technical replicates, not independent strains -- read that\n",
    "   contrast accordingly.\n", sep = "")
cat("  * WAIC scores prediction of one more curve from a strain already in the\n",
    "   model, by which point the strain effect has absorbed the group\n",
    "   difference, so it barely responds to the grouping. loo also flags it as\n",
    "   unreliable here (p_waic > 0.4 for some curves, with only six curves per\n",
    "   strain). It is in the table for completeness; sigma_strain, the LRT and\n",
    "   the ratios are what answer the question.\n", sep = "")

## ---- sampler cross-check ---------------------------------------------------
# The Gibbs sampler and lme4 fit the same model by different routes, so the
# group means should agree closely.
if (requireNamespace("lme4", quietly = TRUE)) {
  d  <- fits[fits$clones %in% MODEL_ORDER, ]
  dd <- data.frame(y = log(d$umax), strain = factor(d$clones),
                   grp = meta$origin[match(d$clones, meta$clone)],
                   w = (1 / d$umax.se^2) / mean(1 / d$umax.se^2))
  ml <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ 0 + grp + (1 | strain), data = dd, weights = w, REML = FALSE)))
  b  <- lme4::fixef(ml)
  bay <- contrasts$ratio_50[contrasts$parameter == "umax" &
                              contrasts$scope == "all (11)"][1]
  cat(sprintf("\nCross-check (umax, ancestor/evolved ratio): Gibbs %.3f, lme4 %.3f\n",
              bay, exp(b[["grpAncestor"]] - b[["grpEvolved"]])))
}

message("\nWrote output/model_comparison.csv and output/model_contrasts.csv")
