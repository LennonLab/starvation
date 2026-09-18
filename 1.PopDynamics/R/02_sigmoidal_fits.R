################################################################################
# 02_sigmoidal_fits.R
#
# Fit each abundance series with a sigmoidal curve on log-log axes.
#
# Input : data/pop_dynamics.csv
# Output: output/sigmoidal_parms.csv   fitted parameters, pooled and per replicate
#         output/sigmoidal_curves.csv  fitted curves for plotting
#
#   log10(N) = b + (a - b) / (1 + exp((m - log10(t)) / w))
#
# a is the level the population settles at, b the level it starts from, m the
# midpoint in log10 days, w the width of the transition, and z the residual SD.
# Fit by maximum likelihood with bbmle::mle2, assuming normal error on log10(N).
#
# Day-zero counts are kept: log10(0) is -Inf, which sends the curve to exactly
# b, so those points inform the starting asymptote and nothing else. They drop
# out of the figures, which are on a log time axis.
################################################################################

## Locate 00_setup.R whether you are in the project root or in R/.
.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the populationDynamics project root (or R/).")
source(.setup[1])

if (!requireNamespace("bbmle", quietly = TRUE)) {
  stop("Package 'bbmle' is required. install.packages(\"bbmle\")")
}

pop <- read.csv(file.path(DATA_DIR, "pop_dynamics.csv"))

# The three series: what to fit, and where sensible starting values sit. The
# starting values are the ones the original script used.
SERIES <- list(
  total = list(column = "total",     label = "Total (S + V)",
               start = list(a = 5.5, b = 8.5, m = 0.7, w = 0.16, z = 0.2)),
  spore = list(column = "spore_est", label = "Spore (S)",
               start = list(a = 6.0, b = 6.0, m = 0.7, w = 0.16, z = 0.2)),
  veg   = list(column = "veg_est",   label = "Vegetative (V)",
               start = list(a = 5.5, b = 7.9, m = 0.7, w = 0.16, z = 0.2))
)

# Bounds on the parameters, wide enough not to bind on a sensible fit but
# tight enough to keep the optimiser in the physical region. Unconstrained,
# it wanders into negative sd (hundreds of "NaNs produced" warnings) and, for
# at least one replicate, off to a settle level of 10^35 CFU/mL.
PARM_LOWER <- c(a = 3,  b = 3,  m = -1.0, w = 0.01, z = 1e-3)
PARM_UPPER <- c(a = 11, b = 11, m =  3.5, w = 2.00, z = 2)

#' Fit one series to one set of rows.
#'
#' @return the mle2 object, or NULL if the optimiser could not converge
fit_series <- function(d, column, start) {
  df <- data.frame(y = log10(d[[column]]), log_t = log10(d$time))
  df <- df[is.finite(df$y), ]          # drop non-positive abundances

  tryCatch(
    bbmle::mle2(y ~ dnorm(mean = sigmoid_log10(log_t, a, b, m, w), sd = z),
                start = start, data = df,
                method = "L-BFGS-B",
                lower = PARM_LOWER[names(start)],
                upper = PARM_UPPER[names(start)]),
    error = function(e) {
      warning("Fit failed for ", column, ": ", conditionMessage(e), call. = FALSE)
      NULL
    })
}

## ---- fit pooled across replicates, then each replicate on its own ----------

targets <- c(list(list(id = "pooled", rows = seq_len(nrow(pop)))),
             lapply(REPS, function(r) list(id = r, rows = which(pop$rep == r))))

parms  <- list()
curves <- list()
fits   <- list()

for (nm in names(SERIES)) {
  cfg <- SERIES[[nm]]
  for (tg in targets) {
    f <- fit_series(pop[tg$rows, ], cfg$column, cfg$start)
    if (is.null(f)) next
    fits[[paste(nm, tg$id, sep = ".")]] <- f

    est <- bbmle::coef(f)
    se  <- sqrt(diag(bbmle::vcov(f)))
    parms[[length(parms) + 1]] <- data.frame(
      series = nm, label = cfg$label, fit = tg$id,
      parameter = names(est), estimate = unname(est), se = unname(se[names(est)])
    )

    if (tg$id == "pooled") {
      log_t <- seq(log10(0.5), log10(max(pop$time)), length.out = 400)
      curves[[length(curves) + 1]] <- data.frame(
        series = nm, label = cfg$label,
        time = 10^log_t,
        fitted = 10^sigmoid_log10(log_t, est[["a"]], est[["b"]],
                                  est[["m"]], est[["w"]])
      )
    }
  }
}

parms  <- bind_rows(parms)
curves <- bind_rows(curves)

write.csv(parms,  file.path(OUT_DIR, "sigmoidal_parms.csv"),  row.names = FALSE)
write.csv(curves, file.path(OUT_DIR, "sigmoidal_curves.csv"), row.names = FALSE)
saveRDS(fits, file.path(OUT_DIR, "sigmoidal_fits.rds"))

## ---- report ----------------------------------------------------------------

cat("\nPooled fits (log10 CFU/mL, time in days)\n\n")
pooled <- parms %>%
  filter(fit == "pooled") %>%
  select(series, parameter, estimate, se) %>%
  mutate(across(c(estimate, se), ~ round(.x, 3)))
print(pooled, row.names = FALSE)

cat("\nStart and settle levels, back on the CFU/mL scale\n\n")
lv <- parms %>%
  filter(fit == "pooled", parameter %in% c("a", "b")) %>%
  select(series, parameter, estimate) %>%
  pivot_wider(names_from = parameter, values_from = estimate) %>%
  mutate(starts_at  = sprintf("%.2e", 10^b),
         settles_at = sprintf("%.2e", 10^a),
         change = ifelse(a >= b,
                         sprintf("%.0f-fold rise", 10^(a - b)),
                         sprintf("%.0f-fold fall", 10^(b - a)))) %>%
  select(series, starts_at, settles_at, change)
print(as.data.frame(lv), row.names = FALSE)

message("\nWrote output/sigmoidal_parms.csv and output/sigmoidal_curves.csv")
