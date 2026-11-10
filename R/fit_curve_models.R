####
# Alternative dose-response curve models
####
#
# Fit functions for apply_fit(), one per curve shape, wrapping the models drc
# already ships. They exist because the single-agent side of the fitting layer
# otherwise offers only fit_drug_response_metrics — the 4-parameter
# log-logistic the pipeline itself runs — so "fit this with a different model"
# had no answer but writing one.
#
# Each returns its own parameters rather than the Metrics schema, which is the
# pattern musyc_fit_fn set and what apply_fit() supports for a custom assay:
# forcing a Gompertz into columns named for a log-logistic would fill them with
# values that do not mean what the column name says. Parameters that play the
# same role across models reuse the vocabulary logisticFit() already uses —
# h, x_inf, x_0, ec50 — and model-specific ones are named for themselves.
#
# Every function also returns r2, converged, n_conc and fit_type, so two models
# on the same data can be compared rather than only read.
#
# The trade: an assay written by these is not the Metrics schema, so nothing
# downstream plots it without work.

#' @keywords internal
.CURVE_MODELS <- list(
  ll3 = list(
    fct = "LL.3",
    names = c("h", "x_0", "ec50"),
    fit_type = "DRC3pLogLogisticLowerFixed",
    note = "3-parameter log-logistic, lower asymptote fixed at 0"
  ),
  ll5 = list(
    fct = "LL.5",
    names = c("h", "x_inf", "x_0", "ec50", "asymmetry"),
    fit_type = "DRC5pLogLogistic",
    note = "5-parameter log-logistic, asymmetric"
  ),
  weibull1 = list(
    fct = "W1.4",
    names = c("h", "x_inf", "x_0", "ec50"),
    fit_type = "DRC4pWeibullType1",
    note = "4-parameter Weibull type 1, asymmetric"
  ),
  weibull2 = list(
    fct = "W2.4",
    names = c("h", "x_inf", "x_0", "ec50"),
    fit_type = "DRC4pWeibullType2",
    note = "4-parameter Weibull type 2, asymmetric"
  ),
  gompertz = list(
    fct = "G.4",
    names = c("h", "x_inf", "x_0", "ec50"),
    fit_type = "DRC4pGompertz",
    note = "4-parameter Gompertz, asymmetric"
  ),
  hormesis = list(
    # BC.4 fixes the lower asymptote, so it takes four names and has no
    # x_inf. BC.5 keeps both asymptotes but needs more of the curve.
    fct = "BC.4",
    names = c("h", "x_0", "ec50", "hormesis"),
    fit_type = "DRC4pBrainCousens",
    note = "Brain-Cousens, for a low-dose stimulation a log-logistic cannot represent"
  )
)


# One fit, shared by every wrapper below.
#
# Returns a one-row list whatever happens: a model that does not converge is a
# result, not an error. apply_fit()'s on_error = "warn" would otherwise drop
# the cell, and a dropped cell is indistinguishable from one that was never
# asked for.
#' @keywords internal
.fit_curve_model <- function(avg_dt, model, x_col = "x",
                             conc_col = "Concentration") {
  spec <- .CURVE_MODELS[[model]]
  empty <- stats::setNames(
    as.list(rep(NA_real_, length(spec$names))), spec$names
  )
  failed <- c(empty, list(r2 = NA_real_, converged = FALSE,
                          ec50_extrapolated = NA, n_conc = 0L,
                          fit_type = spec$fit_type))

  if (!all(c(x_col, conc_col) %in% names(avg_dt))) {
    return(failed)
  }
  df_ <- data.frame(
    conc = as.numeric(avg_dt[[conc_col]]),
    resp = as.numeric(avg_dt[[x_col]])
  )
  df_ <- df_[is.finite(df_$conc) & is.finite(df_$resp) & df_$conc > 0, ]
  n_conc <- length(unique(df_$conc))
  failed$n_conc <- n_conc

  # Below the number of free parameters the fit is not underdetermined so much
  # as meaningless, and drc will happily return something for it.
  if (n_conc < length(spec$names)) {
    return(failed)
  }

  # drm() writes this argument into the global na.action option exactly as
  # written and never restores it, so an unqualified na.omit is deliberate and
  # the option is put back by hand. See the same note in fit_curves.R.
  old_na_action <- getOption("na.action")
  on.exit(options(na.action = old_na_action), add = TRUE)

  # By name, not by value. drc's model constructors record themselves with
  # match.call()[[1]], which is a closure rather than a symbol when the
  # function object is passed, and llogistic() then dies in as.character().
  fct <- do.call(spec$fct, list(names = spec$names),
                 envir = asNamespace("drc"))
  controls <- drc::drmc(relTol = 1e-04, errorm = FALSE, noMessage = TRUE,
                        rmNA = TRUE)

  # ec50 is a concentration, and every model here evaluates log(dose/ec50), so
  # a negative one makes the model undefined at every dose. Left free, drc's
  # self-starter reached exactly that on 11 of 24 real curves. Bounding it
  # below recovers all 11 and costs nothing on the curves that already fitted:
  # across those it never lowered r2, and raised it on three.
  lowerl <- rep(-Inf, length(spec$names))
  lowerl[match("ec50", spec$names)] <- 1e-12

  # Errors are failures; warnings are not. drc warns routinely on fits that
  # converge perfectly well, and treating a warning as a failure threw away
  # W1.4 on a clean sigmoid.
  fit <- withCallingHandlers(
    tryCatch(
      drc::drm(resp ~ conc, data = df_, logDose = NULL, fct = fct,
               control = controls, lowerl = lowerl, na.action = na.omit),
      error = function(e) NULL
    ),
    warning = function(w) invokeRestart("muffleWarning")
  )
  if (is.null(fit) || is.null(stats::coef(fit))) {
    return(failed)
  }

  co <- stats::coef(fit)
  # drc labels coefficients "<param>:(Intercept)"; keep our names.
  names(co) <- sub(":.*$", "", names(co))
  out <- as.list(co[spec$names])
  names(out) <- spec$names

  # A model can finish on parameters that make its own prediction undefined:
  # a log-logistic landing on a negative ec50 evaluates log(dose/ec50) and
  # returns NaN at every dose. coef() still returns numbers, so the residuals
  # are the only place it shows, and summing them with na.rm = TRUE scores an
  # all-NaN fit as a perfect one. Seen on real RV data, where it reported
  # r2 = 1 alongside ec50 = -1.2.
  #
  # The warning is still muffled rather than treated as the signal: drc warns
  # "NaNs produced" during optimisation on fits whose final residuals are all
  # finite.
  resid <- withCallingHandlers(
    stats::residuals(fit),
    warning = function(w) invokeRestart("muffleWarning")
  )
  if (length(resid) == 0L || !all(is.finite(resid)) ||
      !all(vapply(out, is.finite, logical(1)))) {
    return(failed)
  }

  # The bound keeps ec50 positive; it cannot keep it meaningful. On a response
  # too flat to have a midpoint the optimiser runs away instead, reaching 1e19
  # on a 7-point curve here while the residuals stayed finite. A thousandfold
  # beyond the tested range is far wider than any honest extrapolation — the
  # genuine ones land within about ten times the top dose — so past it the
  # parameter was not estimated from these data.
  conc_rng <- range(df_$conc)
  if (out$ec50 < conc_rng[1] / 1e3 || out$ec50 > conc_rng[2] * 1e3) {
    return(failed)
  }

  ss_res <- sum(resid^2)
  ss_tot <- sum((df_$resp - mean(df_$resp, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_

  # Between "estimated" and "runaway" sits a real and common case: a drug too
  # weak to reach its own midpoint over the doses tested. The fit describes
  # the data well — on this panel, r2 up to 0.98 with ec50 nine times the top
  # dose — but the ec50 was extrapolated, not measured, and read as a potency
  # it would rank that drug against others on a number no plate supports.
  # Refusing it would throw away a good fit; reporting it silently would not.
  extrapolated <- out$ec50 < conc_rng[1] || out$ec50 > conc_rng[2]

  c(out, list(r2 = r2, converged = TRUE, ec50_extrapolated = extrapolated,
              n_conc = n_conc, fit_type = spec$fit_type))
}


#' Alternative dose-response curve models
#'
#' Fit functions for \code{\link{apply_fit}}, each fitting one curve shape from
#' \code{drc} to the averaged response of a single-agent experiment.
#'
#' @details
#' These complement \code{\link{fit_drug_response_metrics}}, which fits the
#' 4-parameter log-logistic the processing pipeline uses. Reach for one of
#' these when that shape is the wrong one: \code{ll3_fit_fn} when a curve is
#' too sparse or too shallow to support four free parameters,
#' \code{ll5_fit_fn}, \code{weibull1_fit_fn}, \code{weibull2_fit_fn} and
#' \code{gompertz_fit_fn} when the response is asymmetric, and
#' \code{hormesis_fit_fn} for a low-dose stimulation that a log-logistic
#' cannot represent at all.
#'
#' Each returns the model's own parameters rather than the \code{Metrics}
#' schema. Parameters playing the same role across models share the names
#' \code{logisticFit()} uses — \code{h} for the slope, \code{x_0} and
#' \code{x_inf} for the asymptotes, \code{ec50} for the inflection point — and
#' the rest are named for what they are. Write them to an assay of their own:
#' an assay of these is not the \code{Metrics} schema, so nothing reading
#' \code{Metrics} will understand it.
#'
#' All of them also return \code{r2}, \code{converged}, \code{n_conc} and
#' \code{fit_type}, so two models fitted to the same data can be compared, and
#' \code{ec50_extrapolated}, which is \code{TRUE} when the fitted \code{ec50}
#' lies outside the concentrations tested. That happens whenever a drug is too
#' weak to reach its own midpoint on the plate: the fit can still be a good
#' one, but the \code{ec50} is not a measurement and should not be read as a
#' potency.
#'
#' A fit that does not converge returns a row of \code{NA} parameters with
#' \code{converged = FALSE}, rather than raising. A fit that failed is a
#' result; dropping the cell would make it indistinguishable from one that was
#' never fitted. \code{converged = FALSE} also covers a fit that finished on
#' parameters that say nothing about the data: one whose own model cannot
#' evaluate them, and one whose \code{ec50} ran a thousandfold outside the
#' tested concentrations, which happens on a response too flat to have a
#' midpoint.
#'
#' \code{ec50} is bounded below at a positive value, because every model here
#' evaluates \code{log(dose/ec50)} and is undefined without it.
#'
#' \code{ll3_fit_fn} and \code{hormesis_fit_fn} fix the lower asymptote at
#' zero, so neither suits a GR-normalized response that goes below it: the
#' model cannot reach the data and the misfit shows up as a poor \code{r2}
#' rather than as a failure.
#'
#' @param avg_dt \code{data.table} for one drug, cell line and normalization
#'   type, as \code{\link{apply_fit}} supplies it.
#' @param x_col string; the response column.
#' @param conc_col string; the concentration column.
#'
#' @return Named list of one row: the model's parameters, then \code{r2},
#'   \code{converged}, \code{ec50_extrapolated}, \code{n_conc} and
#'   \code{fit_type}.
#'
#' @examples
#' dt <- data.frame(Concentration = 10^(-3:3), x = c(1, .98, .9, .5, .2, .1, .08))
#' ll3_fit_fn(dt)
#' gompertz_fit_fn(dt)
#'
#' @seealso \code{\link{apply_fit}}, \code{\link{fit_drug_response_metrics}}
#' @name curve_models
#' @keywords metrics
NULL


#' @rdname curve_models
#' @export
ll3_fit_fn <- function(avg_dt, x_col = "x", conc_col = "Concentration") {
  .fit_curve_model(avg_dt, "ll3", x_col, conc_col)
}

#' @rdname curve_models
#' @export
ll5_fit_fn <- function(avg_dt, x_col = "x", conc_col = "Concentration") {
  .fit_curve_model(avg_dt, "ll5", x_col, conc_col)
}

#' @rdname curve_models
#' @export
weibull1_fit_fn <- function(avg_dt, x_col = "x", conc_col = "Concentration") {
  .fit_curve_model(avg_dt, "weibull1", x_col, conc_col)
}

#' @rdname curve_models
#' @export
weibull2_fit_fn <- function(avg_dt, x_col = "x", conc_col = "Concentration") {
  .fit_curve_model(avg_dt, "weibull2", x_col, conc_col)
}

#' @rdname curve_models
#' @export
gompertz_fit_fn <- function(avg_dt, x_col = "x", conc_col = "Concentration") {
  .fit_curve_model(avg_dt, "gompertz", x_col, conc_col)
}

#' @rdname curve_models
#' @export
hormesis_fit_fn <- function(avg_dt, x_col = "x", conc_col = "Concentration") {
  .fit_curve_model(avg_dt, "hormesis", x_col, conc_col)
}


# Descriptions attached once, so the provenance and the registry report a model
# name rather than "a function".
#
# The attribute is set directly rather than through set_fit_fn_info(). R sources
# a package in collation order, and this file sorts before fit_provenance.R
# where that helper is defined, so calling it here fails the package load with
# "could not find function". gDRcore's musyc_fit_fn.R can call it only because
# gDRutils is a different, already-loaded package by then.
for (.nm in names(.CURVE_MODELS)) {
  .fn_name <- paste0(.nm, "_fit_fn")
  .fn <- get(.fn_name)
  attr(.fn, "fit_fn_info") <- list(
    model = .CURVE_MODELS[[.nm]]$fct,
    description = .CURVE_MODELS[[.nm]]$note,
    reference = paste("Ritz et al. 2015, PLoS ONE 10(12):e0146021,",
                      "doi:10.1371/journal.pone.0146021")
  )
  assign(.fn_name, .fn)
}
rm(.nm, .fn_name, .fn)
