####
# The pipeline's fit, plus the concentration at other response levels
####
#
# fit_drug_response_metrics() reports xc50 and nothing else on that axis, so
# IC90 has had no answer except inverting the curve by hand. This wraps the
# same fit and reads further points off the same curve: the model, its
# parameters and every existing column are unchanged.
#
# levels are inhibition percentages, which is the one thing to get right here.
# See the note on .xc_at_levels().

# Concentrations at given inhibition percentages, from an already-fitted row.
#
# levels are INHIBITION percentages, so level 90 asks for the concentration at
# which the response reaches 0.1. That is the convention CIScore_80 already
# uses on the combination side (gDRcore computes it as pick(0.2)); it is NOT
# the reading under which the number is the response level. xc50 satisfies both
# only because 50 is its own complement, which is exactly why this has to be
# written down rather than inferred.
#
# Capping goes through cap_xc50() so a level that extrapolates past the plate
# gets the same +/-Inf treatment xc50 already gets, rather than a second rule.
#' @keywords internal
.xc_at_levels <- function(metrics, levels, capping_fold = 5) {
  nms <- paste0("xc", levels)

  # xc50 is already in the fit's own output and means the same thing. Emitting
  # it again would put two elements of the same name in one row, so the fit's
  # value wins and the request is satisfied rather than duplicated.
  keep <- !nms %in% names(metrics)
  nms <- nms[keep]
  levels <- levels[keep]
  if (length(levels) == 0L) {
    return(list())
  }

  # Read the parameters positively rather than trusting every return path of
  # the fit to carry them: a missing element would be NULL, rep(NULL, n) is
  # NULL, and predict_conc_from_efficacy()'s length assertion would fail on a
  # row that merely represents a failed fit.
  pars <- lapply(c("x_inf", "x_0", "ec50", "h"), function(k) metrics[[k]])
  if (!all(vapply(pars, function(v) length(v) == 1L, logical(1)))) {
    return(stats::setNames(as.list(rep(NA_real_, length(nms))), nms))
  }

  n <- length(levels)
  conc <- predict_conc_from_efficacy(
    efficacy = 1 - levels / 100,
    x_inf = rep(pars[[1L]], n),
    x_0   = rep(pars[[2L]], n),
    ec50  = rep(pars[[3L]], n),
    h     = rep(pars[[4L]], n)
  )

  max_log <- metrics[["maxlog10Concentration"]]
  max_conc <- if (length(max_log) == 1L) 10^max_log else NA_real_
  out <- vapply(conc, function(v) {
    # cap_xc50() asserts a single finite-or-infinite number, so NA from a
    # failed fit and Inf from an unreachable level both bypass it unchanged.
    if (!is.finite(v) || !is.finite(max_conc)) {
      return(as.numeric(v))
    }
    cap_xc50(xc50 = v, max_conc = max_conc, capping_fold = capping_fold)
  }, numeric(1))

  stats::setNames(as.list(out), nms)
}


#' The pipeline's fit, with the concentration at further response levels
#'
#' Runs \code{\link{fit_drug_response_metrics}} and appends the concentration
#' at each requested inhibition level — \code{xc90} for IC90, \code{xc75} for
#' IC75, and so on.
#'
#' @details
#' The curve is not changed. This calls the pipeline's own fit with the
#' arguments it is given, then inverts the resulting model with
#' \code{\link{predict_conc_from_efficacy}} at the requested levels. Use it
#' when the fit is right and only \code{xc50} is too few points to read off it;
#' use one of the alternative curve models when the *shape* is wrong.
#'
#' \strong{\code{levels} are inhibition percentages.} \code{levels = 90} asks
#' for the concentration at which the response reaches \code{0.1}, not
#' \code{0.9}. This matches \code{CIScore_80} on the combination side.
#' \code{xc50} is consistent with either reading only because 50 is its own
#' complement, so the distinction matters at every other level.
#'
#' \code{xc50} is never recomputed: it is already produced by the fit and means
#' the same thing, so asking for level 50 returns the fit's own value rather
#' than a second column of the same name.
#'
#' The \code{xc} prefix is deliberate and stays normalization-agnostic, so the
#' same column name is honest on the GR scale, where "IC" is not.
#'
#' Levels that the curve never reaches come back as \code{Inf}, and levels that
#' extrapolate more than \code{capping_fold} beyond the tested concentrations
#' are capped to \code{±Inf} by \code{\link{cap_xc50}} — the same treatment
#' \code{xc50} already receives, not a parallel rule.
#'
#' The fit-control arguments are passed through so a dataset processed with
#' non-default settings can be matched exactly. \code{fit_SE()} records those
#' settings in the provenance recipe, so they can be read back with
#' \code{\link{get_fit_provenance}} and supplied here.
#'
#' @section Not registered for averaging:
#' The new columns are concentrations, but they are not in
#' \code{get_header("metric_average_fields")$geometric_mean}. Do not average or
#' merge them as metrics: \code{\link{average_biological_replicates_dt}}
#' derives its field list from every numeric column and would take an
#' arithmetic mean of a quantity spanning orders of magnitude. Nothing inside
#' gDR averages a custom assay, and \code{merge_SE()} fills with \code{NA}
#' rather than averaging, so this is a limitation to respect rather than a
#' trap already sprung.
#'
#' @param avg_dt \code{data.table} for one drug, cell line and normalization
#'   type, as \code{\link{apply_fit}} supplies it.
#' @param x_col string; the response column.
#' @param levels numeric vector of inhibition percentages between 1 and 99.
#' @param capping_fold,range_conc,pcutoff,n_point_cutoff,force_fit,cap passed
#'   to \code{\link{fit_drug_response_metrics}} unchanged.
#'
#' @return The named list \code{\link{fit_drug_response_metrics}} returns, plus
#'   one element per requested level named \code{xc<level>}.
#'
#' @examples
#' dt <- data.table::data.table(
#'   Concentration = 10^(-3:3),
#'   x = c(1, .98, .9, .5, .2, .1, .08),
#'   normalization_type = "RV"
#' )
#' fit_drug_response_metrics_xcx(dt, levels = c(75, 90))
#'
#' @seealso \code{\link{fit_drug_response_metrics}},
#'   \code{\link{predict_conc_from_efficacy}}, \code{\link{cap_xc50}}
#' @keywords fit_models
#' @export
fit_drug_response_metrics_xcx <- function(avg_dt,
                                          x_col = "x",
                                          levels = c(50, 90),
                                          capping_fold = 5,
                                          range_conc = c(5e-3, 5),
                                          pcutoff = 0.05,
                                          n_point_cutoff = 4L,
                                          force_fit = FALSE,
                                          cap = 0.1) {
  # Integer levels only: they become column names, and xc33.33333 is not one.
  # 0 and 100 are excluded because both are degenerate - 0% inhibition is the
  # zero dose and 100% is never reached by a curve with a non-zero asymptote.
  checkmate::assert_integerish(
    levels, lower = 1, upper = 99,
    any.missing = FALSE, unique = TRUE, min.len = 1L
  )

  metrics <- fit_drug_response_metrics(
    avg_dt,
    x_col = x_col,
    capping_fold = capping_fold,
    range_conc = range_conc,
    pcutoff = pcutoff,
    n_point_cutoff = n_point_cutoff,
    force_fit = force_fit,
    cap = cap
  )

  c(metrics, .xc_at_levels(metrics, levels, capping_fold = capping_fold))
}


# Described so the provenance and the catalogue name a model rather than "a
# function". set_fit_fn_info() can be called here, unlike in
# fit_curve_models.R: with no Collate field R sources alphabetically, and
# fit_xcx.R sorts after fit_provenance.R where the helper is defined.
fit_drug_response_metrics_xcx <- set_fit_fn_info(
  fit_drug_response_metrics_xcx,
  model = "3-parameter Hill, inverted at further response levels",
  description = paste(
    "The pipeline's own fit, plus the concentration at each requested",
    "inhibition level (xc90 = IC90, the concentration at which the response",
    "reaches 0.1). The curve is unchanged."
  ),
  reference = "Hafner et al. 2016, Nat Methods 13(6):521-527, doi:10.1038/nmeth.3853"
)
