####
# Hill model maths and metric setters
####
#
# The Hill fit itself is fit_drug_response_metrics() (fit_utils.R), with its optimiser in
# fit_solvers.R. This file holds the closed-form model used around it - predictions in
# both directions, the combination smoothing - and the setters other packages use to
# build Metrics rows.


#' Deprecated: fit curves with the old interface
#'
#' \code{fit_curves()} and \code{logisticFit()} were a second copy of the gDR Hill fit,
#' with hard-coded priors. They are now thin wrappers over
#' \code{\link{fit_drug_response_metrics}}, the single gDR Hill fit, whose priors and
#' bounds come from \code{inst/extdata/fit_profiles.json} and whose optimiser is a
#' parameter (\code{\link{solver_drc}} by default). They are deprecated and will be
#' defunct in the next release: call \code{fit_drug_response_metrics(avg_dt, x_0 = )}
#' once per normalization type instead.
#'
#' @param df_ \code{data.table} or \code{DFrame} with the concentrations in
#'   \code{series_identifiers}, the responses in \code{x}, optionally their standard
#'   deviation in \code{x_std}, and \code{normalization_type}.
#' @param series_identifiers string, the concentration column of \code{df_}.
#' @param e_0,GR_0 upper asymptote of the RV and of the GR fit.
#' @param n_point_cutoff,range_conc,force_fit,pcutoff,cap,capping_fold passed to
#'   \code{fit_drug_response_metrics()}.
#' @param normalization_type normalization types to fit, among \code{"GR"} and
#'   \code{"RV"}.
#' @param concs concentrations, not logged.
#' @param norm_values responses, the same length as \code{concs}.
#' @param std_norm_values their standard deviation.
#' @param x_0 upper asymptote, or \code{NA} to fit it.
#' @param priors,lower ignored: the start values and bounds are those of the
#'   normalization type's fit profile. A \code{lower} value below zero for
#'   \code{x_inf} still selects the GR profile, as the GR priors did.
#'
#' @return \code{fit_curves()}: a \code{data.table} with one row per normalization type,
#'   \code{fit_source} \code{"gDR"}, and rownames \code{"<type>_gDR"}.
#'   \code{logisticFit()}: a one-row \code{data.table} of metrics and fit parameters.
#'
#' @examples
#' suppressWarnings(logisticFit(
#'   c(0.001, 0.00316227766016838, 0.01, 0.0316227766016838),
#'   c(0.9999964000144, 0.999964001439942, 0.999640143942423, 0.996414342629482),
#'   rep(0.1, 4)
#' ))
#'
#' @keywords fit_curves
#' @name fit_curves-deprecated
NULL

#' @rdname fit_curves-deprecated
#' @export
fit_curves <- function(df_,
                       series_identifiers,
                       e_0 = 1,
                       GR_0 = 1,
                       n_point_cutoff = 4,
                       range_conc = c(5e-3, 5),
                       force_fit = FALSE,
                       pcutoff = 0.05,
                       cap = 0.1,
                       normalization_type = c("GR", "RV")) {
  .Deprecated("fit_drug_response_metrics", package = "gDRutils",
              msg = paste("fit_curves() is deprecated and will be defunct in the next release:",
                          "fit each normalization type with",
                          "fit_drug_response_metrics(avg_dt, x_0 = ...)."))
  checkmate::assert_string(series_identifiers)
  checkmate::assert_subset(normalization_type, c("GR", "RV"))
  stopifnot(inherits(df_, "data.table") || inherits(df_, "DFrame"))
  df_ <- data.table::as.data.table(df_)
  if (!all(c(series_identifiers, "x") %in% names(df_))) {
    stop(sprintf("missing one of the required fields: '%s'",
                 toString(c(series_identifiers, "x"))))
  }
  conc_id <- get_env_identifiers("concentration")
  has_std <- "x_std" %in% names(df_)
  # RV first, then GR, as the old copy returned them
  out <- lapply(intersect(c("RV", "GR"), normalization_type), function(nt) {
    sel <- df_$normalization_type == nt
    avg <- data.table::data.table(
      conc = df_[[series_identifiers]][sel],
      x = df_$x[sel],
      x_std = if (has_std) df_$x_std[sel] else NA_real_,
      normalization_type = nt
    )
    data.table::setnames(avg, "conc", conc_id)
    data.table::as.data.table(fit_drug_response_metrics(
      avg, x_0 = if (nt == "RV") e_0 else GR_0, n_point_cutoff = n_point_cutoff,
      range_conc = range_conc, force_fit = force_fit, pcutoff = pcutoff, cap = cap
    ))
  })
  out <- data.table::rbindlist(out, fill = TRUE)
  out$fit_source <- "gDR"
  rownames(out) <- paste0(out$normalization_type, "_", out$fit_source)
  out
}

#' @rdname fit_curves-deprecated
#' @export
logisticFit <- function(concs,
                        norm_values,
                        std_norm_values = NA,
                        x_0 = 1,
                        priors = NULL,
                        lower = NULL,
                        range_conc = c(5e-3, 5),
                        force_fit = FALSE,
                        pcutoff = 0.05,
                        cap = 0.1,
                        n_point_cutoff = 4,
                        capping_fold = 5) {
  .Deprecated("fit_drug_response_metrics", package = "gDRutils",
              msg = paste("logisticFit() is deprecated and will be defunct in the next release:",
                          "use fit_drug_response_metrics(avg_dt, x_0 = ...)."))
  if (length(concs) != length(norm_values)) {
    stop("unequal vector lengths for 'conc' and 'norm_values'")
  }
  # The old GR priors were the ones with x_inf allowed below zero.
  nt <- if (length(lower) >= 2L && isTRUE(lower[[2]] < 0)) "GR" else "RV"
  avg <- data.table::data.table(conc = concs, x = norm_values,
                                x_std = std_norm_values, normalization_type = nt)
  data.table::setnames(avg, "conc", get_env_identifiers("concentration"))
  data.table::as.data.table(fit_drug_response_metrics(
    avg, x_0 = x_0, capping_fold = capping_fold, range_conc = range_conc,
    pcutoff = pcutoff, n_point_cutoff = n_point_cutoff, force_fit = force_fit, cap = cap
  ))
}

#' Predict efficacy values given fit parameters and a concentration.
#'
#' Predict efficacy values given fit parameters and a concentration.
#'
#' @param c Numeric vector representing concentrations to predict efficacies for.
#' @param x_inf Numeric vector representing the asymptotic value of the sigmoidal fit to the dose-response
#'  data as concentration goes to infinity.
#' @param x_0 Numeric vector representing the asymptotic metric value corresponding to a concentration of 0
#'  for the primary drug.
#' @param ec50 Numeric vector representing the drug concentration at half-maximal effect.
#' @param h Numeric vector representing the hill coefficient of the fitted curve, which reflects how steep
#'  the dose-response curve is.
#' @keywords fit_curves
#'
#' @return Numeric vector representing predicted efficacies from given concentrations and fit parameters.
#'
#' @details The inverse of this function is \code{predict_conc_from_efficacy}.
#' @seealso predict_conc_from_efficacy
#'
#' @examples
#' predict_efficacy_from_conc(c = 1, x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
#'
#' @export
predict_efficacy_from_conc <- function(c, x_inf, x_0, ec50, h) {
  checkmate::assert_numeric(c)
  checkmate::assert_numeric(x_inf)
  checkmate::assert_numeric(x_0)
  checkmate::assert_numeric(ec50)
  checkmate::assert_numeric(h)
  assert_equal_input_len(outlier = c, x_inf, x_0, ec50, h)

  # avoid issues with c=0 for DRCConstantFitResult
  as.numeric(ifelse(c > 0,
                    x_inf + (x_0 - x_inf) / (1 + (c / ec50) ^ h),
                    x_0))
}

#' Predict a smoothed response for a drug combination
#'
#' Predicts a 'smooth' value for a single pair of drug concentrations by
#' snapping to the nearest available models from a metrics table and averaging
#' their predictions. This is the combination equivalent of
#' 'predict_efficacy_from_conc'.
#'
#' @param conc_1 A single numeric value for the desired concentration of the first drug.
#' @param conc_2 A single numeric value for the desired concentration of the second drug.
#' @param metrics_merged A data.table containing all pre-calculated curve fit parameters.
#'   Expects columns: 'dilution_drug', 'cotrt_value', 'ratio', 'ec50', 'h', 'x_inf', 'x_0'.
#'   The co-dilution fit used is the one whose 'ratio' equals
#'   \code{round_concentration(conc_2 / conc_1, ndigit = 1)} for the snapped concentrations.
#'
#' @return A single numeric value for the predicted 'smooth' response.
#'
#' @examples
#' mae <- get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se[1, 1], "Metrics")[normalization_type == "RV"]
#' predict_smooth_from_combo(conc_1 = 1.2, conc_2 = 9.8, metrics_merged = dt_metrics)
#'
#' @keywords fit_curves
#' @export
#'
predict_smooth_from_combo <- function(conc_1, conc_2, metrics_merged) {
  checkmate::assert_number(conc_1, lower = 0)
  checkmate::assert_number(conc_2, lower = 0)
  checkmate::assert_data_table(metrics_merged)
  checkmate::assert_names(
    colnames(metrics_merged),
    must.include = c("dilution_drug", "cotrt_value", "ratio", "ec50", "h", "x_inf", "x_0")
  )

  available_cotrt_1 <- unique(metrics_merged[dilution_drug == "drug_2"]$cotrt_value)
  available_cotrt_2 <- unique(metrics_merged[dilution_drug == "drug_1"]$cotrt_value)

  snapped_conc_1 <- .snap_conc_to_model(conc_1, available_cotrt_1)
  snapped_conc_2 <- .snap_conc_to_model(conc_2, available_cotrt_2)

  message(sprintf("Requested: (%.2f, %.2f) ==> Using models for nearest concentrations: (%.2f, %.2f)",
                  conc_1, conc_2, snapped_conc_1, snapped_conc_2))


  col_params <- metrics_merged[dilution_drug == "drug_1" & cotrt_value == snapped_conc_2, ]
  col_value <- if (NROW(col_params) == 1) {
    predict_efficacy_from_conc(conc_1,
                               col_params$x_inf,
                               col_params$x_0,
                               col_params$ec50,
                               col_params$h)
    } else {
      NA
    }

  row_params <- metrics_merged[dilution_drug == "drug_2" & cotrt_value == snapped_conc_1, ]
  row_value <- if (NROW(row_params) == 1) {
    predict_efficacy_from_conc(conc_2,
                               row_params$x_inf,
                               row_params$x_0,
                               row_params$ec50,
                               row_params$h)
    } else {
      NA
    }

  codil_value <- NA
  if (!is.na(snapped_conc_1) && !is.na(snapped_conc_2) && snapped_conc_1 > 0 && snapped_conc_2 > 0) {
    # the ratio as the co-dilution fits store it; compared outside [.data.table, where
    # `ratio` would be the column itself and the condition always TRUE
    target_ratio <- round_concentration(snapped_conc_2 / snapped_conc_1, ndigit = 1)
    is_codil <- metrics_merged$dilution_drug %in% "codilution" & !is.na(metrics_merged$ratio) &
      abs(metrics_merged$ratio - target_ratio) <= 1e-8 * target_ratio
    codil_params <- metrics_merged[is_codil]
    codil_value <- if (NROW(codil_params) == 1) {
      predict_efficacy_from_conc(conc_1 + conc_2,
                                 codil_params$x_inf,
                                 codil_params$x_0,
                                 codil_params$ec50,
                                 codil_params$h)
      } else {
        NA
      }
  }

  predicted_values <- c(col_value, row_value, codil_value)
  final_prediction <- mean(predicted_values, na.rm = TRUE)

  if (is.nan(final_prediction)) {
    NA_real_
  } else {
    final_prediction
  }
}

#' Predict a concentration for a given efficacy with fit parameters.
#'
#' Predict a concentration for a given efficacy with fit parameters.
#'
#' @details The inverse of this function is \code{predict_efficacy_from_conc}.
#'
#' @param efficacy Numeric vector representing efficacies to predict concentrations for.
#' @param x_inf Numeric vector representing the asymptotic value of the sigmoidal fit to the dose-response
#'  data as concentration goes to infinity.
#' @param x_0 Numeric vector representing the asymptotic metric value corresponding to a concentration of 0
#'  for the primary drug.
#' @param ec50 Numeric vector representing the drug concentration at half-maximal effect.
#' @param h Numeric vector representing the hill coefficient of the fitted curve, which reflects how steep
#' @keywords fit_curves
#'
#' @return Numeric vector representing predicted concentrations from given efficacies and fit parameters.
#'
#' @examples
#' predict_conc_from_efficacy(efficacy = c(1, 1.5), x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
#'
#' @seealso predict_efficacy_from_conc .calculate_x50
#' @export
predict_conc_from_efficacy <- function(efficacy, x_inf, x_0, ec50, h) {
  assert_equal_input_len(outlier = efficacy, x_inf, x_0, ec50, h)
  ifelse(efficacy > x_0,
    0,
    ifelse(efficacy < x_inf,
      Inf,
      ec50 * ((x_0 - x_inf) / (efficacy - x_inf) - 1) ^ (1 / h)
    )
  )
}

#' Snap a concentration to the nearest available model concentration
#'
#' Finds the value in a vector of available concentrations that is closest
#' (on a log scale) to a user-specified concentration. This is an internal
#' helper function.
#'
#' @param user_conc A single numeric value for the desired concentration.
#' @param available_concs A numeric vector of concentrations for which a model exists.
#'
#' @return A single numeric value from 'available_concs'.
#' @keywords internal
#'
.snap_conc_to_model <- function(user_conc, available_concs) {
  checkmate::assert_number(user_conc, lower = 0)
  checkmate::assert_numeric(available_concs)

  if (length(available_concs) == 0 || is.na(user_conc)) {
    return(NA_real_)
  }
  # find the index of the concentration that is closest on a log scale.
  closest_idx <- which.min(abs(log10(available_concs) - log10(user_conc)))
  available_concs[closest_idx]
}

#' Cap XC50 value.
#'
#' Set IC50/GR50 value to \code{Inf} or \code{-Inf} based on upper and lower limits.
#'
#' @details
#' Note: \code{xc50} and \code{max_conc} should share the same units.
#' Ideally, the \code{lower_cap} should be based on the lowest tested concentration.
#' However, since we don't record that, it is set 5 orders of magnitude below the highest dose.
#'
#' @param xc50 Numeric value of the IC50/GR50 to cap.
#' @param max_conc Numeric value of the highest concentration in a dose series used to calculate the \code{xc50}.
#' @param min_conc Numeric value of the lowest concentration in a dose series used to calculate the \code{xc50}.
#' If \code{NA} (default), using \code{max_conc/1e5} instead.
#' @param capping_fold Integer value of the fold number to use for capping. Defaults to \code{5}.
#' @keywords fit_curves
#'
#' @return Capped IC50/GR50 value.
#'
#' @examples
#' cap_xc50(xc50 = 1, max_conc = 2)
#' cap_xc50(xc50 = 2, max_conc = 5, min_conc = 1)
#' cap_xc50(xc50 = 26, max_conc = 5, capping_fold = 5)
#'
#' @export
cap_xc50 <- function(xc50, max_conc, min_conc = NA, capping_fold = 5) {
  checkmate::assert_numeric(capping_fold)
  checkmate::assert_number(xc50)
  checkmate::assert_number(max_conc)
  checkmate::assert_number(min_conc, na.ok = TRUE)

  upper_cap <- max_conc * capping_fold
  lower_cap <- if (!is.na(min_conc)) {
    min_conc / capping_fold
  } else {
    max_conc / (capping_fold * 1e5)
  }
  if (xc50 > upper_cap) {
    xc50 <- Inf
  } else if (xc50 < lower_cap) {
    xc50 <- -Inf
  }
  xc50
}

#' @keywords fit_curves
#' @export
.setup_metric_output <- function() {
  template <- .metric_output_template()
  lapply(template, function(x) NA)
}

#' Build a named list of NA values for all response metric columns.
#' @keywords internal
.metric_output_template <- function() {
  resp_metric_all_cols <- get_header("response_metrics")
  resp_metric_cols <- resp_metric_all_cols[!endsWith(resp_metric_all_cols, "_sd")]
  out <- as.list(rep(NA, length(resp_metric_cols)))
  names(out) <- resp_metric_cols
  out
}

# Sign of an undefined IC/GR50 from the response level: Inf above 0.5, -Inf at or below.
# Used by the metric setters below; the fit itself takes its threshold per normalization
# type from fit_config in inst/extdata/fit_profiles.json.
#' @keywords internal
.estimate_xc50 <- function(param) {
  if (all(is.na(param))) {
    NA
  } else if (all(param > 0.5, na.rm = TRUE)) {
    Inf
  } else if (all(param <= 0.5, na.rm = TRUE)) {
    -Inf
  } else {
    NA
  }
}

#' @importFrom checkmate assert_number
#' @keywords fit_curves
#' @export
.set_mean_params <- function(out, mean_norm_value) {
  assert_number(mean_norm_value)

  out$xc50 <- .estimate_xc50(mean_norm_value)

  out$x_0 <- out$x_inf <- out$x_mean <- mean_norm_value
  out$x_AOC_range <- out$x_AOC <- .calculate_complement(mean_norm_value)
  out
}

#' Set fit parameters for a constant fit.
#'
#' Replace values for flat fits: ec50 = 0, h = 0.0001 and xc50 = +/- Inf
#'
#' @param out Named list of fit parameters.
#' @param mean_norm_value Numeric value that be used to set all parameters
#' that can be calculated from the mean.
#' @keywords fit_curves
#' @return Modified named list of fit parameters.
#'
#' @examples
#' na <- list(x_0 = NA)
#' set_constant_fit_params(na, mean_norm_value = 0.6)
#'
#' @export
set_constant_fit_params <- function(out, mean_norm_value) {
  out$fit_type <- "DRCConstantFitResult"
  out$ec50 <- 0
  out$h <- 0.0001
  out$r2 <- 0

  out <- .set_mean_params(out, mean_norm_value)
  out
}

#' Set fit parameters for an invalid fit.
#' @param out Named list of fit parameters.
#' @param norm_values Numeric vector used to estimate an \code{xc50} value.
#' @keywords fit_curves
#' @return Modified named list of fit parameters.
#'
#' @examples
#' .set_invalid_fit_params(list(), norm_values = rep(0.3, 6))
#'
#' @export
.set_invalid_fit_params <- function(out, norm_values) {
  out$fit_type <- "DRCInvalidFitResult"
  out$r2 <- NA
  out$xc50 <- .estimate_xc50(norm_values)
  out
}

#' @export
#' @keywords fit_curves
.calculate_complement <- function(x) {
  1 - x
}
