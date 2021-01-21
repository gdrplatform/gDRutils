# TODO: Delete me once no longer in use. 
#' @export
#'
RVGRfits <- function(df_,
		      e_0 = 1,
		      GR_0 = 1,
		      n_point_cutoff = 4,
		      range_conc = c(5e-3, 5),
		      force = FALSE, 
                      pcutoff = 0.05) {
  .Deprecated("fit_curves", package="gDRutils")
  fit_curves(df_ = df_, e_0 = e_0, GR_0 = GR_0, n_point_cutoff = n_point_cutoff, 
    range_conc = range_conc, force = force, pcutoff = pcutoff)
}


#' Fit curves
#'
#' Fit GR and RV curves from a data.frame.
#'
#' @param df_ data.frame containing data to fit. See details.
#' @param e_0 numeric value representing the \code{x_0} value for the RV curve.
#' Defaults to \code{1}.
#' @param GR_0 numeric value representing the \code{x_0} value for the GR curve.
#' Defaults to \code{1}.
#' @param n_point_cutoff integer of how many points should be considered the minimum required to try to fit a curve.
#' @param range_conc numeric vector of length 2 indicating the lower and upper concentration ranges. See details.
#' @param force boolean indicating whether or not to force a constant fit.
#' @param pcutoff numeric of pvalue significance threshold above which to use a constant fit.
#'
#' @return data.frame of RV and GR fit parameters.
#'
#' @details
#' The \code{df_} expects the following columns:
#' \itemize{
#'  \item{\code{Concentration}}{unlogged concentration values}, 
#'  \item{\code{RelativeViability}}{normalized values},
#'  \item\code{std_RelativeViability}{standardized normalized values}
#' }
#'
#' The \code{range_conc} is used to calculate the \code{x_AOC_range} statistic.
#' The purpose of this statistic is to enable comparison across different experiments with slightly 
#' different concentration ranges. 
#'
#' @export
#'
fit_curves <- function(df_,
		       e_0 = 1,
		       GR_0 = 1,
		       n_point_cutoff = 4,
		       range_conc = c(5e-3, 5),
		       force = FALSE, 
                       pcutoff = 0.05) {

  checkmate::assert_class(df_, "data.frame")

  med_conc <- stats::median(df_$Concentration)
  min_conc <- min(df_$Concentration)

  df_RV <- gDRutils::logisticFit(
    df_$Concentration,
    df_$RelativeViability,
    df_$std_RelativeViability,
    x_0 = e_0,
    priors = c(2, 0.4, 1, med_conc),
    lower = c(0.1, 0, 0, min_conc / 10),
    range_conc = range_conc,
    force = force,
    pcutoff = pcutoff, 
    n_point_cutoff = n_point_cutoff
  )

  df_GR <- gDRutils::logisticFit(
    df_$Concentration,
    df_$GRvalue,
    df_$std_GRvalue,
    x_0 = GR_0,
    priors = c(2, 0.1, 1, med_conc),
    lower = c(0.1, -1, -1, min_conc / 10),
    range_conc = range_conc,
    force = force,
    pcutoff = pcutoff, 
    n_point_cutoff = n_point_cutoff
  )

  df_metrics <- rbind(df_RV, df_GR)
  rownames(df_metrics) <- c("RV", "GR")

  return(df_metrics)
}


#' Logistic fit
#'
#' Fit a logistic curve to drug response data.
#'
#' @param concs concentrations that have not been transformed into log space.
#' @param norm_values normalized response values (Untreated = 1).
#' @param std_norm_values std of values.
#' @param x_0 upper limit.
#' Defaults to \code{1}. For co-treatments, this value should be set to \code{NA}.
#' @param priors numeric vector containing starting values for all.
#' mean parameters in the model. Overrules any self starter function.
#' @param lower numeric vector of lower limits for all parameters in a 4-param model.
#' @param range_conc range of concentration for calculating AOC_range.
#' @param force boolean indicating whether or not to force a parameter-based fit.
#' @param pcutoff numeric of pvalue significance threshold above which to use a constant fit.
#' @param cap numeric value capping \code{norm_values} to stay below (\code{x_0} + cap).
#' @param n_point_cutoff integer indicating number of unique concentrations required to fit curve.
#'
#' @return data.frame with metrics and fit parameters.
#'
#' @details
#' Implementation of the genedata approach for curve fit: 
#' https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
#'
#' Parameter names correspond to the following definitions:
#' \itemize{
#'  \item{x_mean}{The mean of a given dose-response metric}
#'  \item{x_AOC_range}{The range of the area over the curve}
#'  \item{x_AOC}{The area over the GR curve or, respectively, under the relative cell count curve, averaged over the range of concentration values}
#'  \item{xc50}{The concentration at which the effect reaches a value of 0.5 based on interpolation of the fitted curve}
#'  \item{x_max}{The maximum effect of the drug}
#'  \item{c50}{The drug concentration at half-maximal effect}
#'  \item{x_inf}{The asymptotic value of the sigmoidal fit to the dose-response data as concentration goes to infinity}
#'  \item{x_0}{The asymptotic metric value corresponding to a concentration of 0 for the primary drug}
#'  \item{h}{The hill coefficient of the fitted curve, which reflects how steep the dose-response curve is}
#'  \item{r2}{The goodness of the fit}
#'  \item{x_sd_avg}{The standard deviation of GR/IC}
#'  \item{fit_type}{This will be given by one of the following:
#'    \item{"DRC4pHillFitModel"}{Successfully fit with a 4-parameter model}
#'    \item{"DRC3pHillFitModelFixS0"}{Successfully fit with a 3-parameter model}
#'    \item{"DRCConstantFitResult"}{Successfully fit with a constant fit}
#'    \item{"DRCTooFewPointsToFit"}{Not enough points to run a fit}
#'    \item{"DRCInvalidFitResult"}{Fit was attempted but failed}
#'  }
#'  \item{maxlog10Concentration}{The highest log10 concentration}
#'  \item{N_conc}{Number of unique concentrations}
#' }
#'
#' @export
#'
logisticFit <-
  function(concs,
           norm_values,
           std_norm_values,
           x_0 = 1,
           priors = NULL,
           lower = NULL,
           range_conc = c(5e-3, 5),
           force = FALSE,
           pcutoff = 0.05,
           cap = 0.1,
           n_point_cutoff = 4) {

    if (length(concs) != length(norm_values)) {
      stop("unequal vector lengths for 'conc' and 'norm_values'")
    }

    resp_metric_cols <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
    out <- as.data.frame(matrix(NA, ncol = length(resp_metric_cols)))
    names(out) <- resp_metric_cols

    # Check that values have not been logged yet. 
    if (any(concs < 0)) {
      stop("function accepts only unlogged concentrations, negative concentrations are detected")
    }

    # Cap norm_values at (x_0 + cap) so as not to throw off the fit.
    norm_values <- pmin(norm_values, (ifelse(is.na(x_0), 1, x_0) + cap))

    ## Calculate metrics that do not require fitting.
    out$maxlog10Concentration <- log10(max(concs))
    out$N_conc <- length(unique(concs))
    out$x_sd_avg <- mean(std_norm_values, na.rm = TRUE)

    df_ <- data.frame(concs = concs,
                      norm_values = norm_values)

    freq <- table(df_$concs)
    xAvg <- df_
    if (any(dups <- freq != 1L)) {
      warning(sprintf("duplicate concentrations were found: '%s', averaging values across a single concentration", 
        paste0(names(freq)[dups], collapse = ", ")))
      xAvg <- stats::aggregate(
	list(norm_values = norm_values),
	by = list(conc = df_$concs),
	FUN = function(x) {mean(x, na.rm = TRUE)}
      )
    }

    mean_norm_value <- mean(xAvg$norm_values, na.rm = TRUE)
    # Set temp values if fit fails.
    out$x_mean <- mean_norm_value
    out$x_AOC <- 1 - mean_norm_value

    ## 'x_max' can be considered either the lowest readout (max efficacy) 
    ## or the efficacy at the max concentration. We take the min 
    ## of the two highest concentrations as a compromise.
    l <- nrow(xAvg) # NAs have been removed. 
    out$x_max <- min(xAvg$norm_values[c(l, l - 1)], na.rm = TRUE)

    # Replace values for flat fits: c50 = 0, h = 0.0001 and xc50 = +/- Inf
    .set_constant_fit_params <- function(out) {
      out$fit_type <- 'DRCConstantFitResult'
      out$c50 <- 0
      out$h <- 0.0001
      out$xc50 <- .estimate_xc50(mean_norm_value)
      out$x_0 <- out$x_inf <- out$x_mean <- mean_norm_value
      out$x_AOC_range <- out$x_AOC <- 1 - mean_norm_value
      if (!is.na(x_0)) {
        warning(sprintf("overriding original x_0 argument '%s' with '%s'", x_0, mean_norm_value))
      }
      out
    }

    non_na_avg_norm <- !is.na(xAvg$norm_values)

    ## Test for known failure points.
    if (length(unique(xAvg$norm_values[non_na_avg_norm])) == 1L) {
      out <- .set_constant_fit_params(out)
      return(out)
    }

    if (sum(non_na_avg_norm) < n_point_cutoff) {
      out$fit_type <- 'DRCTooFewPointsToFit'
      out$xc50 <- .estimate_xc50(norm_values)
      return(out)
    }

    # Set fit parameters and boundaries.
    fit_param <- c("h", "x_inf", "x_0", "c50")
    controls <- drc::drmc(relTol = 1e-06, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)

    ## Perform a 3-param or 4-param fit. 
    ## Fit type is determined based on number of free variables available.
    tryCatch({
      if (!is.na(x_0)) { # For co-treatments.
	# Override existing params for removal of x_0 parameter. 
	fit_param <- fit_param[-3]
	priors <- priors[-3]
	lower <- lower[-3]

	output_model_new <- try(drc::drm(
	  norm_values ~ concs,
	  data = df_,
	  logDose = NULL,
	  fct = drc::LL.3u(upper = x_0, names = fit_param),
	  start = priors,
	  lowerl = lower,
	  upperl = c(5, min(x_0 + cap, 1), max(concs) * 10),
	  control = controls,
	  na.action = na.omit
	))
	out$fit_type <- 'DRC3pHillFitModelFixS0'
	out$x_0 <- x_0

      } else {
	output_model_new <- try(drc::drm(
	  norm_values ~ concs,
	  data = df_,
	  logDose = NULL,
	  fct = drc::LL.4(names = fit_param),
	  start = priors,
	  lowerl = lower,  
	  upperl = c(5, 1, 1 + cap, max(concs) * 10), 
	  control = controls,
	  na.action = na.omit
	))
	out$fit_type <- 'DRC4pHillFitModel'
      }
    }, error = function(e) {
      out$r2 <- 0
      out$fit_type <- 'DRCInvalidFitResult'
      out$xc50 <- .estimate_xc50(df_$norm_values)
      warning(sprintf("fitting failed with error: '%s'", e))
      return(out)
    })

    # Successful fit.
    for (p in fit_param) { 
      ## drm will output model with the ":(Intercept)" term concatenated at end. 
      out[[p]] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
    }

    out$x_mean <- mean(stats::predict(output_model_new, 
        data.frame(concs = seq(log10(min(df_$concs)), log10(max(df_$concs)), (log10(max(df_$concs)) - log10(min(df_$concs)))/100))), 
        na.rm = TRUE)
    out$x_AOC <- 1 - out$x_mean
    out$x_AOC_range <- 1 - mean(stats::predict(output_model_new, 
      data.frame(concs = seq(log10(range_conc[1]), log10(range_conc[2]), ((log10(range_conc[2]) - log10(range_conc[1]))/100))), 
      na.rm = TRUE))

    # F-test for the significance of the sigmoidal fit.
    RSS2 <- sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
    RSS1 <- sum((df_$norm_values - mean(df_$norm_values,
				      na.rm = TRUE)) ^ 2, na.rm = TRUE)

    out$r2 <- 1 - RSS2 / RSS1

    nparam <- 3 + (is.na(x_0) * 1) # N of parameters in the growth curve; if (x_0 == NA) {4}
    df1 <- nparam - 1 # (N of parameters in the growth curve) - (F-test for the models)
    df2 <- (length(stats::na.omit(df_$norm_values)) - nparam + 1)

    f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
    f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)

    # analytical solution for ic50
    out$xc50 <- out$c50 * ((out$x_0 - out$x_inf) / (0.5 - out$x_inf) - 1) ^
      (1 / out$h)

    # Test the significance of the fit and replace with flat function if required.
    if ((!force) & ((exists("f_pval") & !is.na(f_pval) & f_pval >= pcutoff) | is.na(out$c50))) {
      out <- .set_constant_fit_params(out)
      return(out)
    }

    # Add xc50 = +/-Inf for any curves that do not reach RelativeViability = 0.5
    if (is.na(out$xc50)) {
      out$xc50 <- .estimate_xc50(out$x_inf)
    }

    return(out)
  }


#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  Vinf + (V0 - Vinf) / (1 + (c / EC50) ^ h)
}


logistic_metrics <- function(c, x_metrics) {
  metrics <- c('x_inf', 'x_0', 'h', 'c50')
  if (all(metrics %in% names(x_metrics))) {
    DRC_metrics <- as.vector(x_metrics[metrics])
  } else {
    gr_metric_cols <- get_header('GR_metrics')
    rv_metric_cols <- get_header('RV_metrics')
    if (all(gr_metric_cols %in% names(x_metrics))) {
      metric_cols <- gr_metric_cols
    } else if (all(rv_metric_cols %in% names(x_metrics))) {
      metric_cols <- rv_metric_cols
    } else {
      stop('wrong input parameters')
    }
    DRC_metrics <- as.vector(x_metrics[metric_cols[metrics]])
    names(DRC_metrics) <- metrics
  }

  DRC_metrics$x_inf + (DRC_metrics$x_0 - DRC_metrics$x_inf) /
                          (1 + (c / DRC_metrics$c50) ^ DRC_metrics$h)
}


##################
# Helper functions
##################

# Estimate values for undefined IC/GR50 values. 
#' @keywords internal
.estimate_xc50 <- function(param) {
  if (all(param > 0.5, na.rm = TRUE)) {
    Inf
  } else if (all(param < 0.5, na.rm = TRUE)) {
    -Inf
  } else {
    NA
  }
}
