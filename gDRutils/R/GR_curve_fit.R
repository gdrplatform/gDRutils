#' Fit curves
#'
#' Fit GR and RV curves from a data.frame.
#'
#' @import reshape2
#' @param df_ data.frame containing \code{Concentration}, normalized values, and standardized normalized values
#' @param e_0
#' Defaults to \code{1}.
#' @param GR_0
#' Defaults to \code{1}.
#' @param force boolean indicating whether or not to force a constant fit
#' @param pcutoff numeric of pvalue significance threshold above which to use a constant fit
#'
#' @return data.frame of RV and GR fit parameters.
#'
#' @examples
#'
#' @export
#'
RVGRfits <- function(df_,
		      e_0 = 1,
		      GR_0 = 1,
		      n_point_cutoff = 4,
		      range_conc = c(5e-3, 5),
		      force = FALSE, 
                      pcutoff = 0.05, 
                      perform_log = TRUE) {

  df_RV <- gDRutils::logisticFit(
    df_$Concentration,
    df_$RelativeViability,
    df_$std_RelativeViability,
    x_0 = e_0,
    curve_type = "RV",
    range_conc = range_conc,
    force = force,
    pcutoff = pcutoff, 
    n_point_cutoff = n_point_cutoff, 
    perform_log = perform_log
  )

  df_GR <- gDRutils::logisticFit(
    df_$Concentration,
    df_$GRvalue,
    df_$std_GRvalue,
    x_0 = GR_0,
    curve_type = "GR",
    range_conc = range_conc,
    force = force,
    pcutoff = pcutoff, 
    n_point_cutoff = n_point_cutoff,
    perform_log = perform_log
  )

  df_metrics <- rbind(df_RV, df_GR)
  rownames(df_metrics) <- c("RV", "GR")

  return(df_metrics)
}


#' Logistic fit
#'
#' \code{logisticFit} returns fit parameters
#'
#' @import reshape2
#' @param concs concentrations
#' @param norm_values normalized response values (Untreated = 1)
#' @param std_norm_values std of values
#' @param x_0 upper limit
#' Defaults to \code{1}. For co-treatments, this value should be set to \code{NA}.
#' @param curve_type response curve: either RV ([0,1]) or GR([-1,1])
#' @param range_conc range of concentration for calculating AOC_range
#' @param force boolean indicating whether or not to force a constant fit
#' @param pcutoff numeric of pvalue significance threshold above which to use a constant fit
#' @param cap numeric value capping \code{norm_values} to stay below (\code{x_0} + cap)
#' @param n_point_cutoff integer indicating number of unique concentrations required to fit curve
#' @param perform_log boolean indicating whether or not to log10 transform the \code{concs}
#'
#' @return Named list with metrics and fit parameters
#'
#' @details
#' Implementation of the genedata approach for curve fit: 
#' https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
#'
#' \itemize{
#'  \item{x_mean}{The mean of a given dose-response metric}
#'  \item{x_AOC_range}{The range of the area over the curve}
#'  \item{x_AOC}{The area over the GR curve or, respectively, under the relative cell count curve, averaged over the range of concentration values}
#'  \item{xc50}{The concentration at which the effect reaches a value of 0.5 based on interpolation of the fitted curve}
#'  \item{x_max}{The maximum effect of the drug}
#'  \item{c50}{The drug concentration at half-maximal effect}
#'  \item{x_inf}{The asymptotic value of the sigmoidal fit to the dose-response data as concentration goes to infinity}
#'  \item{x_0}{The asymptotic metric value corresponding to a concentration of 0 for the primary drug}
#'  \item{h}{}
#'  \item{r2}{}
#'  \item{x_sd_avg}{}
#'  \item{fit_type}{This will be given by one of the following:
#'    \item{"DRC4pHillFitModel"}{}
#'    \item{"DRC3pHillFitModelFixS0"}{}
#'    \item{"DRCConstantFitResult"}{}
#'    \item{"DRCTooFewPointsToFit"}{}
#'    \item{"DRCInvalidFitResult"}{}
#'  }
#'  \item{maxlog10Concentration}{}
#'  \item{N_conc}{Number of unique concentrations}
#' }
#' @export
#'
logisticFit <-
  function(concs,
           norm_values,
           std_norm_values,
           x_0 = 1,
           curve_type = c("RV", "GR"),
           range_conc = c(5e-3, 5),
           force = FALSE,
           pcutoff = 0.05,
           cap = 0.1,
           n_point_cutoff = 4, 
           perform_log = TRUE) {

    resp_metric_cols <- gDRutils::get_header("response_metrics")
    out <- vector("list", length(resp_metric_cols))
    names(out) <- resp_metric_cols

    if (perform_log) {
      log10concs <- log10(concs)
    } else {
      log10concs <- concs
    } 

    # Cap norm_values at (x_0 + cap).
    norm_values <- pmin(norm_values, (ifelse(is.na(x_0), 1, x_0) + cap))

    out$maxlog10Concentration <- max(log10concs)
    out$N_conc <- length(unique(log10concs))
    out$x_sd_avg <- mean(std_norm_values, na.rm = TRUE)

    df_ <- data.frame(log10conc = log10concs,
                      norm_values = norm_values)

    # Calculate metrics that do not require fitting.
    xAvg <- aggregate(
      list(norm_values = norm_values),
      by = list(log10conc = df_$log10conc),
      FUN = function(x) {mean(x, na.rm = TRUE)}
    )
    l <- nrow(xAvg) # NAs have been removed. 

    mean_norm_value <- mean(xAvg$norm_values, na.rm = TRUE)

    ## 'x_max' can be considered either the lowest readout (max efficacy) 
    ## or the efficacy at the max concentration. We take the min 
    ## of the two highest concentrations as a compromise.
    out$x_max <- min(xAvg$norm_values[c(l, l - 1)], na.rm = TRUE)

    # Set temp values if fit fails.
    out$x_mean <- mean_norm_value
    out$x_AOC <- 1 - out$x_mean

    # Best estimate if the data cannot be fit.
    .estimate_xc50 <- function(norm_val) {
      if (all(norm_val > 0.5, na.rm = TRUE)) {
        Inf
      } else if (all(norm_val < 0.5, na.rm = TRUE)) {
        -Inf
      } else {
        NA
      }
    }

    non_na_avg_norm <- !is.na(xAvg$norm_values)

    if (length(unique(xAvg$norm_values[non_na_avg_norm])) == 1L) {
      out$fit_type <- 'DRCConstantFitResult'
      out$x_0 <- x_0
      out$c50 <- 0
      out$h <- 0.0001
      out$xc50 <- ifelse(mean_norm_value > 0.5, Inf, -Inf)
      out$x_inf <- out$x_mean <- mean_norm_value
      out$x_AOC_range <- out$x_AOC <- 1 - mean_norm_value
      return(out)
    }

    if (sum(non_na_avg_norm) < n_point_cutoff) {
      out$fit_type <- 'DRCTooFewPointsToFit'

      out$xc50 <- .estimate_xc50(norm_values)
      return(out)
    }

    # Set fit parameters and boundaries.
    fit_param <- c("h", "x_inf", "x_0", "c50")

    if (curve_type == "RV") {
      priors <- c(2, 0.4, 1, median(concs))
      lower <- c(0.1, 0, 0, min(concs) / 10)
    } else if (curve_type == "GR") {
      priors <- c(2, 0.1, 1, median(concs))
      lower <- c(0.1, -1, -1, min(concs) / 10)
    }

    controls <- drc::drmc(relTol = 1e-06, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)

    # Curve fitting.
    if (!is.na(x_0)) { # For co-treatments.
      # Override existing params for removal of x_0 parameter. 
      fit_param <- fit_param[-3]
      priors <- priors[-3]
      lower <- lower[-3]

      output_model_new <- try(drc::drm(
        norm_values ~ log10conc,
        data = df_,
        logDose = 10,
        fct = drc::LL.3u(upper = x_0, names = fit_param),
        start = priors,
        lowerl = lower,
        upperl = c(5, min(x_0 + 0.1, 1), max(concs) * 10),
        control = controls,
        na.action = na.omit
      ))
      out$fit_type <- 'DRC3pHillFitModelFixS0'
      out$x_0 <- x_0

    } else {
      output_model_new <- try(drc::drm(
        norm_values ~ log10conc,
        data = df_,
        logDose = 10,
        fct = drc::LL.4(names = fit_param),
        start = priors,
        lowerl = lower,  
        upperl = c(5, 1, 1 + cap, max(concs) * 10),
        control = controls,
        na.action = na.omit
      ))
      out$fit_type <- 'DRC4pHillFitModel'
    }

    # Successful fit.
    if (class(output_model_new) != "try-error") {
      for (p in fit_param) { 
        out[[p]] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }

      out$x_mean <- mean(stats::predict(output_model_new, data.frame(
          concs = seq(min(df_$log10conc), max(df_$log10conc), 0.03))), na.rm = TRUE)
      out$x_AOC <- 1 - out$x_mean
      out$x_AOC_range <- 1 - mean(stats::predict(output_model_new, data.frame(
          concs = seq(log10(range_conc[1]), log10(range_conc[2]), 0.03))), na.rm = TRUE)

      # F-test for the significance of the sigmoidal fit.
      RSS2 <- sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <- sum((df_$norm_values - mean(df_$norm_values,
                                        na.rm = TRUE)) ^ 2, na.rm = TRUE)

      out$r2 <- 1 - RSS2 / RSS1

      Npara <- 3 + (is.na(x_0)*1) # N of parameters in the growth curve; if x_0 == NA -> 4
      df1 <- Npara - 1 # (N of parameters in the growth curve) - (F-test for the models)
      df2 <- (length(na.omit(df_$norm_values)) - Npara + 1)

      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)

      # analytical solution for ic50
      out$xc50 <- out$c50 * ((out$x_0 - out$x_inf) / (0.5 - out$x_inf) - 1) ^
        (1 / out$h)

    } else { # Failed fit.
      out$r2 <- 0
      out$fit_type <- 'DRCInvalidFitResult'

      out$xc50 <- .estimate_xc50(df_$norm_values)

      return(out)
    }

    # Test the significance of the fit and replace with flat function if required.
    pcutoff <- ifelse(force, 1, pcutoff)
    out$fit_type <- ifelse((!is.na(f_pval) & f_pval >= pcutoff) | is.na(out$c50), 
			   'DRCConstantFitResult', 
			   out$fit_type)

    # Replace values for flat fits: c50 = 0, h = 0.01 and xc50 = +/- Inf
    if (out$fit_type == 'DRCConstantFitResult') {
      out$c50 <- 0
      out$h <- 0.0001
      # TODO: I think this is missing an x_0 value, even if it is NA.

      out$xc50 <- ifelse(mean_norm_value > 0.5, Inf, -Inf)
      out$x_inf <- out$x_mean <- mean_norm_value
      out$x_AOC_range <- out$x_AOC <- 1 - mean_norm_value
    }

    # Add xc50 = +/-Inf for any curves that do not reach RelativeViability = 0.5
    if (is.na(out$xc50)) {
      out$xc50 <- ifelse(out$x_inf > 0.5, Inf, -Inf)
    }

    return(data.frame(out))
  }


# TODO: replace by drc::LL.4
# logistic function (not used in the file but useful for plotting externally)
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


# TODO: remove this function once there are no dependencies in gDRcore
#' Actual fitting function
#'
#' \code{ICGRlogisticFit} returns fit parameters
#'
#' returns fit parameters
#'
#' @import reshape2
#' @param log10concs concentrations
#' @param RelativeViability values
#' @param GRvalues values
#' @param e_0 =1 by default
#' @param GR_0 =1 by default
#' @param force use signifcance or not
#' @param cap enforce e_0 and GR_0
#' @return vector of parameters
#' @examples
#' @export
ICGRlogisticFit <-
  function(log10concs,
           RelativeViability,
           GRvalues,
           e_0 = 1,
           GR_0 = 1,
           force = FALSE,
           cap = FALSE) {
    # Implementation of the genedata approach for curve fit: https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
    #
    .Deprecated("RVGRfits", package = "gDRutils")

    # define variables and prepare data
    IC_data_exp <-
      data.frame(log10conc = log10concs, RelativeViability = RelativeViability)
    GR_data_exp <-
      data.frame(log10conc = log10concs, GRvalue = GRvalues)
    concs <- 10 ** log10concs
    ICfit_parameters <- c("h_rv", "e_inf", "ec50")
    GRfit_parameters <- c("h_GR", "GRinf", "GEC50")

    out <- array(NA, length(get_header("metrics_results")))
    names(out) <- get_header("metrics_results")
    out["maxlog10Concentration"] <- max(log10concs)
    out["N_conc"] <- length(unique(log10concs))
    out["e_0"] <- e_0
    out["GR_0"] <- GR_0

    # fit parameters and boundaries
    RVpriors <- c(2, 0.4, median(concs))
    GRpriors <- c(2, 0.1, median(concs))
    RVlower <- c(.1, 0, min(concs) / 10)
    GRlower <- c(.1, -1, min(concs) / 10)
    upper <- c(5, 1, max(concs) * 10)

    controls <- drc::drmc()
    controls$relTol <- 1e-06
    controls$errorm <- FALSE
    controls$noMessage <- TRUE
    controls$rmNA <- TRUE

    ######################################
    # RV curve fitting
    output_model_new <- try(drc::drm(
      RelativeViability ~ log10conc,
      data = RV_data_exp,
      logDose = 10,
      fct = drc::LL.3u(upper = e_0, names = RVfit_parameters),
      start = RVpriors,
      lowerl = RVlower,
      upperl = upper,
      control = controls,
      na.action = na.omit
    ))

    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in RVfit_parameters) {
        out[p] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 # N of parameters in the growth curve
      Npara_flat <- 1 # F-test for the models
      RSS2 <-
        sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <-
        sum((
          RV_data_exp$RelativeViability - mean(RV_data_exp$RelativeViability,
                                               na.rm = TRUE)
        ) ^ 2, na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <-
        (length(na.omit(RV_data_exp$RelativeViability)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out["rv_r2"] <- 1 - RSS2 / RSS1
    }


    # non-fitted metrics
    RVavg <-
      aggregate(
        RV_data_exp$RelativeViability,
        by = list(log10conc = RV_data_exp$log10conc),
        FUN = mean
      )
    colnames(RVavg)[2] <- "RelativeViability"
    l <- dim(RVavg)[1]

    out["e_max"] <-
      min(RVavg$RelativeViability[c(l, l - 1)], na.rm = TRUE)

    out["rv_mean"] <- mean(RVavg$RelativeViability)

    # analytical solution for ic50
    out["ic50"] <-
      out["ec50"] * ((e_0 - out["e_inf"]) / (0.5 - out["e_inf"]) - 1) ^ (1 /
                                                                           out["h_rv"])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff <- ifelse(force, 1, .05)
    if (!is.na(f_pval)) {
      out["fit_model_rv"] <- ifelse(f_pval >= pcutoff |
                                     is.na(out["ec50"]), 1, 0)
    } else {
      out["fit_model_rv"] <- ifelse(is.na(out["ec50"]), 1, 0)
    }

    # Replace values for flat fits: ec50 = 0, h_rv = 0.01 and ic50 = +/- Inf
    if (out["fit_model_rv"] == 1) {
      out["ec50"] <- 0
      out["h_rv"] <- 0.0001
      out["ic50"] <-
        ifelse(mean(ICavg$RelativeViability) > .5, Inf, -Inf)
      out["e_inf"] <- mean(ICavg$RelativeViability)
    }

    # Add ic50 = +/-Inf for any curves that don"t reach RelativeViability = 0.5
    if (is.na(out["ic50"])) {
      out["ic50"] <- ifelse(out["e_inf"] > .5, Inf, -Inf)
    }

    ######################################
    # GR curve fitting
    output_model_new <- try(drc::drm(
      GRvalue ~ log10conc,
      data = GR_data_exp,
      logDose = 10,
      fct = drc::LL.3u(upper = GR_0, names = GRfit_parameters),
      start = GRpriors,
      lowerl = GRlower,
      upperl = upper,
      control = controls,
      na.action = na.omit
    ))

    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in GRfit_parameters) {
        out[p] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 # N of parameters in the growth curve
      Npara_flat <- 1 # F-test for the models
      RSS2 <-
        sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <-
        sum((
          GR_data_exp$GRvalue - mean(GR_data_exp$GRvalue, na.rm = TRUE)
        ) ^ 2,
        na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <- (length(na.omit(GR_data_exp$GRvalue)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out["GR_r2"] <- 1 - RSS2 / RSS1
    }

    # non-fitted metrics
    GRavg <-
      aggregate(
        GR_data_exp$GRvalue,
        by = list(log10conc = GR_data_exp$log10conc),
        FUN = mean
      )
    colnames(GRavg)[2] <- "GRvalue"
    l <- dim(GRavg)[1]

    out["GRmax"] <- min(GRavg$GRvalue[c(l, l - 1)], na.rm = TRUE)

    out["GR_AOC"] <-
      mean(1 - GRavg$GRvalue) # use mean for consistency with mean.viability
    # ## Alternative version: Trapezoid rule for integration of GR_AOC
    # diff_vector <- diff(GRavg$log10conc, lag = 1)
    # conc_range <- GRavg$log10conc[l] - GRavg$log10conc[1]
    # out["GR_AOC"] <- sum((1 - (GRavg$GRvalue[1:(l-1)]+GRavg$GRvalue[2:l])/2)*
    #                diff_vector, na.rm = TRUE)/conc_range

    # analytical solution for GR50
    out["GR50"] <-
      out["GEC50"] * ((GR_0 - out["GRinf"]) / (0.5 - out["GRinf"]) - 1) ^ (1 /
                                                                             out["h_GR"])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff <- ifelse(force, 1, 0.05)
    if (!is.na(f_pval)) {
      out["flat_fit_GR"] <- ifelse(f_pval >= pcutoff |
                                     is.na(out["GEC50"]), 1, 0)
    } else {
      out["flat_fit_GR"] <- ifelse(is.na(out["GEC50"]), 1, 0)
    }

    # Replace values for flat fits: GEC50 = 0, h_GR = 0.01 and GR50 = +/- Inf
    if (out["flat_fit_GR"] == 1) {
      out["GEC50"] <- 0
      out["h_GR"] <- 0.0001
      out["GR50"] <- ifelse(mean(GRavg$GRvalue) > 0.5, Inf, -Inf)
      out["GRinf"] <- mean(GRavg$GRvalue)
    }

    # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
    if (is.na(out["GR50"])) {
      out["GR50"] <- ifelse(out["GRinf"] > 0.5, Inf, -Inf)
    }

    data.frame(out)
  }
