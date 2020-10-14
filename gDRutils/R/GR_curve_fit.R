#' Actual fitting function
#'
#' \code{RVGRfits} returns fit parameters
#'
#' returns fit parameters
#'
#' @import reshape2
#' @param log10concs concentrations
#' @param RelativeViability values
#' @param GRvalues values
#' @param e_0
#' Defaults to \code{1}.
#' @param GR_0
#' Defaults to \code{1}.
#' @param force use signifcance or not
#' @param cap enforce e_0 and GR_0
#' @return vector of parameters
#' @examples
#' @importFrom drc drm drmc LL.3u LL.4
#' @export
RVGRfits <- function(df_,
                     e_0 = 1,
                     GR_0 = 1,
                     n_point_cutoff = 4,
                     range_conc = c(5e-3, 5),
                     force = FALSE) {

  df_RV = gDRutils::logisticFit(
    df_$Concentration,
    df_$RelativeViability,
    df_$std_RelativeViability,
    x_0 = e_0,
    curve_type = "RV",
    range_conc = range_conc,
    force = force,
    n_point_cutoff = n_point_cutoff
  )

  df_GR = gDRutils::logisticFit(
    df_$Concentration,
    df_$GRvalue,
    df_$std_GRvalue,
    x_0 = GR_0,
    curve_type = "GR",
    range_conc = range_conc,
    force = force,
    n_point_cutoff = n_point_cutoff
  )

  df_metrics = rbind(df_RV, df_GR)
  rownames(df_metrics) <- c("RV", "GR")

  return(df_metrics)
}


#' Actual fitting function
#'
#' \code{logisticFit} returns fit parameters
#'
#' @param log10concs log10 of concentrations
#' @param normValues normalized response values (Untreated = 1)
#' @param std_normValues std of values
#' @param x_0 upper limit
#' Defaults to \code{1}. 
#' @param curve_type response curve: either RV ([0,1]) or GR([-1,1])
#' @param range_conc range of concentration for calculating AOC_range
#' @param force use signifcance or not
#' @param cap enforce values stay below (x_0+cap)
#' @return DataFrame with metrics and fit parameters
#' @examples
#' @details
#' Implementation of the genedata approach for curve fit: https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
#' @import reshape2
#' @importFrom drc drm drmc LL.3u
#' @export
logisticFit <-
  function(concs,
           normValues,
           std_normValues,
           x_0 = 1,
           curve_type = c("RV", "GR"),
           range_conc = c(5e-3, 5),
           force = FALSE,
           cap = 0.1,
           n_point_cutoff = 4) {

    # define variables and prepare data
    log10concs <- log10(concs)
    df_ <- data.frame(log10conc = log10concs,
                     normValues = pmin(normValues, ifelse(is.na(x_0),1,x_0) + cap))

    out <- matrix(NA, 1, length(gDRutils::get_header("response_metrics")))
    colnames(out) <- gDRutils::get_header("response_metrics")
    # transform to a matrix to allow for string (for fitting model name)
    out = S4Vectors::DataFrame(out)
    out$maxlog10Concentration <- max(log10concs)
    out$N_conc <- length(unique(log10concs))
    out$x_sd_avg <- mean(std_normValues, rm.na=T)

    # non-fitted metrics
    xAvg <- aggregate(
      df_$normValues,
      by = list(log10conc = df_$log10conc),
      FUN = function(x)
        mean(x, na.rm = T)
    )
    colnames(xAvg)[2] <- "normValues"
    l <- dim(xAvg)[1]

    out$x_max <- min(xAvg$normValues[c(l, l - 1)], na.rm = TRUE)
    # temp values if fit  fails
    out$x_mean <- mean(xAvg$normValues, na.rm = T)
    out$x_AOC <- 1 - out$x_mean

    if (length(unique(xAvg$normValues[!is.na(xAvg$normValues)])) == 1) {
      out$fit_type == 'DRCConstantFitResult'
      out$c50 <- 0
      out$h <- 0.0001
      out$xc50 <- ifelse(mean(xAvg$normValues, na.rm = T) > .5, Inf,-Inf)
      out$x_inf <- out$x_mean <- mean(xAvg$normValues, rm.na = T)
      out$x_AOC_range <- out$x_AOC <- 1 - mean(xAvg$normValues, na.rm = T)
      return(out)
    }

    if (sum(!is.na(xAvg$normValues)) < n_point_cutoff) {
      out$fit_type = 'DRCTooFewPointsToFit'
      # best estimate if the data cannot be fit
      out$xc50 <- ifelse(all(df_$normValues > .5, na.rm = T), Inf,
                    ifelse(all(df_$normValues < .5, na.rm = T), -Inf,
                      NA))
      return(out)
    }

    # fit parameters and boundaries
    fit_para <- c("h", "x_inf", "c50")
    if (curve_type == "RV") {
      priors <- c(2, 0.4, median(concs))
      lower <- c(.1, 0, min(concs) / 10)
    } else if (curve_type == "GR") {
      priors <- c(2, 0.1, median(concs))
      lower <- c(.1,-1, min(concs) / 10)
    }

    controls <- drc::drmc()
    controls$relTol <- 1e-06
    controls$errorm <- FALSE
    controls$noMessage <- TRUE
    controls$rmNA <- TRUE

    ######################################
    # curve fitting
    if (!is.na(x_0)) {
      output_model_new = try(drc::drm(
        normValues ~ log10conc,
        data = df_,
        logDose = 10,
        fct = drc::LL.3u(upper = x_0, names = fit_para),
        start = priors,
        lowerl = lower,
        upperl = c(5, min(x_0 + .1, 1), max(concs) * 10),
        control = controls,
        na.action = na.omit
      ))
      out$fit_type = 'DRC3pHillFitModelFixS0'
      out$x_0 <- x_0
    } else {
      fit_para <- c(fit_para[1:2], 'x_0', fit_para[3])
      output_model_new = try(drc::drm(
        normValues ~ log10conc,
        data = df_,
        logDose = 10,
        fct = drc::LL.4(names = fit_para),
        start = c(priors[1:2], 1, priors[3]),
        lowerl = lower[c(1:2, 2, 3)],
        upperl = c(5, 1, 1 + cap, max(concs) * 10),
        control = controls,
        na.action = na.omit
      ))
      out$fit_type = 'DRC4pHillFitModel'
    }
    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in fit_para) {
        out[p] = stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      out$x_mean = mean(stats::predict(output_model_new, data.frame(
          concs = seq(min(df_$log10conc), max(df_$log10conc), .03))), na.rm = T)
      out$x_AOC = 1 - out$x_mean
      out$x_AOC_range = 1 - mean(stats::predict(output_model_new, data.frame(
          concs = seq(log10(range_conc[1]), log10(range_conc[2]), .03))), na.rm = T)
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 + (is.na(x_0)*1) # N of parameters in the growth curve; if x_0 == NA -> 4
      Npara_flat <- 1 # F-test for the models
      RSS2 <- sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <- sum((df_$normValues - mean(df_$normValues,
                                        na.rm = TRUE)) ^ 2, na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <- (length(na.omit(df_$normValues)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out$r2 <- 1 - RSS2 / RSS1

      # analytical solution for ic50
      out$xc50 <- out$c50 * ((out$x_0 - out$x_inf) / (0.5 - out$x_inf) - 1) ^
        (1 / out$h)
    } else {
      # fit error
      out$r2 = 0
      out$fit_type = 'DRCInvalidFitResult'
      # best estimate if the data cannot be fit
      out$xc50 <- ifelse(all(df_$normValues > .5, na.rm = T), Inf,
                    ifelse(all(df_$normValues < .5, na.rm = T), -Inf,
                      NA))
      return(out)
    }

    # testing the significance of the fit and replacing with flat function if required
    pcutoff <- ifelse(force, 1, .05)
    if (!is.na(f_pval)) {
      out$fit_type <- ifelse(f_pval >= pcutoff |
                                 is.na(out$c50), 'DRCConstantFitResult', out$fit_type)
    } else {
      out$fit_type <- ifelse(is.na(out$c50), 'DRCConstantFitResult', out$fit_type)
    }

    # Replace values for flat fits: c50 = 0, h = 0.01 and xc50 = +/- Inf
    if (out$fit_type == 'DRCConstantFitResult') {
      out$c50 <- 0
      out$h <- 0.0001
      out$xc50 <- ifelse(mean(xAvg$normValues, na.rm = T) > .5, Inf,-Inf)
      out$x_inf <- out$x_mean <- mean(xAvg$normValues, na.rm = T)
      out$x_AOC_range <- out$x_AOC <- 1 - mean(xAvg$normValues, na.rm = T)
    }

    # Add xc50 = +/-Inf for any curves that don"t reach RelativeViability = 0.5
    if (is.na(out$xc50)) {
      out$xc50 <- ifelse(out$x_inf > .5, Inf,-Inf)
    }
    return(out)
  }


# TODO: replace by drc::LL.4
# logistic function (not used in the file but useful for plotting externally)
#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  Vinf + (V0 - Vinf) / (1 + (c / EC50) ^ h)
}

logistic_metrics <- function(c, x_metrics) {
  metrics = c('x_inf', 'x_0', 'h', 'c50')
  if (all(metrics %in% names(x_metrics))) DRC_metrics = as.vector(x_metrics[metrics])
  else if (all(get_header('GR_metrics')[metrics] %in% names(x_metrics))) {
    DRC_metrics = as.vector(x_metrics[get_header('GR_metrics')[metrics]])
    names(DRC_metrics) = metrics
  } else if (all(get_header('IC_metrics')[metrics] %in% names(x_metrics))) {
    DRC_metrics = as.vector(x_metrics[get_header('IC_metrics')[metrics]])
    names(DRC_metrics) = metrics
  } else stop('wrong input parameters')

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
#' @param log10concs concentrations
#' @param RelativeViability values
#' @param GRvalues values
#' @param e_0 =1 by default
#' @param GR_0 =1 by default
#' @param force use signifcance or not
#' @param cap enforce e_0 and GR_0
#' @return vector of parameters
#' @examples
#' sum(1:10)
#' @import reshape2
#' @importFrom drc drm drmc LL.3u
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
    ICpriors <- c(2, 0.4, median(concs))
    GRpriors <- c(2, 0.1, median(concs))
    IClower <- c(.1, 0, min(concs) / 10)
    GRlower <- c(.1, -1, min(concs) / 10)
    upper <- c(5, 1, max(concs) * 10)

    controls <- drc::drmc()
    controls$relTol <- 1e-06
    controls$errorm <- FALSE
    controls$noMessage <- TRUE
    controls$rmNA <- TRUE

    ######################################
    # IC curve fitting
    output_model_new <- try(drc::drm(
      RelativeViability ~ log10conc,
      data = IC_data_exp,
      logDose = 10,
      fct = drc::LL.3u(upper = e_0, names = ICfit_parameters),
      start = ICpriors,
      lowerl = IClower,
      upperl = upper,
      control = controls,
      na.action = na.omit
    ))

    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in ICfit_parameters) {
        out[p] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 # N of parameters in the growth curve
      Npara_flat <- 1 # F-test for the models
      RSS2 <-
        sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <-
        sum((
          IC_data_exp$RelativeViability - mean(IC_data_exp$RelativeViability,
                                               na.rm = TRUE)
        ) ^ 2, na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <-
        (length(na.omit(IC_data_exp$RelativeViability)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out["rv_r2"] <- 1 - RSS2 / RSS1
    }


    # non-fitted metrics
    ICavg <-
      aggregate(
        IC_data_exp$RelativeViability,
        by = list(log10conc = IC_data_exp$log10conc),
        FUN = mean
      )
    colnames(ICavg)[2] <- "RelativeViability"
    l <- dim(ICavg)[1]

    out["e_max"] <-
      min(ICavg$RelativeViability[c(l, l - 1)], na.rm = TRUE)

    out["rv_mean"] <- mean(ICavg$RelativeViability)

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
    pcutoff <- ifelse(force, 1, .05)
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
      out["GR50"] <- ifelse(mean(GRavg$GRvalue) > .5, Inf, -Inf)
      out["GRinf"] <- mean(GRavg$GRvalue)
    }

    # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
    if (is.na(out["GR50"])) {
      out["GR50"] <- ifelse(out["GRinf"] > .5, Inf, -Inf)
    }

    out
  }
