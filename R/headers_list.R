## The following function utilizes the get_env_identifiers() function which can be
## changed at run time, which is why it needs to be wrapped in a function. 

#' @keywords internal
.getHeadersList <- function() {
  HEADERS_LIST <- list(
    manifest = get_env_identifiers(c("barcode", "template", "duration"), simplify = FALSE),
    raw_data = .getRawDataList(),
    normalized_results = .getNormalizedResultsList(),
    averaged_results = .getAveragedResultsList(),
    response_metrics = .getResponseMetricsList(),
    metric_average_fields = .getMetricAverageFields(),
    # corresponds to the field "celllinename", "primarytissue", "doublingtime" from gneDB CLIDs
    add_clid = get_env_identifiers(c("cellline_name", "cellline_tissue",
                                     "cellline_parental_identifier",
                                     "cellline_subtype", "cellline_ref_div_time"), simplify = FALSE)
  )
  metrics_names <- .getMetricNamesList()
  colnames(metrics_names) <- HEADERS_LIST[["response_metrics"]]
  HEADERS_LIST[["metrics_names"]] <- metrics_names

  HEADERS_LIST[["metrics_results"]] <- c(
    "maxlog10Concentration",
    "N_conc",
    "cotrt_value",
    "source",
    HEADERS_LIST[["response_metrics"]],
    as.character(HEADERS_LIST[["metrics_names"]])
  )

  HEADERS_LIST[["controlled"]] <- c(
    get_env_identifiers("cellline", simplify = TRUE),
    HEADERS_LIST[["manifest"]],
    get_env_identifiers("drug", simplify = TRUE),
    "Concentration",
    paste0(get_env_identifiers("drug", simplify = TRUE), "_", 2:10),
    paste0("Concentration_", 2:10)
  )

  HEADERS_LIST[["reserved"]] <- c(
    HEADERS_LIST[["add_clid"]],
    get_env_identifiers("drug_name", simplify = TRUE),
    get_env_identifiers("masked_tag", simplify = TRUE),
    paste0(get_env_identifiers("drug_name", simplify = TRUE), "_", 2:10),
    paste0(get_env_identifiers("drug_moa", simplify = TRUE), "_", 2:10),
    HEADERS_LIST[["raw_data"]],
    HEADERS_LIST[["normalized_results"]],
    HEADERS_LIST[["averaged_results"]],
    HEADERS_LIST[["metrics_results"]],
    get_env_identifiers("well_position", simplify = TRUE)
  )
  
  HEADERS_LIST[["ordered_1"]] <- .orderHeaderList(HEADERS_LIST, 1)
  HEADERS_LIST[["ordered_2"]] <- .orderHeaderList(HEADERS_LIST, 2)
  
  HEADERS_LIST[["id"]] <- c("rId", "cId")
  
  HEADERS_LIST[["combo"]] <- c("normalization_type",
                               "iso_level",
                               "pos_x",
                               "pos_y",
                               "pos_x_ref",
                               "pos_y_ref",
                               "log2_CI",
                               "log10_ratio_conc")
  
  HEADERS_LIST[["obsolete"]] <- c("RV",
                                  "GR",
                                  "Excess")
  

  HEADERS_LIST
}

#' @keywords internal
.getRawDataList <- function() {
  c(
    "ReadoutValue",
    "BackgroundValue",
    "UntrtReadout",
    "Day0Readout",
    get_env_identifiers("masked_tag", simplify = TRUE)
  )
}

#' @keywords internal
.getNormalizedResultsList <- function() {
  c(
    "x",
    "CorrectedReadout",
    "GRvalue",
    "RelativeViability",
    "DivisionTime",
    "RefGRvalue",
    "RefRelativeViability"
  )
}
  
#' @keywords internal
.getAveragedResultsList <- function() {
  c(
    "x",
    "x_std",
    "std_RelativeViability",
    "std_GRvalue"
  )
}

#' @keywords internal
.getResponseMetricsList <- function() {
  c(
    "x_mean",
    "x_AOC",
    "x_AOC_range",
    "xc50",
    "x_max",
    "ec50",
    "x_inf",
    "x_0",
    "h",
    "r2",
    "x_sd_avg",
    "fit_type"
  )
}

#' @keywords internal
.getMetricNamesList <- function() {
  rbind(
    RV = c(
      "RV_mean",
      "RV_AOC",
      "RV_AOC_range",
      "IC50",
      "E_max",
      "EC50",
      "E_inf",
      "E_0",
      "h_RV",
      "RV_r2",
      "RV_sd_avg",
      "fit_type_RV"
    ),
    GR = c(
      "GR_mean",
      "GR_AOC",
      "GR_AOC_range",
      "GR50",
      "GR_max",
      "GEC50",
      "GR_inf",
      "GR_0",
      "h_GR",
      "GR_r2",
      "GR_sd_avg",
      "fit_type_GR"
    )
  )
}

#' @keywords internal
.getMetricAverageFields <- function() {
  list(
    mean = c(
      "x_mean", 
      "x_AOC", 
      "x_AOC_range", 
      "x_max", 
      "x_inf", 
      "x_0"
    ),
    geometric_mean = c(
      "xc50", 
      "ec50"
    )
  )
}

#' @keywords internal
.orderHeaderList <- function(headers_list, type) {
  type <- match.arg(as.character(type), c(1, 2))
  if (type == 1) {
    c(
      get_env_identifiers("cellline_name", simplify = TRUE),
      get_env_identifiers("cellline_tissue", simplify = TRUE),
      get_env_identifiers("duration", simplify = TRUE),
      get_env_identifiers("drug_name", simplify = TRUE),
      "Concentration",
      paste0(c(paste0(get_env_identifiers("drug_name", simplify = TRUE), "_"), "Concentration_"), 
             rep(2:10, each = 2))
    )
  } else {
    c(
      headers_list[["normalized_results"]],
      headers_list[["averaged_results"]],
      headers_list[["metrics_results"]],
      headers_list[["raw_data"]],
      get_env_identifiers("cellline_ref_div_time", simplify = TRUE),
      get_env_identifiers("cellline", simplify = TRUE),
      get_env_identifiers("drug", simplify = TRUE),
      paste0(get_env_identifiers("drug", simplify = TRUE), "_", 2:10),
      headers_list[["manifest"]],
      get_env_identifiers("well_position", simplify = TRUE)
    )
  }
}
