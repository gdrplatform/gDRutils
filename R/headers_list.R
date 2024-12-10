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
    "maxlog10Concentration_sd",
    "N_conc",
    "N_conc_sd", 
    "cotrt_value",
    "source",
    "count",
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
  
  
  HEADERS_LIST[["iso_position"]] <- c("iso_level",
                                      "pos_x",
                                      "pos_y",
                                      "pos_x_ref",
                                      "pos_y_ref")
  
  HEADERS_LIST[["excess"]] <- names(get_combo_excess_field_names())
  
  
  HEADERS_LIST[["scores"]] <- names(get_combo_score_field_names())
  
  HEADERS_LIST[["isobolograms"]] <- c("normalization_type",
                                      HEADERS_LIST[["iso_position"]],
                                      "log2_CI",
                                      "log10_ratio_conc")
  
  HEADERS_LIST[["fit_source"]] <- "fit_source"
  
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
    "std_GRvalue",
    # after averaging for biological replicates
    "count",
    "x_sd",
    "x_std_sd"
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
    "p_value",
    "rss",
    "x_sd_avg",
    "fit_type",
    "x_mean_sd",
    "x_AOC_sd",
    "x_AOC_range_sd",
    "xc50_sd",
    "x_max_sd",
    "ec50_sd",
    "x_inf_sd",
    "x_0_sd",
    "h_sd",
    "r2_sd",
    "x_sd_avg_sd"
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
      "RV_p_value",
      "RV_rss",
      "RV_sd_avg",
      "fit_type_RV",
      "RV_mean_sd",
      "RV_AOC_sd",
      "RV_AOC_range_sd",
      "IC50_sd",
      "E_max_sd",
      "EC50_sd",
      "E_inf_sd",
      "E_0_sd",
      "h_RV_sd",
      "RV_r2_sd",
      "RV_sd_avg_sd"
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
      "GR_p_value",
      "GR_rss",
      "GR_sd_avg",
      "fit_type_GR",
      "GR_mean_sd",
      "GR_AOC_sd",
      "GR_AOC_range",
      "GR50_sd",
      "GR_max_sd",
      "GEC50_sd",
      "GR_inf_sd",
      "GR_0_sd",
      "h_GR_sd",
      "GR_r2_sd",
      "GR_sd_avg_sd"
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
      "ec50",
      "GR50",
      "GEC50",
      "IC50",
      "EC50",
      "GR_xc50",
      "RV_xc50",
      "GR_ec50",
      "RV_ec50"
    ),
    fit_type = c(
      "fit_type",
      "Fit Type",
      "Fit Type RV",
      "Fit Type GR",
      "RV_fit_type",
      "GR_fit_type"
   ),
   # due to the fact that there is some freedom in what values are in individual fields, 
   # in order to avoid duplicates in the application we have to exclude some fields from 
   # recognizing duplicates in averaging
   blacklisted = c(
     # tissue
     "cellline_tissue",
     "Tissue", # sometimes this field is missing
     # reference division time
     "cellline_ref_div_time",
     "Reference Division Time", # sometimes this field has `NA`s
     "ReferenceDivisionTime",
     # parental identifier
     "cellline_parental_identifier",
     "Parental Identifier", # sometimes suffixes incorrectly differentiate this field
     "parental_identifier" # sometimes suffixes incorrectly differentiate this field
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
