## The following function utilizes the get_identifier() function which can be
## changed at run time, which is why it needs to be wrapped in a function. 

#' @keywords internal
.getHeadersList <- function() {
  HEADERS_LIST <- list(
    manifest = c(
      "Barcode", 
      "Template", 
      get_identifier("duration")
    ),

    raw_data = c(
      "ReadoutValue",
      "BackgroundValue",
      "UntrtReadout",
      "Day0Readout",
      get_identifier('masked_tag')
    ),

    normalized_results = c(
      "CorrectedReadout",
      "GRvalue",
      "RelativeViability",
      "DivisionTime",
      "RefGRvalue",
      "RefRelativeViability"
    ),

    averaged_results = c(
      "std_GRvalue", 
      "std_RelativeViability"
    ),

    response_metrics = c(
      "x_mean",
      "x_AOC",
      "x_AOC_range",
      "xc50",
      "x_max",
      "c50",
      "x_inf",
      "x_0",
      "h",
      "r2",
      "x_sd_avg",
      "fit_type"
    ),

    # corresponds to the field "celllinename", "primarytissue", "doublingtime" from gneDB CLIDs
    add_clid = c(
      "CellLineName",
      "Tissue",
      "ReferenceDivisionTime"
    )
  )

  rv_metrics <- c(
    "RV_mean",
    "RV_AOC",
    "RV_AOC_range",
    "ic50",
    "e_max",
    "ec50",
    "e_inf",
    "e_0",
    "h_RV",
    "RV_r2",
    "RV_sd_avg",
    "fit_type_RV"
  )
  names(rv_metrics) <- HEADERS_LIST[["response_metrics"]]
  HEADERS_LIST[["RV_metrics"]] <- rv_metrics

  gr_metrics <- c(
    "mean_GR",
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
    "flat_fit_GR"
  )
  names(gr_metrics) <- HEADERS_LIST[["response_metrics"]]
  HEADERS_LIST[["GR_metrics"]] <- gr_metrics

  HEADERS_LIST[["metrics_results"]] <- c(
    "maxlog10Concentration",
    "N_conc",
    HEADERS_LIST[["response_metrics"]],
    HEADERS_LIST[["RV_metrics"]],
    HEADERS_LIST[["GR_metrics"]]
  )

  HEADERS_LIST[["controlled"]] <- c(
    get_identifier("cellline"),
    HEADERS_LIST[["manifest"]],
    get_identifier("drug"),
    "Concentration",
    paste0(get_identifier("drug"), "_", 2:10),
    paste0("Concentration_", 2:10)
  )

  HEADERS_LIST[["reserved"]] <- c(
    HEADERS_LIST[["add_clid"]],
    get_identifier("drugname"),
    get_identifier("masked_tag"),
    paste0(get_identifier("drugname"), "_", 2:10),
    HEADERS_LIST[["raw_data"]],
    HEADERS_LIST[["normalized_results"]],
    HEADERS_LIST[["averaged_results"]],
    HEADERS_LIST[["metrics_results"]],
    "WellRow",
    "WellColumn"
  )

  HEADERS_LIST[["ordered_1"]] <- c(
    HEADERS_LIST[["add_clid"]][1:2],
    get_identifier("duration"),
    get_identifier("drugname"),
    "Concentration",
    paste0(c(paste0(get_identifier("drugname"), "_"), "Concentration_"), 
      sort(c(2:10, 2:10)))
  )

  HEADERS_LIST[["ordered_2"]] <- c(
    HEADERS_LIST[["normalized_results"]],
    HEADERS_LIST[["averaged_results"]],
    HEADERS_LIST[["metrics_results"]],
    HEADERS_LIST[["raw_data"]],
    HEADERS_LIST[["add_clid"]][-2:-1],
    get_identifier("cellline"),
    get_identifier("drug"),
    paste0(get_identifier("drug"), "_", 2:10),
    HEADERS_LIST[["manifest"]],
    "WellRow",
    "WellColumn"
  )

  HEADERS_LIST
}
