#' @export
get_identifier <- function(x = NULL) {
  identifiersList <- list(
    duration = "Duration",

    cellline = "clid",

    drug = "Gnumber",
    drugname = "DrugName",
    # corresponds to the fieLd  'gcsi_drug_name' from gCellGenomics::getDrugs()

    untreated_tag = c("untreated", "vehicle"),
    # flag to identify control treatments

    WellPosition = c("WellRow", "WellColumn")
  )
  if (!is.null(x) &&
      x %in% names(identifiersList))
    return(identifiersList[[x]])
  else
    return(identifiersList)
}

#######-------------------------------------------------------
# these should not be changed and are protected field names
#' @export
get_header <- function(x = NULL) {
  headersList <- list(
    manifest = c("Barcode", "Template", get_identifier("duration")),
    raw_data = c(
      "ReadoutValue",
      "BackgroundValue",
      "UntrtReadout",
      "Day0Readout"
    ),
    normalized_results = c(
      "CorrectedReadout",
      "GRvalue",
      "RelativeViability",
      "DivisionTime",
      "RefGRvalue",
      "RefRelativeViability"
    ),
    averaged_results = c("std_GRvalue", "std_RelativeViability"),
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
    add_clid = c("CellLineName", "Tissue", "ReferenceDivisionTime")
    # corresponds to the fieLd  "celllinename", "primarytissue", "doublingtime" from gneDB CLIDs
  )
  headersList[["RV_metrics"]] <-
    array(
      c(
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
      ),
      dimnames = headersList["response_metrics"]
    )
  headersList[["GR_metrics"]] <-
    array(
      c(
        "mean_GR",
        "GR_AOC",
        "RG_AOC_range",
        "GR50",
        "GR_max",
        "GEC50",
        "GR_inf",
        "GR_0",
        "h_GR",
        "GR_r2",
        "GR_sd_avg",
        "flat_fit_GR"
      ),
      dimnames = headersList["response_metrics"]
    )
  headersList[["metrics_results"]] <-
    c("maxlog10Concentration",
      "N_conc",
      headersList[["response_metrics"]],
      headersList[["RV_metrics"]],
      headersList[["GR_metrics"]])
  headersList[["controlled"]] <- c(
    get_identifier("cellline"),
    headersList[["manifest"]],
    get_identifier("drug"),
    "Concentration",
    paste0(get_identifier("drug"), "_", 2:10),
    paste0("Concentration_", 2:10)
  )
  headersList[["reserved"]] <-
    c(
      headersList[["add_clid"]],
      get_identifier("drugname"),
      paste0(get_identifier("drugname"), "_", 2:10),
      headersList[["raw_data"]],
      headersList[["normalized_results"]],
      headersList[["averaged_results"]],
      headersList[["metrics_results"]],
      "WellRow",
      "WellColumn"
    )

  headersList[["ordered_1"]] <- c(
    headersList[["add_clid"]][1:2],
    get_identifier("duration"),
    get_identifier("drugname"),
    "Concentration",
    paste0(c(
      paste0(get_identifier("drugname"), "_"), "Concentration_"
    ),
    sort(c(2:10, 2:10)))
  )
  headersList[["ordered_2"]] <- c(
    headersList[["normalized_results"]],
    headersList[["averaged_results"]],
    headersList[["metrics_results"]],
    headersList[["raw_data"]],
    headersList[["add_clid"]][-2:-1],
    get_identifier("cellline"),
    get_identifier("drug"),
    paste0(get_identifier("drug"), "_", 2:10),
    headersList[["manifest"]],
    "WellRow",
    "WellColumn"
  )


  if (!is.null(x)) { 
    stopifnot(x %in% names(headersList))
    return(headersList[[x]])
  } else { 
    return(headersList)
  }
}
