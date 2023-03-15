## The following function utilizes the get_env_identifiers() function which can be
## changed at run time, which is why it needs to be wrapped in a function. 

#' @keywords internal
.getHeadersList <- function() {
  HEADERS_LIST <- list(
    manifest = get_env_identifiers(c("barcode", "template", "duration"), simplify = FALSE),

    raw_data = c(
      "ReadoutValue",
      "BackgroundValue",
      "UntrtReadout",
      "Day0Readout",
      get_env_identifiers("masked_tag", simplify = TRUE)
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
      "ec50",
      "x_inf",
      "x_0",
      "h",
      "r2",
      "x_sd_avg",
      "fit_type"
    ),

    # corresponds to the field "celllinename", "primarytissue", "doublingtime" from gneDB CLIDs
    add_clid = get_env_identifiers(c("cellline_name", "cellline_tissue", "cellline_ref_div_time", "cellline_subtype"), simplify = FALSE)
  )

  metrics_names <- rbind(
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
  colnames(metrics_names) <- HEADERS_LIST[["response_metrics"]]
  HEADERS_LIST[["metrics_names"]] <- metrics_names

  HEADERS_LIST[["metrics_results"]] <- c(
    "maxlog10Concentration",
    "N_conc",
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
    HEADERS_LIST[["raw_data"]],
    HEADERS_LIST[["normalized_results"]],
    HEADERS_LIST[["averaged_results"]],
    HEADERS_LIST[["metrics_results"]],
    get_env_identifiers("well_position", simplify = TRUE)
  )

  HEADERS_LIST[["ordered_1"]] <- c(
    get_env_identifiers("cellline_name", simplify = TRUE),
    get_env_identifiers("cellline_tissue", simplify = TRUE),
    get_env_identifiers("duration", simplify = TRUE),
    get_env_identifiers("drug_name", simplify = TRUE),
    "Concentration",
    paste0(c(paste0(get_env_identifiers("drug_name", simplify = TRUE), "_"), "Concentration_"), 
      rep(2:10, each = 2))
  )

  HEADERS_LIST[["ordered_2"]] <- c(
    HEADERS_LIST[["normalized_results"]],
    HEADERS_LIST[["averaged_results"]],
    HEADERS_LIST[["metrics_results"]],
    HEADERS_LIST[["raw_data"]],
    get_env_identifiers("cellline_ref_div_time", simplify = TRUE),
    get_env_identifiers("cellline", simplify = TRUE),
    get_env_identifiers("drug", simplify = TRUE),
    paste0(get_env_identifiers("drug", simplify = TRUE), "_", 2:10),
    HEADERS_LIST[["manifest"]],
    get_env_identifiers("well_position", simplify = TRUE)
  )

  HEADERS_LIST
}
