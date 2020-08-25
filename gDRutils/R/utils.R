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

    masked_tag = 'masked',
    # flag for masked wells

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
      get_identifier("masked_tag"),
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

#' assay_to_df
#'
#' Convert SE assay to data.frame
#'
#' @param se  SummarizedExperiment object with dose-response data
#' @param assay_name string name of the assay
#' @param merge_metrics a logical indicating whether the metrics should be merged
#'
#' @return data.frame with dose-reponse data
#'
#' @export
assay_to_df <- function(se, assay_name, merge_metrics = FALSE) {
  stopifnot(any("SummarizedExperiment" %in% class(se)))
  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assertTRUE(checkmate::test_count(assay_name) || checkmate::test_string(assay_name))
  checkmate::assert_logical(merge_metrics)

  # define data.frame with data from rowData/colData
  ids <- expand.grid(rownames(SummarizedExperiment::rowData(se)), rownames(SummarizedExperiment::colData(se)))
  colnames(ids) <- c("rId", "cId")
  ids[] <- lapply(ids, as.character)
  rData <- data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)
  rData$rId <- rownames(rData)
  cData <- data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  cData$cId <- rownames(cData)
  annotTbl <-
    dplyr::left_join(ids, rData, by = "rId")
  annotTbl <-
    dplyr::left_join(annotTbl, cData, by = "cId")

  #merge assay data with data from colData/rowData
  SE_assay = SummarizedExperiment::assay(se, assay_name)
  asL <- lapply(1:nrow(SummarizedExperiment::colData(se)), function(x) {
    myL <- SE_assay[, x]
    # if only one row (nrow==1), the name of the row is not kept which results in a bug
    names(myL) = rownames(SE_assay) # this line is not affecting results if now>1

    # in some datasets there might be no data for given drug/cell_line combination
    # under such circumstances DataFrame will be empty
    # and should be filtered out
    # example: testdata 7 - SummarizedExperiment::assay(seL[[7]],"Normalized")[5:8,1]
    myL <- myL[vapply(myL, nrow, integer(1)) > 0]

    myV <- vapply(myL, nrow, integer(1))
    rCol <- rep(names(myV), as.numeric(myV))
    # there might be DataFrames with different number of columns
    # let's fill with NAs where necessary
    if (length(unique(sapply(myL, ncol))) > 1) {
      df <- do.call(plyr::rbind.fill, lapply(myL, data.frame))
    } else {
      df <- data.frame(do.call(rbind, myL))
    }

    df$rId <- rCol
    df$cId <- rownames(SummarizedExperiment::colData(se))[x]
    full.df <- dplyr::left_join(df, annotTbl, by = c("rId", "cId"))
  })
  asDf <- data.frame(do.call(rbind, asL))
  if (assay_name == "Metrics") {
    asDf$dr_metric <- c("IC", "GR")
    if (merge_metrics) {
      colnames_IC <- gDRutils::get_header("RV_metrics")
      diff_RV_columns <- setdiff(names(colnames_IC), colnames(asDf))
      if (!is.null(diff_RV_columns)) {
        warning(paste("missing column(s) in SE:", paste(diff_RV_columns, collapse = ", ")))
        colnames_IC <- colnames_IC[!names(colnames_IC) %in% diff_RV_columns]
      }
      colnames_GR <- gDRutils::get_header("GR_metrics")
      diff_GR_columns <- setdiff(names(colnames_GR), colnames(asDf))
      if (!is.null(diff_GR_columns)) {
        warning(paste("missing column(s) in SE:", paste(diff_GR_columns, collapse = ", ")))
        colnames_GR <- colnames_GR[!names(colnames_GR) %in% diff_GR_columns]
      }

      Df_IC <- subset(asDf, dr_metric == "IC", select = c("rId", "cId", names(colnames_IC)))
      Df_GR <- subset(asDf, dr_metric == "GR") %>% dplyr::select(-dr_metric)

      data.table::setnames(Df_IC,
                           old = names(colnames_IC),
                           new = unname(colnames_IC))
      data.table::setnames(Df_GR,
                           old = names(colnames_GR),
                           new = unname(colnames_GR))
      asDf <- dplyr::full_join(Df_IC, Df_GR, by = c("rId", "cId"))
    }
  }
  asDf
}
