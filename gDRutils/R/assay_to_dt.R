#### AUXILIARY FUNCTIONS ####

#' Transform a SummarizedExperiment assay to a long data.table
#'
#' Transform a SummarizedExperiment assay to a long data.table with a single entry for each row and column combination.
#' @import reshape2
#' @param se \linkS4class{SummarizedExperiment} object with dose-response data.
#' @param assay_name String of name of the assay or index of the assay in the \code{se}.
#' @param merge_metrics Logical indicating whether the metrics should be merged.
#' Defaults to \code{FALSE}.
#' @param include_controls Logical indicating whether the controls should be included.
#' Defaults to \code{FALSE}.
#'
#' @return data.table with dose-response data
#'
#' @export
# TODO: we should rename the function
# - with the verb and more precise
# - maybe 'convert_assay_to_dt' ?
# - alternative naming conventions:
# assay_to_dt => convert_se_to_dt and convert_assay_data_to_dt => convert_assay_to_dt
assay_to_dt <- function(se, assay_name, merge_metrics = FALSE, include_controls = FALSE) {
  # check arguments
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assertTRUE(checkmate::test_count(assay_name) ||
                          checkmate::test_string(assay_name))
  checkmate::assert_flag(merge_metrics)
  checkmate::assert_flag(include_controls)
  
  # define data.table with data from rowData/colData
  rData <- SummarizedExperiment::rowData(se)
  rData$rId <- rownames(rData)
  
  cData <- SummarizedExperiment::colData(se)
  cData$cId <- rownames(cData)
  
  ids <-
    expand.grid(rData$rId,
                cData$cId)
  colnames(ids) <- c("rId", "cId")
  ids[] <- lapply(ids, as.character)
  
  
  annotTbl <-
    merge(merge(ids, rData, by = "rId", all.x = TRUE),
          cData,
          by = "cId",
          all.x = TRUE)
  
  # use method to convert assay data to data_table
  as_dt <-
    convert_assay_data_to_dt(SummarizedExperiment::assay(se, assay_name))
  
  # empty case
  if (nrow(as_dt) == 0) {
    return(as_dt)
  }
 
  as_dt <- if ((assay_name == "Metrics") || 
      (is.numeric(assay_name) && names(SummarizedExperiment::assays(se))[assay_name] == "Metrics")) {
    as_dt <-
      data.table::as.data.table(merge(
        as_dt,
        annotTbl,
        by = c("rId", "cId"),
        all.x = TRUE
      ))
    as_dt$dr_metric <-  rep_len(c("RV", "GR"), nrow(as_dt))
    if (merge_metrics) {
      colnames_RV <- get_header("RV_metrics")
      diff_RV_columns <-
        setdiff(names(colnames_RV), colnames(as_dt))
      if (length(diff_RV_columns) > 0) {
        warning(paste(
          "missing column(s) in SE:",
          paste(diff_RV_columns, collapse = ", ")
        ))
        colnames_RV <-
          colnames_RV[!names(colnames_RV) %in% diff_RV_columns]
      }
      colnames_GR <- get_header("GR_metrics")
      diff_GR_columns <-
        setdiff(names(colnames_GR), colnames(as_dt))
      if (length(diff_GR_columns) > 0) {
        warning(paste(
          "missing column(s) in SE:",
          paste(diff_GR_columns, collapse = ", ")
        ))
        colnames_GR <-
          colnames_GR[!names(colnames_GR) %in% diff_GR_columns]
      }
      
      vars <- c("rId", "cId", names(colnames_RV))
      Df_RV <- as_dt[dr_metric == "RV", ..vars]
      Df_GR <- as_dt[dr_metric == "GR", - "dr_metric"]
      
      data.table::setnames(Df_RV,
                           old = names(colnames_RV),
                           new = unname(colnames_RV))
      data.table::setnames(Df_GR,
                           old = names(colnames_GR),
                           new = unname(colnames_GR))
      merge(Df_RV,
            Df_GR,
            by = c("rId", "cId"),
            all = TRUE)
    } else {
      as_dt
    }
  } else {
    as_dt <- data.table::as.data.table(merge(
      as_dt,
      annotTbl,
      by = c("rId", "cId"),
      all.x = TRUE
    ))
    if (include_controls) {

      as_dt_ctrl <-
        convert_assay_data_to_dt(SummarizedExperiment::assay(se, ifelse(assay_name == "Normalized" ||
          (is.numeric(assay_name) && names(SummarizedExperiment::assays(se))[assay_name] %in% "Normalized"),
            "Controls", "Avg_Controls")))
      as_dt_ctrl <- merge(
        as_dt_ctrl,
        annotTbl,
        by = c("rId", "cId"),
        all.x = TRUE
      )
      as_dt_ctrl <- data.table::as.data.table(as_dt_ctrl)

      as_dt_ctrl$Gnumber <- gDRutils::get_identifier("untreated_tag")[2]
      as_dt_ctrl$DrugName <- gDRutils::get_identifier("untreated_tag")[2]
      as_dt_ctrl$Concentration <- 0
      as_dt_ctrl$std_GRvalue <- NA
      as_dt_ctrl$std_RelativeViability <- NA

      data.table::setnames(as_dt_ctrl,
                           old = c("RefRelativeViability", "RefGRvalue", "RefReadout"),
                           new = c("RelativeViability", "GRvalue", "CorrectedReadout"))
      as_dt_ctrl[, c("UntrtReadout", "DivisionTime", "Day0Readout") := NULL]
      as_dt <- rbind(as_dt, as_dt_ctrl, fill = TRUE)
    }
    
    as_dt
  }
}

#' Convert assay data into data.table
#' @param object An object comprising assay in SummarizedExperiment
#' @return data.table with assay data
convert_assay_data_to_dt <- function(object) {
  UseMethod("convert_assay_data_to_dt")
}

#' @rdname convert_assay_data_to_dt
convert_assay_data_to_dt.default <- function(object) {
  stop(
    paste(
      "Unable to find an inherited method for function 'convert_assay_data_to_dt' for signature",
      shQuote(class(object))
    )
  )
}

convert_assay_data_to_dt.BumpyMatrix <-
  function(object) {
      BumpyMatrix::unsplitAsDataFrame(object,
                                      row.field = "rId",
                                      column.field = "cId")
  }

convert_assay_data_to_dt.matrix <- function(object) {
  
  # we expect matrix object to be the list of DFrame(s) or NULL(s)
  checkmate::assertTRUE(all(lapply(object, class) %in% c("DFrame", "NULL")))
  
  asL <-
    lapply(seq_len(ncol(object)), function(x) {
      myL <- object[, x]
      # if only one row (nrow==1), the name of the row is not kept which results in a bug
      names(myL) <- rownames(object) # this line is not affecting results if now>1
      
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
        df <-
          do.call(plyr::rbind.fill,
                  lapply(myL, data.table::as.data.table))
      } else {
        df <-
          data.table::rbindlist(lapply(myL, data.table::as.data.table), fill = TRUE)
      }
      if (nrow(df) == 0)
        return()
      df$rId <- rCol
      df$cId <- colnames(object)[x]
      df
    })
  as_df <- data.table::rbindlist(asL)
}

#' .get_untreated_conditions
#'
#' Get untreated conditions
#'
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with untreated conditions
#'
.get_untreated_conditions <-
  function(drug_data) {
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    .untreated_tag_patterns <- vapply(get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
    .untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")
    drugnames <- tolower(as.data.frame(drug_data)[, get_identifier("drugname")])
    drug_data[grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }

#' .get_treated_conditions
#'
#' Get treated conditions
#'
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with treated conditions
#'
.get_treated_conditions <-
  function(drug_data) {
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    .untreated_tag_patterns <- vapply(get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
    .untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")
    drugnames <- tolower(data.table::as.data.table(drug_data)[, get_identifier("drugname")])
    drug_data[!grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }
