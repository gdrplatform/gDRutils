#' Transform a SummarizedExperiment assay to a long data.table
#'
#' Transform a SummarizedExperiment assay to a long data.table with a single entry for each row and column combination.
#' @import reshape2
#' @param se SummarizedExperiment object with dose-response data.
#' @param assay_name String of name of the assay in the /code{se}.
#' @param merge_metrics Logical indicating whether the metrics should be merged.
#' Defaults to \code{FALSE}.
#'
#' @return data.table with dose-response data
#'
#' @export
assay_to_dt <- function(se, assay_name, merge_metrics = FALSE) {
  # check arguments
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assertTRUE(checkmate::test_count(assay_name) || checkmate::test_string(assay_name))
  checkmate::assert_logical(merge_metrics)

  # define data.table with data from rowData/colData
  ids <- expand.grid(rownames(SummarizedExperiment::rowData(se)), rownames(SummarizedExperiment::colData(se)))
  colnames(ids) <- c("rId", "cId")
  ids[] <- lapply(ids, as.character)
  rData <- SummarizedExperiment::rowData(se)
  rData$rId <- rownames(rData)
  cData <- SummarizedExperiment::colData(se)
  cData$cId <- rownames(cData)
  annotTbl <- merge(ids, rData, by = "rId", all.x = TRUE)
  annotTbl <- merge(annotTbl, cData, by = "cId", all.x = TRUE)
  
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
      df <- do.call(plyr::rbind.fill, lapply(myL, data.table::as.data.table))
    } else {
      df <- data.table::as.data.table(do.call(rbind, myL))
    }
    if(nrow(df)==0) return()
    df$rId <- rCol
    df$cId <- rownames(SummarizedExperiment::colData(se))[x]
    full.df <- merge(df, annotTbl, by = c("rId", "cId"), all.x = TRUE)
  })
  asDf <- data.table::rbindlist(asL)
  if (assay_name == "Metrics") {
    asDf <- asDf[, dr_metric := rep_len(c("RV", "GR"), .N)]
    if (merge_metrics) {
      
      colnames_RV <- gDRutils::get_header("RV_metrics")
      diff_RV_columns <- setdiff(names(colnames_RV), colnames(asDf))
      if (length(diff_RV_columns) > 0) {
        warning(paste("missing column(s) in SE:", paste(diff_RV_columns, collapse = ", ")))
        colnames_RV <- colnames_RV[!names(colnames_RV) %in% diff_RV_columns]
      }
      colnames_GR <- gDRutils::get_header("GR_metrics")
      diff_GR_columns <- setdiff(names(colnames_GR), colnames(asDf))
      if (length(diff_GR_columns) > 0) {
        warning(paste("missing column(s) in SE:", paste(diff_GR_columns, collapse = ", ")))
        colnames_GR <- colnames_GR[!names(colnames_GR) %in% diff_GR_columns]
      }
      
      Df_RV <- subset(asDf, dr_metric == "RV", select = c("rId", "cId", names(colnames_RV)))
      Df_GR <- subset(asDf, dr_metric == "GR", select = -c(dr_metric))
      
      data.table::setnames(Df_RV,
                           old = names(colnames_RV),
                           new = unname(colnames_RV))
      data.table::setnames(Df_GR,
                           old = names(colnames_GR),
                           new = unname(colnames_GR))
      asDf <- merge(Df_RV, Df_GR, by = c("rId", "cId"), all = TRUE)
    }
  }
  asDf
}
