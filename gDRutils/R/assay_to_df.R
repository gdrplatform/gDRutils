#' Transform a SummarizedExperiment assay to a long data.frame.
#'
#' Transform a SummarizedExperiment assay to a long data.frame with a single entry for each row and column combination.
#' @import reshape2
#' @param se SummarizedExperiment object with dose-response data.
#' @param assay_name String of name of the assay in the /code{se}.
#' @param merge_metrics Logical indicating whether the metrics should be merged.
#' Defaults to \code{FALSE}.
#'
#' @return data.frame with dose-response data
#'
#' @export
assay_to_df <- function(se, assay_name, merge_metrics = FALSE) {
  # check arguments
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
    if(nrow(df)==0) return()
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
      if (length(diff_RV_columns) > 0) {
        warning(paste("missing column(s) in SE:", paste(diff_RV_columns, collapse = ", ")))
        colnames_IC <- colnames_IC[!names(colnames_IC) %in% diff_RV_columns]
      }
      colnames_GR <- gDRutils::get_header("GR_metrics")
      diff_GR_columns <- setdiff(names(colnames_GR), colnames(asDf))
      if (length(diff_GR_columns) > 0) {
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