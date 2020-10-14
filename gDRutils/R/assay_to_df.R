#' Transform a SummarizedExperiment assay to a long data.frame.
#'
#' Transform a SummarizedExperiment assay to a long data.frame with a single entry for each row and column combination.
#'
#' @param se SummarizedExperiment object with dose-response data.
#' @param assay_name String of name of the assay in the /code{se}.
#' @param merge_metrics Logical indicating whether the metrics should be merged.
#' Defaults to \code{FALSE}.
#'
#' @return data.frame with dose-response data
#'
#' @importFrom checkmate assert_class assertTRUE assert_logical test_count test_string
#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom dplyr left_join select full_join
#' @importFrom plyr rbind.fill
#' @importFrom data.table set.names
#' @export
assay_to_df <- function(se, assay_name, merge_metrics = FALSE) {
  # check arguments
  assert_class(se, "SummarizedExperiment")
  assertTRUE(test_count(assay_name) || test_string(assay_name))
  assert_logical(merge_metrics)

  # define data.frame with data from rowData and colData
  rData <- data.frame(rowData(se), stringsAsFactors = FALSE)
  rData$rId <- rownames(rData)

  cData <- data.frame(colData(se), stringsAsFactors = FALSE)
  cData$cId <- rownames(cData)

  ids <- expand.grid(rData$rId, cData$cId)
  ids[] <- lapply(ids, as.character)
  colnames(ids) <- c("rId", "cId")

  annotTbl <- left_join(ids, rData, by = "rId")
  annotTbl <- left_join(annotTbl, cData, by = "cId")

  # merge assay data with data from colData/rowData
  SE_assay <- assay(se, assay_name)
  asL <- lapply(seq_len(ncol(SE_assay)), function(x) {
    myL <- SE_assay[, x, drop = FALSE]

    ## filter out empty DataFrames for given drug/cell_line combination
    ## example: testdata 7 - assay(seL[[7]], "Normalized")[5:8,1]
    myL <- myL[vapply(myL, nrow, integer(1)) > 0L]

    myV <- vapply(myL, nrow, integer(1))
    rCol <- rep(names(myV), as.numeric(myV))

    # Fill NAs for DataFrames with different number of columns
    if (length(unique(sapply(myL, ncol))) > 1L) {
      df <- do.call(rbind.fill, lapply(myL, data.frame))
    } else {
      df <- data.frame(do.call(rbind, myL))
    }

    if (nrow(df) == 0L) {
        return()
    }

    df$rId <- rCol
    df$cId <- rownames(colData(se))[x]
    full.df <- left_join(df, annotTbl, by = c("rId", "cId"))
  })

  asDf <- data.frame(do.call(rbind, asL))

  if (assay_name == "Metrics") {
    asDf$dr_metric <- c("IC", "GR")
    if (merge_metrics) {
      colnames_IC <- gDRutils::get_header("RV_metrics")
      diff_RV_columns <- setdiff(names(colnames_IC), colnames(asDf))
      if (length(diff_RV_columns) > 0L) {
        warning(paste("missing column(s) in SE:", paste(diff_RV_columns, collapse = ", ")))
        colnames_IC <- colnames_IC[!names(colnames_IC) %in% diff_RV_columns]
      }

      colnames_GR <- gDRutils::get_header("GR_metrics")
      diff_GR_columns <- setdiff(names(colnames_GR), colnames(asDf))
      if (length(diff_GR_columns) > 0L) {
        warning(paste("missing column(s) in SE:", paste(diff_GR_columns, collapse = ", ")))
        colnames_GR <- colnames_GR[!names(colnames_GR) %in% diff_GR_columns]
      }

      Df_IC <- subset(asDf, dr_metric == "IC", select = c("rId", "cId", names(colnames_IC)))
      Df_GR <- subset(asDf, dr_metric == "GR") %>% select(-dr_metric)

      data.table::setnames(Df_IC,
                           old = names(colnames_IC),
                           new = unname(colnames_IC))

      data.table::setnames(Df_GR,
                           old = names(colnames_GR),
                           new = unname(colnames_GR))
      asDf <- full_join(Df_IC, Df_GR, by = c("rId", "cId"))
    }
  }

  asDf
}
