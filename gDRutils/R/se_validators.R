#' Check whether or not an assay exists in a SummarizedExperiment object.
#'
#' Check for the presence of an assay in a SummarizedExperiment object.
#'
#' @param se A \linkS4class{SummarizedExperiment} object.
#' @param name String of name of the assay to validate.
#'
#' @return \code{NULL} invisibly if the assay name is valid.
#' Throws an error if the assay is not valid.
#'
#' @export
#'
is_valid_se_assay_name <- function(se, name) {
  if (!name %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("'%s' is not on of the available assays: '%s'",
      name, paste0(SummarizedExperiment::assayNames(se), collapse = ", ")))
  }
  invisible(NULL)
}



#' Validate SummarizedExperiment object
#' 
#'
#' @param SE SummarizedExperiment object 
#' produced by the gDR pipeline
#'
#' @return
#' @export
#'
#' @examples
validate_SE <- function(SE) {
  checkmate::assert(
    checkmate::check_class(SE, "SummarizedExperiment"),
    checkmate::check_set_equal(SummarizedExperiment::assayNames(SE), c("RawTreated", "Controls", 
                                                                       "Normalized", "RefGRvalue", 
                                                                       "RefRelativeViability", "DivisionTime", 
                                                                       "Averaged", "Metrics")),
    checkmate::check_set_equal(names(S4Vectors::metadata(SE)),
                               c("experiment_metadata", "df_",
                                 "Keys", "df_raw_data",
                                 "fit_parameters", "drug_combinations")),
    combine = "and"
  )
  unsplittedBM <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(
      SE, "Metrics"))
  checkmate::assert(
    checkmate::assert("Metrics assay must have only 2 entries",
                      all(table(unsplittedBM$row, unsplittedBM$column) == 2)),
    checkmate::check_names(names(gDRutils::convert_se_assay_to_dt(SE, "Metrics")),
                           must.include = c("metric_type", "fit_source")),
    
    combine = "and")
}

