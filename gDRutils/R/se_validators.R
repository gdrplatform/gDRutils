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
validate_se_assay_name <- function(se, name) {
  if (!name %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("'%s' is not on of the available assays: '%s'",
      name, paste0(SummarizedExperiment::assayNames(se), collapse = ", ")))
  }
  invisible(NULL)
}



#' Validate SummarizedExperiment object
#' 
#'
#' @param se SummarizedExperiment object 
#' produced by the gDR pipeline
#'
#' @return
#' @export
#'
#' @examples
validate_SE <- function(SE) {
  checkmate::assert(
    checkmate::check_class(SE, "SummarizedExperiment"),
    checkmate::check_subset(SummarizedExperiment::assayNames(SE), c("RawTreated", "Controls", 
                                                                       "Normalized", "RefGRvalue", 
                                                                       "RefRelativeViability", "DivisionTime", 
                                                                       "Averaged", "Metrics")),
    checkmate::check_set_equal(names(S4Vectors::metadata(SE)),
                               c("experiment_metadata", "df_",
                                 "Keys", "df_raw_data",
                                 "fit_parameters", 
                                 #"drug_combinations",
                                 ".internals")),
    checkmate::check_names(names(gDRutils::convert_se_assay_to_dt(SE, "Metrics")),
                           must.include = c("normalization_type", "fit_source")),
    checkmate::check_true(identical(dimnames(SE), dimnames(SummarizedExperiment::assay(SE, "Normalized")))),
    combine = "and")
  coldata <- colData(SE)
  rowdata <- rowData(SE)
  checkmate::assert(
    checkmate::check_true(expect_equal(gsub("_.*", "", rownames(SE)), rowdata$Gnumber)),
    checkmate::check_true(expect_equal(gsub("_.*", "", colnames(SE)), coldata$clid)),
    checkmate::check_true(nrow(coldata) == nrow(unique(coldata))),
    checkmate::check_true(nrow(rowdata) == nrow(unique(rowdata))),
    combine = "and")
}

