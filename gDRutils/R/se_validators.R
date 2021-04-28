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
validate_SE <- function(se) {
  checkmate::assert(
    checkmate::check_class(se, "SummarizedExperiment"),
    checkmate::check_subset(SummarizedExperiment::assayNames(se), c("RawTreated", "Controls", 
                                                                       "Normalized", "RefGRvalue", 
                                                                       "RefRelativeViability", "DivisionTime", 
                                                                       "Averaged", "Metrics")),
    checkmate::check_set_equal(names(S4Vectors::metadata(se)),
                               c("experiment_metadata", "df_",
                                 "Keys", "df_raw_data",
                                 "fit_parameters", 
                                 #"drug_combinations",
                                 ".internals")),
    checkmate::check_names(names(gDRutils::convert_se_assay_to_dt(se, "Metrics")),
                           must.include = c("normalization_type", "fit_source")),
    checkmate::check_true(identical(dimnames(se), dimnames(SummarizedExperiment::assay(se, "Normalized")))),
    combine = "and")
}

