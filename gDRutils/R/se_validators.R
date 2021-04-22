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

