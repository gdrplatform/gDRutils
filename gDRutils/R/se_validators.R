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
  if (!name %in% assayNames(se)) {
    stop(sprintf("'%s' is not on of the available assays: '%s'",
      name, paste0(assayNames(se), collapse = ", ")))
  }
  invisible(NULL)
}



#' Validate SummarizedExperiment object
#' 
#' Function validates correctness of SE by checking multiple cases:
#' - detection of duplicated rowData/colData,
#' - incompatibility of rownames/colnames,
#' - occurrence of necessary assays,
#' - detection of mismatch of CLIDs inside colData and colnames (different order),
#' - correctness of metadata names.
#'
#' @param se SummarizedExperiment object 
#' produced by the gDR pipeline
#'
#' @return \code{NULL} invisibly if the SummarizedExperiment is valid.
#' Throws an error if the SummarizedExperiment is not valid.
#' @export
#'
validate_SE <- function(se) {
  checkmate::assert(
    checkmate::check_class(se, "SummarizedExperiment"),
    checkmate::check_subset(assayNames(se), c("RawTreated", "Controls", 
                                                                       "Normalized", "RefGRvalue", 
                                                                       "RefRelativeViability", "DivisionTime", 
                                                                       "Averaged", "Metrics")),
    checkmate::check_set_equal(names(S4Vectors::metadata(se)),
                               c("experiment_metadata", "df_",
                                 "Keys", "df_raw_data",
                                 "fit_parameters", 
                                 #"drug_combinations",
                                 ".internals")),
    checkmate::check_names(names(convert_se_assay_to_dt(se, "Metrics")),
                           must.include = c("normalization_type", "fit_source")),
    checkmate::check_true(identical(dimnames(se), dimnames(assay(se, "Normalized")))),
    checkmate::check_true(identical(dimnames(se), dimnames(assay(se, "Averaged")))),
    checkmate::check_true(identical(dimnames(se), dimnames(assay(se, "Metrics")))),
    combine = "and")
  coldata <- colData(se)
  rowdata <- rowData(se)
  checkmate::assert(
    checkmate::check_true(expect_equal(gsub("_.*", "", rownames(se)), rowdata$Gnumber)),
    checkmate::check_true(expect_equal(gsub("_.*", "", colnames(se)), coldata$clid)),
    checkmate::check_true(nrow(coldata) == nrow(unique(coldata))),
    checkmate::check_true(nrow(rowdata) == nrow(unique(rowdata))),
    combine = "and")
  invisible(NULL)
}

