#' Validate dimnames
#'
#' Assure that dimnames of two objects are the same
#'
#' @param obj first object with dimnames to compare
#' @param obj2 second object with dimnames to compare
#' @param skip_empty a logical indicating whether to skip comparison
#'        if a given dimname is empty in the case of both objects
#' @return \code{NULL}
#'
validate_dimnames <- function(obj, obj2, skip_empty = TRUE) {

  dn1 <- dimnames(obj)
  dn2 <- dimnames(obj2)

  checkmate::assert_true(identical(length(dn1), length(dn2)))

  for (idx in seq_along(dn1)) {
    found_non_empty <- sum(length(dn1[[idx]]), length(dn2[[idx]])) > 0
    if (!skip_empty || (skip_empty && found_non_empty)) {
      checkmate::assert_true(identical(dn1[idx], dn2[idx]))
    }
  }
}

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
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small") 
#' se <- mae[[1]]
#' validate_se_assay_name(se, "RawTreated")
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
#' @param expect_single_agent a logical indicating if the function
#' should check whether the SummarizedExperiment is single-agent data
#'
#' @return \code{NULL} invisibly if the SummarizedExperiment is valid.
#' Throws an error if the SummarizedExperiment is not valid.
#' @export
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' validate_SE(se)
#'
validate_SE <- function(se,
                        expect_single_agent = FALSE) {
  # Validate the SE structure, assays and metadata, as well as dimnames of assays
  checkmate::assert_class(se, "SummarizedExperiment")
  exp_assay_names <- c("Normalized", "Averaged")
  if (expect_single_agent) {
    exp_assay_names <- c(exp_assay_names, "Metrics")
  }
  checkmate::assert_subset(exp_assay_names, assayNames(se))
  exp_metadata_names <- c("experiment_metadata", "Keys")
  checkmate::assert_true(all(exp_metadata_names %in% names(S4Vectors::metadata(se))))
  validate_dimnames(se, assay(se, "Normalized"))
  validate_dimnames(se, assay(se, "Averaged"))
  if (expect_single_agent) {
    validate_dimnames(se, assay(se, "Metrics"))
    checkmate::assert_true(all(c("normalization_type", "fit_source") %in%
                                 names(convert_se_assay_to_dt(se, "Metrics"))))
  }
  coldata <- colData(se)
  rowdata <- rowData(se)
  # Validate the correctness of rowData and colData
  # building gsub expression by finding location between _x_ for drug identifier
  #nolint start
  pattern <- paste0("^(?:[^_]+_){",
                    as.character(which(names(rowdata) == get_env_identifiers("drug")) - 1),
                    "}([^_]+).*")
  #nolint end
  checkmate::assert_true(nrow(coldata) == nrow(unique(coldata)))
  checkmate::assert_true(nrow(rowdata) == nrow(unique(rowdata)))

  # Validate non-empty values in rowData and colData
  checkmate::assert_false(any(stats::na.omit(unlist(coldata)) == ""))
  checkmate::assert_false(any(stats::na.omit(unlist(rowdata)) == ""))

  # Validate the correctness of single-agent data
  drug_name2 <- get_env_identifiers("drug_name2")
  concentration2 <- get_env_identifiers("concentration2")
  vars_cotreatment <- intersect(c(drug_name2, concentration2), names(rowdata))
  if (expect_single_agent && length(vars_cotreatment) > 0) {
    if (drug_name2 %in% names(rowdata)) {
      checkmate::assert_subset(rowdata[[drug_name2]], get_SE_identifiers(se, "untreated_tag", simplify = TRUE))
    }
    if (concentration2 %in% names(rowdata)) {
      checkmate::assert_subset(rowdata[[concentration2]], 0)
    }
  }
  invisible(NULL)
}

#' Validate MultiAssayExperiment object
#'
#' Function validates correctness of SE included in MAE by checking multiple cases:
#' - detection of duplicated rowData/colData,
#' - incompatibility of rownames/colnames,
#' - occurrence of necessary assays,
#' - detection of mismatch of CLIDs inside colData and colnames (different order),
#' - correctness of metadata names.
#'
#' @param mae MultiAssayExperiment object
#' produced by the gDR pipeline
#'
#' @return \code{NULL} invisibly if the MultiAssayExperiment is valid.
#' Throws an error if the MultiAssayExperiment is not valid.
#' @export
#'
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small") 
#' validate_MAE(mae)
#' 
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
validate_MAE <- function(mae) {
  # Validate the SE structure, assays and metadata, as well as dimnames of assays
  checkmate::assert_class(mae, "MultiAssayExperiment")
  experiments <- names(mae)
  checkmate::assert_subset(experiments, c("single-agent", "co-dilution",
                                           "matrix", "cotreatment", "other"))
  for (experiment in experiments) {
    if (experiment %in% get_experiment_groups("single-agent")[["single-agent"]]) {
      expect_single_agent <- TRUE
    } else {
      expect_single_agent <- FALSE
    }
    validate_SE(mae[[experiment]], expect_single_agent = expect_single_agent)
  }
  invisible(NULL)
}

