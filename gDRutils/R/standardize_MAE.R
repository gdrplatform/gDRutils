#' get optional colData fields
#' 
#' @param se a SummarizedExperiment object with drug-response data generate by gDR pipeline
#'
#' @return a charvec containing the names of the optional identifiers in the SE colData
#'
get_optional_coldata_fields <- function(se) {
  idfs <- get_SE_identifiers(se)
  
  as.character(idfs["cellline_tissue"])
}

#' get optional rowData fields
#' 
#' @param se a SummarizedExperiment object with drug-response data generate by gDR pipeline
#'
#' @return a charvec containing the names of the optional identifiers in the SE rowData
#'
get_optional_rowdata_fields <- function(se) {
  idfs <- get_SE_identifiers(se)
  
  out <- c(idfs["drug_moa"])
  
  if (!is.null(idfs["drug2"])) {
    out <- c(out, idfs["drug_moa2"])
  }
  if (!is.null(idfs["drug3"])) {
    out <- c(out, idfs["drug_moa3"])
  }
  
  as.character(out)
}

#' refine colData
#' 
#' current improvements done on the colData as a standardization step:
#' - set default value for optional colData fields
#' 
#' @param cd DataFrame with colData
#' @param se a SummarizedExperiment object with drug-response data generate by gDR pipeline
#' @param default_v string with default value for optional columns in colData
#'
#' @return refined colData
#' @export
#'
refine_coldata <- function(cd, se, default_v = "Undefined") {
  
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_class(cd, "DataFrame")
  checkmate::assert_string(default_v)
  
  undef_fields <- setdiff(get_optional_coldata_fields(se), colnames(cd))
  
  if (length(undef_fields)) {
    cd[, undef_fields] <- default_v
  }
  cd
}

#' refine rowData
#' 
#' current improvements done on the rowData as a standardization step:
#' - set default value for optional rowData fields
#'
#' @param rd DataFrame with rowData
#' @param se a SummarizedExperiment object with drug-response data generate by gDR pipeline
#' @param default_v string with default value for optional columns in rowData
#' 
#' @return refined rowData
#' @export
#'

refine_rowdata <- function(rd, se, default_v = "Undefined") {
  
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_class(rd, "DataFrame")
  checkmate::assert_string(default_v)
  
  undef_fields <- setdiff(get_optional_rowdata_fields(se), colnames(rd))
  
  if (length(undef_fields)) {
    rd[, undef_fields] <- default_v
  }
  rd
}
