
#' @keywords internal
.extract_matching_identifiers <- function(default, se_identifiers) {
  matching_default_idfs <- default[names(se_identifiers)]
  # Drop non-matching identifiers
  matching_default_idfs <- matching_default_idfs[lengths(matching_default_idfs) != 0]
  idfs_se <- se_identifiers[names(matching_default_idfs)]
  list(default = matching_default_idfs,
       se = idfs_se)
}

#' @keywords internal
.extract_diff_identifiers <- function(default,
                                      se_identifiers) {
  
  diff_names <- names(which(vapply(names(se_identifiers),
                                   function(x) !identical(se_identifiers[[x]], default[[x]]),
                                   FUN.VALUE = logical(1))))
  list(default = default[diff_names],
       se = se_identifiers[diff_names])
}

#' @keywords internal
.replace_untreated_tag <- function(row_data,
                                   default,
                                   se_identifiers) {
  untreated_tag <- data.frame(x = default[["untreated_tag"]],
                              y = se_identifiers[["untreated_tag"]])
  untreated_tag_mapping <- untreated_tag$x
  names(untreated_tag_mapping) <- untreated_tag$y
  S4Vectors::DataFrame(lapply(row_data, function(x) stringr::str_replace_all(x, untreated_tag_mapping)),
                       row.names = rownames(row_data))
}

#' Rename DFrame
#'
#' @param df a DFrame object
#' @param mapping_vector a named vector for mapping old and new values.
#' The names of the character vector indicate the source names, and the
#' corresponding values the destination names. 
#'
#' @return a renamed DFrame object
#' @export
#'
rename_DFrame <- function(df, mapping_vector) {
  checkmate::assert_class(df, "DFrame")
  names(df)[names(df) %in% names(mapping_vector)] <- 
    mapping_vector[names(df)[names(df) %in% names(mapping_vector)]]
  df
}

#' Rename BumpyMatrix
#'
#' @param bumpy a BumpyMatrix object
#' @param mapping_vector a named vector for mapping old and new values.
#' The names of the character vector indicate the source names, and the
#' corresponding values the destination names. 
#'
#' @return a renamed BumpyMatrix object
#' @export
rename_bumpy <- function(bumpy, mapping_vector) {
  checkmate::assert_class(bumpy, "BumpyMatrix")
  bumpy_cols <- BumpyMatrix::commonColnames(bumpy)
  mapping <- mapping_vector[bumpy_cols]
  mapping[is.na(mapping)] <- bumpy_cols[is.na(mapping)]
  BumpyMatrix::commonColnames(bumpy) <- mapping
  bumpy
}


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
#' @param rd DataFrame with colData
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
