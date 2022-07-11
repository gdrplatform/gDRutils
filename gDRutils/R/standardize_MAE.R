#' Standardize MAE by switching from custom identifiers into gDR-default
#'
#' @param mae a MultiAssayExperiment object with drug-response data generate by gDR pipeline
#'
#' @return mae a MultiAssayExperiment with default gDR identifiers
#' @export
#'
standardize_MAE <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  reset_env_identifiers()
  idfs <- get_env_identifiers()
  experiments <- names(mae)
  for (experiment in experiments) {
    se <- mae[[experiment]]
    idfs_se <- get_SE_identifiers(se)
    matching_default_idfs <- idfs[names(idfs_se)]
    
    matching_default_idfs <- matching_default_idfs[lengths(matching_default_idfs) != 0]
    idfs_se <- idfs_se[names(matching_default_idfs)]
    idf_diff <- idfs_se[unlist(lapply(names(idfs_se),
                                      function(x) !identical(idfs_se[[x]], matching_default_idfs[[x]])))]
    default_idfs_diff <- matching_default_idfs[names(idf_diff)]
    if (length(idf_diff) == 0) next
    mapping_vector_list <- lapply(names(default_idfs_diff), function(x) {
      if (length(default_idfs_diff[[x]]) != length(idf_diff[[x]])) {
        min_val <- min(length(default_idfs_diff[[x]]),  length(idf_diff[[x]]))
        c(name = idf_diff[[x]][seq_len(min_val)],
          value = default_idfs_diff[[x]][seq_len(min_val)])
      } else {
        c(name = idf_diff[[x]],
          value = default_idfs_diff[[x]])
      }
    })
    names(mapping_vector_list) <- names(default_idfs_diff)
    
    mapping_vector <- lapply(mapping_vector_list, "[[", "value")
    names(mapping_vector) <- unlist(lapply(mapping_vector_list, "[[", "name"))
    mapping_vector <- unlist(mapping_vector)
    rowData(se) <- rename_DFrame(rowData(se), mapping_vector)
    colData(se) <- rename_DFrame(colData(se), mapping_vector)
    assayList <- lapply(assays(se), function(x) {
      rename_bumpy(x, mapping_vector)
    })
    assays(se) <- assayList
    idfs_new <- idfs
    idfs_new[names(mapping_vector_list)] <- lapply(mapping_vector_list, "[[", "value")
    se <- set_SE_identifiers(se, idfs_new)
    mae[[experiment]] <- se
  }
  mae
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
