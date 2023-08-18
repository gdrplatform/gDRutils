#' Standardize MAE by switching from custom identifiers into gDR-default
#'
#' @param mae a MultiAssayExperiment object with drug-response data generate by gDR pipeline
#' @param use_default boolean indicating whether or not to use default
#' identifiers for standardization
#' 
#' @return mae a MultiAssayExperiment with default gDR identifiers
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' S4Vectors::metadata(mae[[1]])$identifiers$drug <- "druug"
#' standardize_mae(mae)
#' 
#' @export
standardize_mae <- function(mae, use_default = TRUE) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  experiments <- names(mae)
  for (experiment in experiments) {
    mae[[experiment]] <- standardize_se(mae[[experiment]], use_default = use_default)
  }
  mae
}

#' Standardize SE by switching from custom identifiers into gDR-default
#'
#' @param se a SummarizedExperiment object with drug-response data generate by gDR pipeline
#' @param use_default boolean indicating whether or not to use default
#' identifiers for standardization
#'
#' @return se a SummarizedExperiment with default gDR identifiers
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' se <- mae[[1]]
#' S4Vectors::metadata(se)$identifiers$drug <- "druug"
#' standardize_se(se)
#' 
#' @export
standardize_se <- function(se, use_default = TRUE) {
  checkmate::assert_class(se, "SummarizedExperiment")
  
  reset_env_identifiers()
  idfs <- get_default_identifiers()
  idfs_se <- get_SE_identifiers(se)
  
  if (use_default) {
    from_idfs <- idfs_se
    to_idfs <- idfs
  } else {
    from_idfs <- idfs
    to_idfs <- idfs_se
  }
  
  # Extract matching names of identifiers
  matching_idfs <- .extract_matching_identifiers(to_idfs,
                                                 from_idfs)
  
  # Extract changed identifiers
  diff_identifiers <- .extract_diff_identifiers(matching_idfs$default,
                                                matching_idfs$se)
  
  diff_names <- unique(unlist(lapply(diff_identifiers, names)))
  
  if ("untreated_tag" %in% diff_names) {
    rowData(se) <- .replace_untreated_tag(rowData(se),
                                          to_idfs,
                                          from_idfs)
  }
  
  # Create mapping vector
  mapping_df <- do.call(rbind,
                        lapply(seq_along(diff_names),
                               function(x) {
                                 data.table::data.table(x = unlist(diff_identifiers$default[x]),
                                                        y = unlist(diff_identifiers$se[x]))
                               }
                        )
  )
  
  if (length(mapping_df)) {
    mapping_vector <- mapping_df$x
    names(mapping_vector) <- mapping_df$y
    
    # Replace rowData, colData and assays
    rowData(se) <- rename_DFrame(rowData(se), mapping_vector)
    colData(se) <- rename_DFrame(colData(se), mapping_vector)
    
    assayList <- lapply(assays(se), function(x) {
      rename_bumpy(x, mapping_vector)
    })
    assays(se) <- assayList
  }
  se
}

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
  untreated_tag <- data.table::data.table(x = default[["untreated_tag"]],
                                          y = se_identifiers[["untreated_tag"]])
  untreated_tag_mapping <- untreated_tag$x
  names(untreated_tag_mapping) <- untreated_tag$y
  S4Vectors::DataFrame(lapply(row_data, function(x) {
    if (is.character(x)) {
      stringr::str_replace_all(x, untreated_tag_mapping)
    } else {
      x
    }
  }), row.names = rownames(row_data), check.names = FALSE)
}

#' Rename DFrame
#'
#' @param df a DFrame object
#' @param mapping_vector a named vector for mapping old and new values.
#' The names of the character vector indicate the source names, and the
#' corresponding values the destination names.
#'
#' @return a renamed DFrame object
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' rename_DFrame(SummarizedExperiment::rowData(mae[[1]]), c("Gnumber" = "Gnumber1"))
#' 
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
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' se <- mae[[1]]
#' assay <- SummarizedExperiment::assay(se)
#' rename_bumpy(assay, c("Concentration" = "conc"))
#' 
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
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' refine_coldata(SummarizedExperiment::colData(mae[[1]]), mae[[1]])
#' 
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
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' refine_rowdata(SummarizedExperiment::colData(mae[[1]]), mae[[1]])
#' 
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
