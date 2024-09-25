#' Standardize MAE by switching from custom identifiers into gDR-default
#'
#' @param mae a MultiAssayExperiment object with drug-response data generate by gDR pipeline
#' @param use_default boolean indicating whether or not to use default
#' identifiers for standardization
#' @keywords standardize_MAE
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
#' @keywords standardize_MAE
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
#' @keywords standardize_MAE
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
#' @keywords standardize_MAE
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
#' @keywords standardize_MAE
#'
#' @return a charvec containing the names of the optional identifiers in the SE colData
#'
get_optional_coldata_fields <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  idfs <- get_SE_identifiers(se)
  
  as.character(idfs["cellline_tissue"])
}

#' get optional rowData fields
#'
#' @param se a SummarizedExperiment object with drug-response data generate by gDR pipeline
#' @keywords standardize_MAE
#'
#' @return a charvec containing the names of the optional identifiers in the SE rowData
#'
get_optional_rowdata_fields <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  idfs <- get_SE_identifiers(se)
  rowdata <- SummarizedExperiment::rowData(se)
  
  out <- c(idfs["drug_moa"])
  
  if (!is.null(rowdata[[idfs[["drug2"]]]])) {
    out <- c(out, idfs["drug_moa2"])
  }
  if (!is.null(rowdata[[idfs[["drug3"]]]])) {
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
#' @keywords standardize_MAE
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
#' @keywords standardize_MAE
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


#' Set Unique Parental Identifiers
#'
#' This function sets the `CellLineName` field in 
#' `colData` to be unique by appending the `clid` in parentheses for duplicates.
#'
#' @param se A SummarizedExperiment object.
#' @return A SummarizedExperiment object with unique `CellLineName` in `colData`.
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = matrix(1:4, ncol = 2)),
#'   colData = S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2"))
#' )
#' se <- set_unique_cl_names(se)
#' @export
#' @keywords standardize_MAE
set_unique_cl_names <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  
  col_data <- SummarizedExperiment::colData(se)
  col_data_new <- set_unique_cl_names_dt(col_data)
  SummarizedExperiment::colData(se) <- col_data_new
  
  se
}

#' Set Unique Parental Identifiers in table
#'
#' This function sets the `CellLineName` field in 
#' `colData` to be unique by appending the `clid` in parentheses for duplicates.
#'
#' @param col_data data.table or DFrame with row data
#' @param sep string with separator added before suffix
#' @return fixed input table with unique `CellLineName` in `colData`.
#' @examples
#' col_data <- S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2"))
#' col_data <- set_unique_cl_names_dt(col_data)
#' @export
#' @keywords standardize_MAE
#' 
set_unique_cl_names_dt <- function(col_data, sep = " ") {
  stopifnot(any(inherits(col_data, "data.table") || inherits(col_data, "DFrame")))
  
  cellline_name <- get_env_identifiers("cellline_name")
  clid <- get_env_identifiers("cellline")
  
  if (!is.null(col_data[[cellline_name]])) {
    if (data.table::is.data.table(col_data)) {
      # for data.table check multiple columns; and update cl names only if adding clid suffix can make them different
      unique_col_names <- c(unlist(get_default_identifiers()[
        c("cellline_name", "drug_name", "drug_name2",
          "concentration2", "duration", "data_source")
      ]), "normalization_type")
      unique_col_names <- intersect(unique_col_names, names(col_data))
      duplicated_ids <- col_data[[cellline_name]][duplicated(col_data, by = unique_col_names)]
      unique_col_names_clid <- c(unique_col_names, get_default_identifiers()$cellline)
      duplicated_ids_with_clid <- col_data[[cellline_name]][duplicated(col_data, by = unique_col_names_clid)]
      duplicated_ids <- setdiff(duplicated_ids, duplicated_ids_with_clid)
    } else {
      duplicated_ids <- col_data[[cellline_name]][duplicated(col_data[[cellline_name]])]
    }
    
    if (length(duplicated_ids) > 0) {
     for (dup_id in unique(duplicated_ids)) {
      dup_indices <- which(col_data[[cellline_name]] == dup_id)
      col_data[[cellline_name]][dup_indices] <- 
        paste0(
          col_data[[cellline_name]][dup_indices],
          sep,
          "(", 
          col_data[[clid]][dup_indices], 
          ")"
        )
      } 
    }
  }
  
  col_data
}

#' Set Unique Drug Names
#'
#' This function sets the `DrugName`, `DrugName_2`, and `DrugName_3` fields in `rowData`
#' to be unique by appending the corresponding `Gnumber`, `Gnumber_2`, and `Gnumber_3` in parentheses for duplicates.
#'
#' @param se A SummarizedExperiment object.
#' @return A SummarizedExperiment object with unique `DrugName` fields in `rowData`.
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = matrix(1:9, ncol = 3)),
#'   rowData = S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA", "DrugB"),
#'   Gnumber = c("G1", "G2", "G5"),
#'   DrugName_2 = c("DrugC", "DrugC", "DrugD"),
#'   Gnumber_2 = c("G3", "G4", "G5")
#' ))
#' se <- set_unique_drug_names(se)
#' @export
#' @keywords standardize_MAE
set_unique_drug_names <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  
  row_data <- SummarizedExperiment::rowData(se)
  row_data_new <- set_unique_drug_names_dt(row_data)
  
  SummarizedExperiment::rowData(se) <- row_data_new
  se
}


#' Set Unique Drug Names in table
#'
#' This function sets the `DrugName`, `DrugName_2`, and `DrugName_3` fields in 
#' `rowData` to be unique by appending the corresponding `Gnumber`, `Gnumber_2`, 
#' and `Gnumber_3` in parentheses for duplicates.
#'
#' @param row_data data.table or DFrame with row data
#' @param sep string with separator added before suffix
#' @return fixed input table with unique `DrugName` fields in `rowData`.
#' @examples
#' row_data <- S4Vectors::DataFrame(
#'   DrugName = c("DrugA", "DrugA", "DrugB"),
#'   Gnumber = c("G1", "G2", "G5"),
#'   DrugName_2 = c("DrugC", "DrugC", "DrugD"),
#'   Gnumber_2 = c("G3", "G4", "G5")
#' )
#' row_data <- set_unique_drug_names_dt(row_data)
#' @export
#' @keywords standardize_MAE
#' 
set_unique_drug_names_dt <- function(row_data, sep = " ") {
  stopifnot(any(inherits(row_data, "data.table") || inherits(row_data, "DFrame")))
  
  drug_columns <- intersect(unlist(get_env_identifiers(c("drug_name", "drug_name2", "drug_name3"), simplify = FALSE)),
                            names(row_data))
  gnumber_columns <- intersect(unlist(get_env_identifiers(c("drug", "drug2", "drug3"), simplify = FALSE)),
                               names(row_data))
  
  for (i in seq_along(drug_columns)) {
    drug_col <- drug_columns[i]
    gnumber_col <- gnumber_columns[i]
    
    if (!is.null(row_data[[drug_col]])) {
      # Find duplicated drug names
      duplicated_drugs <- row_data[[drug_col]][duplicated(row_data[[drug_col]])]
      unique_drugs <- unique(duplicated_drugs)
      
      for (dup_drug in unique_drugs) {
        dup_indices <- which(row_data[[drug_col]] == dup_drug)
        if (length(unique(row_data[[gnumber_col]][dup_indices])) > 1) {
          row_data[[drug_col]][dup_indices] <- 
            paste0(
              row_data[[drug_col]][dup_indices],
              sep,
              "(", 
              row_data[[gnumber_col]][dup_indices], 
              ")"
            )
        }
      }
    }
  }
  
  row_data
}

#' Set Unique Identifiers in MultiAssayExperiment
#'
#' This function sets the `CellLineName` in `colData` and `DrugName` fields in `rowData`
#' to be unique for each `SummarizedExperiment` in a `MultiAssayExperiment`.
#'
#' @param mae A MultiAssayExperiment object.
#' @return A MultiAssayExperiment object with unique identifiers.
#' @examples
#' se1 <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = matrix(1:4, ncol = 2)),
#'   colData = S4Vectors::DataFrame(parental_identifier = c("ID1", "ID1"), clid = c("C1", "C2")),
#'   rowData = S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA"), Gnumber = c("G1", "G2"))
#' )
#' rownames(SummarizedExperiment::colData(se1)) <- c("cl1", "cl2")
#' rownames(SummarizedExperiment::rowData(se1)) <- c("g1", "g")
#' se2 <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = matrix(5:8, ncol = 2)),
#'   colData = S4Vectors::DataFrame(parental_identifier = c("ID2", "ID2"), clid = c("C3", "C4")),
#'   rowData = S4Vectors::DataFrame(DrugName = c("DrugB", "DrugB"), Gnumber = c("G3", "G4"))
#' )
#' rownames(SummarizedExperiment::colData(se2)) <- c("cl3", "cl4")
#' rownames(SummarizedExperiment::rowData(se2)) <- c("g3", "g4")
#' mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(se1 = se1, se2 = se2))
#' mae <- set_unique_identifiers(mae)
#' @export
#' @keywords standardize_MAE
set_unique_identifiers <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  
  for (name in names(MultiAssayExperiment::experiments(mae))) {
    se <- MultiAssayExperiment::experiments(mae)[[name]]
    se <- set_unique_cl_names(se)
    se <- set_unique_drug_names(se)
    MultiAssayExperiment::experiments(mae)[[name]] <- se
  }
  
  return(mae)
}
