#' split_SE_components
#'
#' Divide the columns of an input data.table into treatment metadata, condition metadata,
#' experiment metadata, and assay data for further analysis. This will most commonly be used
#' to identify the different components of a \linkS4class{SummarizedExperiment} object.
#'
#' @param df_ data.table with drug-response data
#' @param nested_keys character vector of keys to exclude from the row or column metadata,
#' and to instead nest within an element of the matrix. See details.
#' @param combine_on integer value of \code{1} or \code{2}, indicating whether unrecognized columns
#' should be combined on row or column respectively.
#' Defaults to \code{1}.
#' @keywords SE_operators
#'
#' @examples
#' split_SE_components(data.table::data.table(clid = "CL1", Gnumber = "DrugA"))
#'
#' @return named list containing different elements of a \linkS4class{SummarizedExperiment};
#' see details.
#'
#' @details
#' Named list containing the following elements:
#' \itemize{
#'  \item{"treatment_md": }{treatment metadata}
#'  \item{"condition_md": }{condition metadata}
#'  \item{"data_fields": }{all data.table column names corresponding to fields nested within a BumpyMatrix cell}
#'  \item{"experiment_md": }{metadata that is constant for all entries of the data.table}
#'  \item{"identifiers_md": }{key identifier mappings}
#' }
#'
#' The \code{nested_keys} provides the user the opportunity to specify that they would not
#' like to use that metadata field as a differentiator of the treatments, and instead, incorporate it
#' into a nested \code{DataFrame} in the BumpyMatrix matrix object.
#'
#' In the event that if any of the \code{nested_keys} are constant throughout the whole data.table,
#' they will still be included in the DataFrame of the BumpyMatrix and not in the experiment_metadata.
#'
#' Columns within the \code{df_} will be identified through the following logic:
#' First, the known data fields and any specified \code{nested_keys} are extracted.
#' Following that, known cell and drug metadata fields are detected,
#' and any remaining columns that represent constant metadata fields across all rows are extracted.
#' Next, any cell line metadata will be heuristically extracted.
#' Finally, all remaining columns will be combined on either the rows or columns as specified by
#' \code{combine_on}.
#' 
#' @export
#'
split_SE_components <- function(df_, nested_keys = NULL, combine_on = 1L) {
  stopifnot(any(inherits(df_, "data.table"), inherits(df_, "DataFrame")))
  checkmate::assert_character(nested_keys, null.ok = TRUE)
  checkmate::assert_choice(combine_on, c(1, 2))
  nested_keys <- .clean_key_inputs(nested_keys, colnames(df_))
  identifiers_md <- get_env_identifiers(simplify = TRUE)
  identifiers_md$nested_keys <- nested_keys

  df_ <- S4Vectors::DataFrame(df_, check.names = FALSE)
  all_cols <- colnames(df_)
  # Identify known data fields.
  data_fields <- c(get_header("raw_data"), get_header("normalized_results"),
                   get_header("averaged_results"),
    get_header("metrics_results"), get_env_identifiers("concentration", simplify = TRUE),
    identifiers_md$well_position, identifiers_md$template, nested_keys,
    get_header("scores"), get_header("excess"), get_header("isobolograms"))
  data_fields <- unique(data_fields)
  data_cols <- data_fields[data_fields %in% all_cols]
  md_cols <- setdiff(all_cols, data_cols) 
  md <- unique(df_[, md_cols]) 
  colnames_list <- .extract_colnames(identifiers_md, md_cols)
  remaining_cols <- colnames_list$remaining_cols
  cell_cols <- colnames_list$cell_cols
  drug_cols <- colnames_list$drug_cols

  singletons <- vapply(remaining_cols,
    function(x) {
      nrow(unique(md[, x, drop = FALSE])) == 1L
      },
    logical(1))
  # Get experiment columns.
  constant_cols <- remaining_cols[singletons]
  exp_md <- unique(md[, constant_cols, drop = FALSE])
  remaining_cols <- remaining_cols[!singletons]
  # Identify cellline properties by checking what columns have only a 1:1 mapping for each cell line.
  cl_entries <- identify_linear_dependence(md[c(unname(unlist(cell_cols)), remaining_cols)], identifiers_md$cellline)
  remaining_cols <- setdiff(remaining_cols, cl_entries)
  md_list <- .combine_drug_and_trt_cols(md, drug_cols, cell_cols, combine_on, cl_entries, remaining_cols)
  out <- list(
    condition_md = md_list$condition_md,
    treatment_md = md_list$treatment_md,
    data_fields = data_cols,
    experiment_md = exp_md,
    identifiers_md = identifiers_md
  )
  out
}


#' @keywords internal
.extract_colnames <- function(identifiers_md, md_cols) {
  req_identifier_names <- c("drug", "drug_name", "drug_moa", "duration")
  req_identifiers_idx <- unique(unlist(lapply(req_identifier_names, function(x) grep(x, names(identifiers_md)))))
  drug_md <- identifiers_md[req_identifiers_idx]
  drug_cols <- intersect(drug_md, md_cols)
  cell_id <- identifiers_md$cellline
  cell_fields <- c(cell_id, get_header("add_clid"), identifiers_md["replicate"])
  cell_cols <- cell_fields[cell_fields %in% md_cols]
  
  remaining_cols <- setdiff(md_cols, c(drug_cols, cell_cols))
  list(
    remaining_cols = remaining_cols,
    cell_cols = cell_cols,
    drug_cols = drug_cols
  )
}

#' @keywords internal
.combine_drug_and_trt_cols <- function(md, drug_cols, cell_cols, combine_on, cl_entries, remaining_cols) {
  trt_cols <- drug_cols
  cell_cols <- unique(c(cell_cols, cl_entries))
  if (combine_on == 1L) {
    trt_cols <- c(trt_cols, remaining_cols)
  } else if (combine_on == 2L) {
    cell_cols <- c(cell_cols, remaining_cols)
  } else {
    stop(sprintf("combine_on input: '%s' of class: '%s' is not supported",
                 combine_on, class(combine_on)))
  }
  
  list(
    condition_md = add_rownames_to_metadata(md, cell_cols),
    treatment_md = add_rownames_to_metadata(md, trt_cols)
  )
}


#' @keywords internal
identify_linear_dependence <- function(df, identifier) {
  if (!identifier %in% names(df)) {
    stop(sprintf("identifier: '%s' not found, but is required for calculating linear dependence", identifier))
  }
  entries <- identifier
  for (j in setdiff(colnames(df), identifier)) {
    prop <- split(df[[j]], as.factor(df[, identifier]))
    prop <- lapply(prop, function(grp) {
      length(unique(grp)) == 1L
      })
    if (all(unlist(prop))) {
      entries <- c(entries, j)
    }
  }
  entries
}


#' @keywords internal
add_rownames_to_metadata <- function(md, cols) {
  md <- unique(md[, unname(unlist(cols)), drop = FALSE])
  rownames(md) <- apply(md, 1, function(x) {
    paste(x, collapse = "_")
    })
  md <- md[! names(md) %in% c("unique_id")]
  md
}
