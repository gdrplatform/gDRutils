#' df_to_bm_assay
#'
#' Convert data.table with dose-response data into a BumpyMatrix assay.
#'
#' The 'assay' is simply a BumpyMatrix object with rownames being the treatment ids,
#' colnames being the ids of the cell lines and values with
#' dose-response data for given cell lines under given conditions.
#'
#' @param data data.table with drug-response data
#' @param discard_keys a vector of keys that should be discarded
#' @param precomputed_metadata optional named list from a prior
#'   \code{\link{split_SE_components}} call to avoid recomputation
#' @keywords convert
#'
#' @examples
#' df_to_bm_assay(data.table::data.table(Gnumber = 2, clid = "A"))
#'
#' @return BumpyMatrix object
#'
#' @export
df_to_bm_assay <-
  function(data,
           discard_keys = NULL,
           precomputed_metadata = NULL) {
    stopifnot(any(inherits(data, "data.table"), checkmate::test_character(data)))
    checkmate::assert_character(discard_keys, null.ok = TRUE)
    checkmate::assert_list(precomputed_metadata, null.ok = TRUE)

    required <- c("condition_md", "treatment_md", "data_fields")
    if (!is.null(precomputed_metadata) &&
        checkmate::test_names(names(precomputed_metadata),
                              must.include = required)) {
      allMetadata <- precomputed_metadata
    } else {
      allMetadata <- split_SE_components(data, nested_keys = discard_keys)
    }

    seColData <- data.table::as.data.table(allMetadata$condition_md)
    seRowData <- data.table::as.data.table(allMetadata$treatment_md)

    # Create row_id and col_id directly in data.table for speed
    seColData[, col_id := .I]
    seRowData[, row_id := .I]

    cl_entries <- setdiff(colnames(seColData), c("col_id", "name_"))
    cond_entries <- setdiff(colnames(seRowData), c("row_id", "name_"))

    # Use data.table's cross-join for speed
    complete <- data.table::CJ(col_id = seColData$col_id, row_id = seRowData$row_id)
    complete[, factor_id := .I]

    # Use data.table merge for speed
    completeMerged <- merge(complete, seColData, by = "col_id", all.x = TRUE)
    completeMerged <- merge(completeMerged, seRowData, by = "row_id", all.x = TRUE)

    # Merge with original data using data.table merge
    data_assigned <- merge(data, completeMerged, by = c(cond_entries, cl_entries), all.x = TRUE)

    # Order by factor_id using data.table's setorder
    data.table::setorder(data_assigned, factor_id)

    bm <- BumpyMatrix::splitAsBumpyMatrix(data_assigned[, allMetadata$data_fields, with = FALSE],
                                          row = data_assigned$row_id,
                                          column = data_assigned$col_id)

    # Check if the order of rows/cols are correct
    stopifnot(!is.unsorted(as.numeric(rownames(bm))))
    stopifnot(!is.unsorted(as.numeric(colnames(bm))))

    bm
  }


#' Build a SummarizedExperiment from a long table
#'
#' Turns a long table, one row per measurement, into an SE with a single
#' BumpyMatrix assay. Meant for assays that do not come out of the gDR pipeline,
#' e.g. compounds x targets in a biochemical screen; register a profile with
#' \code{\link{register_fit_profile}} and pass the result to \code{\link{apply_fit}}.
#'
#' Several \code{row_cols} or \code{col_cols} are pasted with \code{"_"} into the
#' row or column names. Combinations absent from \code{dt} come out as empty cells.
#'
#' @param dt data.table with one row per measurement.
#' @param row_cols character vector of columns identifying an SE row.
#' @param col_cols character vector of columns identifying an SE column.
#' @param assay_name string; name of the assay.
#'
#' @return \code{SummarizedExperiment} with \code{row_cols} in \code{rowData},
#'   \code{col_cols} in \code{colData} and the remaining columns nested in
#'   assay \code{assay_name}.
#'
#' @examples
#' dt <- data.table::CJ(compound = c("cA", "cB"), target = "t1", Concentration = c(0.1, 1))
#' dt$Activity <- c(90, 40, 80, 20)
#' convert_dt_to_se(dt, "compound", "target", "Biochemical")
#'
#' @keywords convert
#' @export
convert_dt_to_se <- function(dt, row_cols, col_cols, assay_name) {
  checkmate::assert_data_table(dt, min.rows = 1L)
  checkmate::assert_character(row_cols, min.len = 1L)
  checkmate::assert_character(col_cols, min.len = 1L)
  checkmate::assert_disjunct(row_cols, col_cols)
  checkmate::assert_subset(c(row_cols, col_cols), names(dt))
  checkmate::assert_string(assay_name)

  data_cols <- setdiff(names(dt), c(row_cols, col_cols))
  if (length(data_cols) == 0L) {
    stop("dt has no columns besides row_cols and col_cols")
  }

  row_data <- unique(dt[, row_cols, with = FALSE])
  col_data <- unique(dt[, col_cols, with = FALSE])
  row_names <- .paste_id_cols(row_data)
  col_names <- .paste_id_cols(col_data)
  clashing <- c(row_names[duplicated(row_names)], col_names[duplicated(col_names)])
  if (length(clashing) > 0L) {
    stop(sprintf("identifier columns paste to the same name: %s", toString(unique(clashing))))
  }

  bm <- BumpyMatrix::splitAsBumpyMatrix(
    dt[, data_cols, with = FALSE],
    row = factor(.paste_id_cols(dt[, row_cols, with = FALSE]), levels = row_names),
    column = factor(.paste_id_cols(dt[, col_cols, with = FALSE]), levels = col_names)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(list(bm), assay_name),
    rowData = S4Vectors::DataFrame(row_data, row.names = row_names),
    colData = S4Vectors::DataFrame(col_data, row.names = col_names)
  )
}

#' @keywords internal
.paste_id_cols <- function(id_dt) {
  do.call(paste, c(as.list(id_dt), sep = "_"))
}
