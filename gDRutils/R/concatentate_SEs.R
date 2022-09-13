has_nested_field <- function(asy, nested_field) {
  df <- BumpyMatrix::unsplitAsDataFrame(asy)
  all(nested_field %in% colnames(df))
}

#' Demote a metadata field to be represented as a nested field of the \code{BumpyMatrix} assay.
#'
#' @param se \code{SummarizedExperiment} object.
#' @param fields Character vector of nested fields to demote as nested columns.
#'
#' @return A \code{SummarizedExperiment} object with new dimensions resulting from demoting given \code{fields}.
#'
#' @seealso promote_fields
#' @export
demote_fields <- function(se, fields) {
  rowmd <- colnames(rowData(se)) 
  colmd <- colnames(colData(se))

  if (any(nested_fields <- !fields %in% c(rowmd, colmd))) {
    stop(sprintf("fields '%s' are already demoted fields, perhapy you intended to call 'promote_fields'?",
      paste0(fields[nested_fields], collapse = ", ")))
  }
}

#' Promote a nested field to be represented as a metadata field of the \code{SummarizedExperiment}
#' as either the \code{rowData} or \code{colData}.
#'
#' @param se \code{SummarizedExperiment} object.
#' @param fields Character vector of nested fields to promote.
#' @param MARGIN Numeric of values \code{1} or \code{2} indicating whether to 
#' promote fields to rows or columns respectively.
#'
#' @return A \code{SummarizedExperiment} object with new dimensions resulting from promoting given \code{fields}.
#' @seealso demote_fields
#' @export
promote_fields <- function(se, fields, MARGIN = c(1, 2)) {
  MARGIN <- match.arg(MARGIN)

  rowmd <- colnames(rowData(se)) 
  colmd <- colnames(colData(se))

  if (any(rfields <- fields %in% rowmd) || any(cfields <- fields %in% colmd)) {
    stop(sprintf("fields '%s' are already promoted fields",
      paste0(fields[rfields || cfields], collapse = ", ")))
  }

  if (MARGIN == 1) {
    rowmd <- c(rowmd, fields)
  } else if (MARGIN == 2) {
    colmd <- c(colmd, fields)
  }

  final_assays <- list()
  asynames <- SummarizedExperiment::assayNames(se)
  for (asyname in asynames) {
    asy <- SummarizedExperiment::assay(se, asyname)
    if (!all(has_nested_field(asy, fields))) {
      warning(sprintf("dropping assay '%s' from the final SE", asyname))
    } else {
      df <- S4Vectors::DataFrame(as.data.frame(convert_se_assay_to_dt(se, asyname, include_metadata = TRUE)))
      final_assays[[asyname]] <- .transform_df_to_matrix(df = df,
        row_fields = rowmd,
        column_fields = colmd,
        nested_fields = setdiff(colnames(df), c(rowmd, colmd))) 
    }
  }

  if (length(final_SE) == 0L) {
    stop(sprintf("unable to find nested fields: '%s', perhaps you intended to call 'demote_fields'?",
      paste0(fields, collapse = ", "))) 
  }

  # Ensure order.
  final_assays <- lapply(final_assays, function(x) {x[rownames(rowData), colnames(colData)]})
  final_SE <- SummarizedExperiment::SummarizedExperiment(
    assays = final_assays,
    rowData = ,
    colData = ,
    metadata = S4Vectors::metadata(se)
  )
  final_SE
}

#' Transform a \code{DFrame} to a \code{BumpyMatrix}.
#'
#' This differs from the plain \code{BumpyMatrix::splitAsBumpyMatrix}
#' as it additionally identifies which rows belong in which element of the
#' abstracted matrix based on the arguments \code{row_fields}, \code{column_fields},
#' and \code{nested_fields}.
#'
#' @param df \code{DataFrame} object.
#' @param row_fields Character vector representing the columns that should become \code{rowData}.
#' @param column_fields Character vector representing the columns that should become \code{colData}.
#' @param nested_fields Character vector representing the columns that should become nested in the
#' final matrix object.
#'
#' @return A \code{BumpyMatrix} object.
#' @seealso promote_fields demote_fields 
#'
#' @noRd
#' @importFrom BumpyMatrix splitAsBumpyMatrix
#' @keywords internal
.transform_df_to_matrix <- function(df, row_fields, column_fields, nested_fields) {
  missing <- setdiff(colnames(df), c(row_fields, column_fields, nested_fields))
  if (length(missing) != 0L) {
    stop(sprintf("found columns in 'df' not specified as row, column, or nested fields",
      paste0(missing, collapse = ", ")))
  }

  df$row <- S4Vectors::selfmatch(df[, row_fields])
  df$column <- S4Vectors::selfmatch(df[, column_fields])
  BumpyMatrix::splitAsBumpyMatrix(df[, nested_fields], row = df$row, column = df$column)
}

#' Aggregate a \code{BumpyMatrix} assay by a given aggreation function.
#' 
#' Aggregation can only be performed on nested variables.
#' 
#' @param asy A \code{BumpyMatrix} object.
#' @param by Character vector of the nested fields to aggregate by.
#' @param FUN A function to use to aggregate the data.
#' 
#' @return A \code{BumpyMatrix} object aggregated by \code{FUN}.
#' @importFrom stats aggregate
#' @importFrom BumpyMatrix unsplitAsDataFrame splitAsBumpyMatrix
#' @export
aggregate_assay <- function(asy, by, FUN) {
  checkmate::assert_class(asy, "BumpyMatrix")
  checkmate::assert_class(FUN, "function")
  checkmate::assert_class(by, "character")

  row.field <- "row"
  column.field <- "column"

  df <- BumpyMatrix::unsplitAsDataFrame(asy, row.names = TRUE, column.names = TRUE, row.field = row.field, column.field = column.field)
  if (!all(present <- by %in% setdiff(colnames(df), c(row.field, column.field)))) {
    stop(sprintf("specified 'by' columns: '%s' are not present in 'asy'", paste0(by[!present], collapse = ", ")))
  }

  by <- c(row.field, column.field, by)
  agg <- stats::aggregate(x = data.frame(df[, setdiff(colnames(df), by)]), by = data.frame(df[, by]), FUN = FUN, na.rm = TRUE)
  mat <- BumpyMatrix::splitAsBumpyMatrix(agg[, setdiff(colnames(agg), c(row.field, column.field))], row = agg[[row.field]], column = agg[[column.field]])

  mat[rownames(asy), colnames(asy)] # Ensure order is the same.
}
