has_nested_field <- function(asy, nested_field) {
  checkmate::assertClass(asy, "BumpyMatrix")
  checkmate::assertCharacter(nested_field)

  df <- BumpyMatrix::unsplitAsDataFrame(asy)
  all(nested_field %in% colnames(df))
}

.validate_final_assays <- function(final_assays) {
  checkmate::assert_list(final_assays, types = ("list"))
  if (length(final_assays) > 1L) {
    exp <- final_assays[[1]]
    for (i in 2:length(final_assays)) {
      curr_assay <- final_assays[[i]]
      stopifnot(rownames(exp$rowData) == rownames(curr_assay$rowData))
      stopifnot(rownames(exp$colData) == rownames(curr_assay$colData))
    }
  }
  invisible(NULL)
}

#' Demote a metadata field in the \code{rowData} or \code{colData} of a \code{SummarizedExperiment} object
#' to a nested field of a \code{BumpyMatrix} assay.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param fields Character vector of metadata fields to demote as nested columns.
#' @keywords concatenate_SEs
#'
#' @return A \code{SummarizedExperiment} object with new dimensions resulting from demoting given \code{fields}
#' to nested columns.
#'
#' @seealso promote_fields
#' @details Revert this operation using \code{promote_fields}.
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' se <- promote_fields(se, "ReadoutValue", 2)
#' demote_fields(se, "ReadoutValue")
#' 
#' @export
demote_fields <- function(se, fields) {
  checkmate::assertClass(se, "SummarizedExperiment")
  checkmate::assertCharacter(fields)

  rowmd <- colnames(rowData(se))
  colmd <- colnames(colData(se))

  if (any(are_nested_fields <- !fields %in% c(rowmd, colmd))) {
    stop(sprintf("field(s) '%s' are already demoted fields, perhapy you intended to call 'promote_fields'?",
                 paste0(fields[are_nested_fields], collapse = ", ")))
  }

  rowmd <- setdiff(rowmd, fields)
  colmd <- setdiff(colmd, fields)

  asynames <- SummarizedExperiment::assayNames(se)
  final_assays <- vector("list", length(asynames))
  names(final_assays) <- asynames
  for (asyname in asynames) {
    asy <- SummarizedExperiment::assay(se, asyname)
    df <- S4Vectors::DataFrame(convert_se_assay_to_dt(se, asyname, include_metadata = TRUE))
    df <- df[, setdiff(colnames(df), c("rId", "cId"))]
    final_assays[[asyname]] <- .transform_df_to_matrix(df = df,
                                                       row_fields = rowmd,
                                                       column_fields = colmd,
                                                       nested_fields = setdiff(colnames(df), c(rowmd, colmd)))
  }

  .validate_final_assays(final_assays)

  assay_list <- lapply(final_assays, function(x) {
    x$mat
  })

  SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    rowData = final_assays[[1]]$rowData,
    colData = final_assays[[1]]$colData,
    metadata = S4Vectors::metadata(se)
  )
}

#' Promote a nested field to be represented as a metadata field of the \code{SummarizedExperiment}
#' as either the \code{rowData} or \code{colData}.
#'
#' @param se \code{SummarizedExperiment} object.
#' @param fields Character vector of nested fields in a \code{BumpyMatrix} object to promote to
#' metadata fields of a \code{se}.
#' @param MARGIN Numeric of values \code{1} or \code{2} indicating whether to
#' promote fields to rows or columns respectively.
#' @keywords concatenate_SEs
#'
#' @return A \code{SummarizedExperiment} object with new dimensions resulting from promoting given \code{fields}.
#' @details Revert this operation using \code{demote_fields}.
#' @seealso demote_fields
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' se <- promote_fields(se, "ReadoutValue", 2)
#' 
#' @export
promote_fields <- function(se, fields, MARGIN = c(1, 2)) {
  checkmate::assertClass(se, "SummarizedExperiment")
  checkmate::assertCharacter(fields)
  if (length(MARGIN) != 1L || !MARGIN %in% c(1, 2)) {
    stop("invalid 'MARGIN' argument, must be either '1' or '2'")
  }
  rowmd <- colnames(rowData(se)) 
  colmd <- colnames(colData(se))
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
  } else {
    stop("invalid input for argument 'MARGIN'")
  }
  final_assays <- list()
  asynames <- SummarizedExperiment::assayNames(se)
  for (asyname in asynames) {
    asy <- SummarizedExperiment::assay(se, asyname)
    if (!all(has_nested_field(asy, fields))) {
      warning(sprintf("dropping assay '%s' from the final SE", asyname))
    } else {
      df <- S4Vectors::DataFrame(convert_se_assay_to_dt(se, asyname, include_metadata = TRUE))
      df <- df[, setdiff(colnames(df), c("rId", "cId"))]
      final_assays[[asyname]] <- .transform_df_to_matrix(df = df,
                                                         row_fields = rowmd,
                                                         column_fields = colmd,
                                                         nested_fields = setdiff(colnames(df), c(rowmd, colmd)))
    }
  }
  if (length(final_assays) == 0L) {
    stop(sprintf("unable to find nested fields: '%s' in any assays, perhaps you intended to call 'demote_fields'?",
                 paste0(fields, collapse = ", ")))
  }
  .validate_final_assays(final_assays)
  assay_list <- lapply(final_assays, function(x) {
    x$mat
  })
  SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    rowData = final_assays[[1]]$rowData,
    colData = final_assays[[1]]$colData,
    metadata = S4Vectors::metadata(se)
  )
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
#' @return A named list of \code{mat}, \code{rowData}, and \code{colData}
#' where \code{mat} is a \code{BumpyMatrix} object, ordered by the rownames of
#' \code{rowData} and \code{colData}.
#'
#' @seealso promote_fields demote_fields
#'
#' @noRd
#' @keywords internal
.transform_df_to_matrix <- function(df, row_fields, column_fields, nested_fields) {
  checkmate::assertClass(df, "DFrame")
  checkmate::assertCharacter(row_fields)
  checkmate::assertCharacter(column_fields)
  checkmate::assertCharacter(nested_fields)

  missing <- setdiff(colnames(df), c(row_fields, column_fields, nested_fields))
  if (length(missing) != 0L) {
    stop(sprintf("found columns in 'df' not specified as row, column, or nested fields",
                 paste0(missing, collapse = ", ")))
  }
  if (any(is_duplicated <- duplicated(c(row_fields, column_fields, nested_fields)))) {
    stop(sprintf("fields: '%s' are duplicated across arguments 'row_fields', 'column_fields', 'nested_fields'",
                 c(row_fields, column_fields, nested_fields)[is_duplicated]))
  }


  df$row <- S4Vectors::selfmatch(df[, row_fields])
  rowData <- unique(df[, c("row", row_fields)])
  rownames(rowData) <- rowData$row # Set rownames to maintain association btw rowData and mat.

  df$column <- S4Vectors::selfmatch(df[, column_fields])
  colData <- unique(df[, c("column", column_fields)])
  rownames(colData) <- colData$column # Set rownames to maintain association btw colData and mat.

  mat <- BumpyMatrix::splitAsBumpyMatrix(df[, nested_fields], row = df$row, column = df$column)
  list(mat = mat[rownames(rowData), rownames(colData)],
       rowData = rowData[, setdiff(colnames(rowData), "row")],
       colData = colData[, setdiff(colnames(colData), "column")])
}

#' Aggregate a \code{BumpyMatrix} assay by a given aggreation function.
#'
#' Aggregation can only be performed on nested variables.
#'
#' @param asy A \code{BumpyMatrix} object.
#' @param by Character vector of the nested fields to aggregate by.
#' @param FUN A function to use to aggregate the data.
#' @keywords concatenate_SEs
#'
#' @return A \code{BumpyMatrix} object aggregated by \code{FUN}.
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small") 
#' se <- mae[[1]]
#' assay <- SummarizedExperiment::assay(se)
#' aggregate_assay(assay, FUN = mean, by = c("Barcode"))
#' 
#' @export
aggregate_assay <- function(asy, by, FUN) {
  checkmate::assert_class(asy, "BumpyMatrix")
  checkmate::assert_class(FUN, "function")
  checkmate::assertCharacter(by)

  row.field <- "row"
  column.field <- "column"

  df <- BumpyMatrix::unsplitAsDataFrame(asy,
                                        row.names = TRUE,
                                        column.names = TRUE,
                                        row.field = row.field,
                                        column.field = column.field)
  if (!all(present <- by %in% setdiff(colnames(df), c(row.field, column.field)))) {
    stop(sprintf("specified 'by' columns: '%s' are not present in 'asy'", paste0(by[!present], collapse = ", ")))
  }

  by <- c(row.field, column.field, by)
  dt <- data.table::as.data.table(df)
  agg <- dt[, lapply(.SD, FUN), by = by]
  mat <- BumpyMatrix::splitAsBumpyMatrix(agg[, setdiff(colnames(agg),
                                                       c(row.field, column.field)), with = FALSE],
                                         row = agg[[row.field]],
                                         column = agg[[column.field]])

  mat[rownames(asy), colnames(asy)] # Ensure order is the same.
}
