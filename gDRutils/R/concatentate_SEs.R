#' @param Named list
#' @return A \code{SummarizedExperiment} object.
#' @details Any metadata columns that exist in one experiment and not the other will be
#' appended in a 'sparse' manner, and \code{NA}s filled where the values are missing.
#' @export
concatenate_SEs <- function(ses) {
  checkmate::is_list(ses)
  lapply(ses, function(x) is(x, "SummarizedExperiment"))
  concatenate_assays
}

#' @param asy1
#' @param asy2
#' @export
concatenate_assays <- function(asy1, asy2)  {
  # TODO: Is there a chance we lose some values if they are all NA?
}

aggregate_SEs <- function(ses) {

}

#' Aggregate a \code{BumpyMatrix} assay by a given aggreation function.
#' 
#' Aggregation can only be performed on nested variables.
#' 
#' @param asy A \code{BumpyMatrix} object.
#' 
#' @return A \code{BumpyMatrix} object aggregated by \code{aggregation_fxn}.
#' @export
aggregate_assay <- function(asy, aggregation_fxn, by) {
  checkmate::assert_true(is(asy, "BumpyMatrix"))
  checkmate::assert_true(is(aggregation_fxn, "function"))
  checkmate::assert_true(is(by, "character"))

  row.field <- "row"
  column.field <- "column"

  df <- BumpyMatrix::unsplitAsDataFrame(asy, row.names = TRUE, column.names = TRUE, row.field = row.field, column.field = column.field)
  if (!all(present <- by %in% setdiff(colnames(df), c(row.field, column.field)))) {
    stop(sprintf("specified 'by' columns: '%s' are not present in 'asy'", paste0(by[!present], collapse = ", ")))
  }

  by <- c(row.field, column.field, by)
  agg <- stats::aggregate(x = data.frame(df[, setdiff(colnames(df), by)]), by = data.frame(df[, by]), FUN = aggregation_fxn, na.rm = TRUE)
  mat <- BumpyMatrix::splitAsBumpyMatrix(agg[, setdiff(colnames(agg), c(row.field, column.field))], row = agg[[row.field]], column = agg[[column.field]])

  mat[rownames(asy), colnames(asy)] # Ensure order is the same.
}
