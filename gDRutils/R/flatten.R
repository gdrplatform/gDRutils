#' Flatten a table
#'
#' Flatten a stacked table into a wide format.
#'
#' @param tbl table to flatten.
#' @param groups character vector of column names representing uniquifying groups in expansion.
#' @param wide_cols character vector of column names to flatten.
#' @param sep string representing separator between \code{wide_cols} columns, used in column renaming.
#' Defaults to \code{"_"}.
#'
#' @return table of flattened data as defined by \code{wide_cols}.
#'
#' @details flattened columns will be named with original column names prefixed by \code{wide_cols} columns,
#' concatenated together and separated by \code{sep}.
#'
#' A common use case for this function is 
#' when a flattened version of the \code{"Metrics"} assay is desired.
#'
#' @examples
#' n <- 4
#' m <- 5
#' grid <- expand.grid(normalization_type = c("GR", "RV"),
#'   source = c("GDS", "GDR"))
#' repgrid <- do.call("rbind", rep(list(grid), m))
#' repgrid$wide <- seq(m * n)
#' repgrid$id <- rep(LETTERS[1:m], each = n)
#'
#' groups <- colnames(grid)
#' wide_cols <- c("wide")
#'
#' flatten(repgrid, groups = groups, wide_cols = wide_cols)
#'
#' @seealso convert_se_assay_to_dt
#' @export
#'
flatten <- function(tbl, groups, wide_cols, sep = "_") {

  checkmate::assert_character(groups)
  checkmate::assert_character(wide_cols)
  checkmate::assert_string(sep)
  checkmate::assert_true(
    is(tbl, "data.table") ||
      is(tbl, "data.frame") ||
        is(tbl, "DFrame")
  )

  if (!all(groups %in% colnames(tbl))) {
    stop(sprintf("missing expected uniquifying groups: '%s'",
      paste0(setdiff(groups, colnames(tbl)), collapse = ", ")))
  }	
  
  idx <- which(colnames(tbl) %in% groups)
  uniquifying <- subset(tbl, select = idx)
  uniquifying <- unique(uniquifying)

  out <- split(subset(tbl, select = -idx), subset(tbl, select = idx), sep = sep)
  missing <- setdiff(wide_cols, colnames(tbl))
  if (length(missing) != 0L) {
    warning(sprintf("missing listed wide_cols columns: '%s'", paste0(missing, collapse = ", ")))
  }

  rename <- colnames(out[[1]]) %in% wide_cols 
  for (grp in names(out)) {
    group <- out[[grp]]
    colnames(group)[rename] <- paste0(grp, sep, colnames(group)[rename])
    out[[grp]] <- group
  }

  ## Drop empty elements for successful merge.
  filtered <- out[lapply(out, nrow) > 0L]
  Reduce(function(x, y) merge(x, y, by = intersect(names(x), names(y))), filtered)
}
