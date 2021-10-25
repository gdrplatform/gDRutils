#' @export
.clean_key_inputs <- function(keys, cols) {
  dropped <- setdiff(keys, cols)
  if (length(dropped) != 0L) {
    warning(sprintf("ignoring input keys: '%s' which are not present in data.frame",
      paste0(dropped, collapse = ", ")))
  }
  intersect(keys, cols)
}


#' @noRd
#' @keywords internal
assert_equal_input_len <- function(outlier, ...) {
  first <- list(...)[[1]]
  h <- all(sapply(list(...), length) == length(first))
  if (!h) {
    stop("unequal length objects provided as input")
  }

  contains_length_one <- length(first) == 1L || length(outlier) == 1L
  if (length(first) != length(outlier) && !contains_length_one) {
    stop("unequal lengths detected, either the fit parameters must be length one, or the tested value")
  }

  invisible(NULL)
}
