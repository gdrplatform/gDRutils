#' @export
.clean_key_inputs <- function(keys, cols) {
  dropped <- setdiff(keys, cols)
  if (length(dropped) != 0L) {
    warning(sprintf("ignoring input keys: '%s' which are not present in data.frame",
      paste0(dropped, collapse = ", ")))
  }
  intersect(keys, cols)
}
