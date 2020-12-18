#' @include global_cache.R

#' Get the headers for a header field
#'
#' Get the expected headers for a specific data type
#' 
#' @param x string of field (data type) to return headers for
#' 
#' @return character vector of headers for field \code{x}
#' 
#' @details
#' If \code{get_header} is called with no values, the entire available header list is returned.
#' @examples
#' get_header("manifest")
#' @export
get_header <- function(x = NULL) {
  global_cache$get_head(x)
}
