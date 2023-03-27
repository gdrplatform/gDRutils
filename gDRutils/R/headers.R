#' @title Get or reset headers for one or all header field(s) respectively
#'
#' @description Get the expected header(s) for one field or reset all header fields
#'
#' @return 
#' For \code{get_header} a character vector of headers for field \code{k}.
#'
#' @examples
#' get_header(k = NULL)
#' @name headers
NULL

#' @param k string of field (data type) to return headers for
#' 
#' @details
#' If \code{get_header} is called with no values, the entire available header list is returned.
#' @examples
#' get_header("manifest")
#' @rdname headers
#' @export
#' 
get_header <- function(k = NULL) {
  checkmate::assert_string(k, null.ok = TRUE)
  .get_header(k)
}
