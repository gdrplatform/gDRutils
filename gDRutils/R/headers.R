#' @include global_cache.R

#' @title Get or reset headers for one or all header field(s) respectively
#'
#' @description Get the expected header(s) for one field or reset all header fields
#'
#' @return 
#' For \code{reset_headers} a \code{NULL} invisibly.
#' For \code{get_header} a character vector of headers for field \code{k}.
#' @examples
#' \dontrun{
#' get_header(k = NULL)
#' reset_headers()
#' set_identifier("duration", "X_TIME")
#' get_header(k = NULL)
#' }
#' @name headers
NULL

#' @param k string of field (data type) to return headers for
#' 
#' @details
#' If \code{get_header} is called with no values, the entire available header list is returned.
#' @examples
#' \dontrun{
#' get_header("manifest")
#' }
#' @rdname headers
#' @export
#' 
get_header <- function(k = NULL) {
  checkmate::assert_string(k, null.ok = TRUE)
  global_cache$get_head(k)
}


#' @details
#' \code{reset_headers} should be used if processing multiple files within the same R session
#' that use different file headers.  
#' @rdname headers
#' @export
#' 
reset_headers <- function() {
  global_cache$reset_heads()
  invisible(NULL)
}
