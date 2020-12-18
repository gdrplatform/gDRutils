#' @include identifier-constants.R

## This function maintains a cache of identifiers, headers, and their respective values at run time. 
## It allows users to set identifier values using the 'set_identifier' function.
make_global_cache <- function() {
  identifiers_list <- list()
  headers_list <- list()

  valid_ids <- names(IDENTIFIERS_LIST) 

  .get_id <- function(x = NULL) {
    if (!is.null(x)) {
      if (! x %in% valid_ids) {
	stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	  x, paste0(valid_ids, collapse=", ")))
      }

      v <- identifiers_list[[x]]
      if (is.null(v)) {
	identifiers_list[[x]] <<- IDENTIFIERS_LIST[[x]]
        v <- identifiers_list[[x]]
      }
      return(v)
    } else {
      missing_ids <- setdiff(names(IDENTIFIERS_LIST), names(identifiers_list))
      if (length(missing_ids) != 0L) {
        sapply(missing_ids, get_identifier) 
      }
      return(identifiers_list)
    }
  }

  .set_id <- function(x, v) {
    if (is.null(identifiers_list[[x]])) {
      if (! x %in% valid_ids) {
	stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	  x, paste0(valid_ids, collapse=", ")))
      }
      identifiers_list[[x]] <<- v
    } else {
      stop(sprintf("cannot set '%s' to '%s', already set to '%s'", x, v, identifiers_list[[x]]))
    }
  }

  .get_header <- function(x = NULL) {
    ## The following HEADERS_LIST call is needed inside the function
    ## because otherwise the collation order has cyclical dependencies.

    if (length(headers_list) == 0L) {
      HEADERS_LIST <<- .getHeadersList()
    }

    if (!is.null(x)) {
      valid_headers <- names(HEADERS_LIST) 
      if (! x %in% valid_headers) {
	stop(sprintf("'%s' is not one of the valid headers: '%s'", 
	  x, paste0(valid_headers, collapse=", ")))
      }

      v <- headers_list[[x]]
      if (is.null(v)) {
	headers_list[[x]] <<- HEADERS_LIST[[x]]
        v <- headers_list[[x]]
      }
      return(v)
    } else {
      missing_headers <- setdiff(names(HEADERS_LIST), names(headers_list))
      if (length(missing_headers) != 0L) {
        sapply(missing_headers, get_header) 
      }
      return(headers_list)
    }
  }

  list(get_id=.get_id, set_id=.set_id, get_head=.get_header)
}

global_cache <- make_global_cache()
