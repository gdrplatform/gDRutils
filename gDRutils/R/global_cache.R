#' @include identifiers_list.R

## This function maintains a cache of identifiers, headers, and their respective values at run time. 
## It allows users to set identifier values using the 'set_identifier' function.
make_global_cache <- function() {
  identifiers_list <- IDENTIFIERS_LIST
  headers_list <- list()

  valid_ids <- names(identifiers_list)

  .get_id <- function(x = NULL) {
    if (!is.null(x)) {
      if (! x %in% valid_ids) {
	stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	  x, paste0(valid_ids, collapse=", ")))
      }

      return(identifiers_list[[x]])
    } else {
      return(identifiers_list)
    }
  }

  .set_id <- function(x, v) {
    if (! x %in% valid_ids) {
      stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	x, paste0(valid_ids, collapse=", ")))
    }
    identifiers_list[[x]] <<- v
  }

  .get_header <- function(x = NULL) {
    if (length(headers_list) == 0L) {
      headers_list <<- .getHeadersList()
    }

    if (!is.null(x)) {
      if (! x %in% names(headers_list)) {
	stop(sprintf("'%s' is not one of the valid headers: '%s'", 
	  x, paste0(names(headers_list), collapse=", ")))
      }

      return(headers_list[[x]])
    } else {
      return(headers_list)
    }
  }

  list(get_id=.get_id, set_id=.set_id, get_head=.get_header)
}

global_cache <- make_global_cache()
