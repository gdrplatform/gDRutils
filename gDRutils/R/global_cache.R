#' @include identifiers_list.R

## This function maintains a cache of identifiers, headers, and their respective values at run time. 
## It allows users to set identifier values using the 'set_identifier' function.
make_global_cache <- function() {
  identifiers_list <- IDENTIFIERS_LIST
  headers_list <- list()

  valid_ids <- names(IDENTIFIERS_LIST) 

  .get_id <- function(k = NULL) {
    if (!is.null(k)) {
      checkmate::assert_string(k, null.ok = TRUE)
      checkmate::assert_choice(k, choices = valid_ids)

      if (! k %in% valid_ids) {
	stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	  k, paste0(valid_ids, collapse=", ")))
      }

      return(identifiers_list[[k]])
    } else {
      return(identifiers_list)
    }
  }

  .set_id <- function(k, v) {
    if (! k %in% valid_ids) {
      stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	k, paste0(valid_ids, collapse=", ")))
    }
    identifiers_list[[k]] <<- v
  }

  .get_header <- function(k = NULL) {
    ## The following .getHeadersList() call is inside the function .get_header
    ## to avoid cyclical dependencies and collation order problems.

    if (length(headers_list) == 0L) {
      assign(x = "headers_list", value = .getHeadersList(), pos = parent.frame())
    }

    if (!is.null(k)) {
      valid_headers <- names(HEADERS_LIST) 

      checkmate::assert_string(k, null.ok = TRUE)
      checkmate::assert_choice(k, choices = valid_headers)

      if (! k %in% valid_headers) {
	stop(sprintf("'%s' is not one of the valid headers: '%s'", 
	  k, paste0(valid_headers, collapse=", ")))
      }

      return(headers_list[[k]])
    } else {
      return(headers_list)
    }
  }

  .reset_ids <- function() {
    assign(x = identifiers_list, value = list(), pos = parent.frame())
  }

  .reset_headers <- function() {
    assign(x = headers_list, value = list(), pos = parent.frame())
  }

  list(get_id=.get_id, set_id=.set_id, get_head=.get_header, reset_ids=.reset_ids, reset_heads=.reset_headers)
}

global_cache <- make_global_cache()
