#' @include identifiers_list.R

## This function maintains a cache of identifiers, headers, and their respective values at run time. 
## It allows users to set identifier values using the 'set_identifier' function.
make_global_cache <- function() {
  identifiers_list <- list() 
  headers_list <- list()

  valid_ids <- names(IDENTIFIERS_LIST)

  .get_id <- function(k = NULL) {
    if (length(identifiers_list) == 0L) {
      assign(x = "identifiers_list", value = IDENTIFIERS_LIST, pos = parent.env(environment()))
    }

    if (!is.null(k)) {
      checkmate::assert_string(k, null.ok = TRUE)
      checkmate::assert_choice(k, choices = valid_ids)

      if (! k %in% valid_ids) {
	stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	  k, paste0(valid_ids, collapse=", ")))
      }

      return(identifiers_list[[k]])
    } else {
      missing_ids <- setdiff(names(IDENTIFIERS_LIST), names(identifiers_list))	
      if (length(missing_ids) != 0L) {	
        sapply(missing_ids, get_identifier) 	
      }

      return(identifiers_list)
    }
  }

  .set_id <- function(k, v) {
    checkmate::assert_string(k, null.ok = TRUE)
    checkmate::assert_choice(k, choices = valid_ids)

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
      assign(x = "headers_list", value = .getHeadersList(), pos = parent.env(environment()))
    }

    if (!is.null(k)) {
      valid_headers <- names(headers_list) 

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
    assign(x = "identifiers_list", value = list(), pos = parent.env(environment()))
  }

  .reset_headers <- function() {
    assign(x = "headers_list", value = list(), pos = parent.env(environment()))
  }

  list(get_id=.get_id, set_id=.set_id, get_head=.get_header, reset_ids=.reset_ids, reset_heads=.reset_headers)
}

global_cache <- make_global_cache()
