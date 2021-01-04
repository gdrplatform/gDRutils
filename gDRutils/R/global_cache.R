#' @include identifiers_list.R

## This function maintains a cache of identifiers, headers, and their respective values at run time. 
## It allows users to set identifier values using the 'set_identifier' function.
make_global_cache <- function() {
  cache <- new.env()
  cache$identifiers_list <- list() 
  cache$headers_list <- list()

  valid_ids <- names(IDENTIFIERS_LIST)

  .get_id <- function(k = NULL) {
    if (length(cache$identifiers_list) == 0L) {
      cache$identifiers_list <- IDENTIFIERS_LIST
    }

    if (!is.null(k)) {
      checkmate::assert_string(k, null.ok = TRUE)
      checkmate::assert_choice(k, choices = valid_ids)

      if (! k %in% valid_ids) {
	stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	  k, paste0(valid_ids, collapse=", ")))
      }

      return(cache$identifiers_list[[k]])
    } else {
      missing_ids <- setdiff(names(IDENTIFIERS_LIST), names(cache$identifiers_list))	
      if (length(missing_ids) != 0L) {	
        sapply(missing_ids, get_identifier) 	
      }

      return(cache$identifiers_list)
    }
  }

  .set_id <- function(k, v) {
    checkmate::assert_string(k, null.ok = TRUE)
    checkmate::assert_choice(k, choices = valid_ids)

    if (! k %in% valid_ids) {
      stop(sprintf("'%s' is not one of the valid identifiers: '%s'", 
	k, paste0(valid_ids, collapse=", ")))
    }
    cache$identifiers_list[[k]] <- v
  }


  .get_header <- function(k = NULL) {
    ## The following .getHeadersList() call is inside the function .get_header
    ## to avoid cyclical dependencies and collation order problems.

    if (length(cache$headers_list) == 0L) {
      cache$headers_list <- .getHeadersList()
    }

    if (!is.null(k)) {
      valid_headers <- names(cache$headers_list) 

      checkmate::assert_string(k, null.ok = TRUE)
      checkmate::assert_choice(k, choices = valid_headers)

      if (! k %in% valid_headers) {
	stop(sprintf("'%s' is not one of the valid headers: '%s'", 
	  k, paste0(valid_headers, collapse=", ")))
      }

      return(cache$headers_list[[k]])
    } else {
      return(cache$headers_list)
    }
  }

  .reset_ids <- function() {
    cache$identifiers_list <- list()
  }

  .reset_headers <- function() {
    cache$headers_list <- list()
  }

  list(get_id=.get_id, set_id=.set_id, get_head=.get_header, reset_ids=.reset_ids, reset_heads=.reset_headers)
}

global_cache <- make_global_cache()
