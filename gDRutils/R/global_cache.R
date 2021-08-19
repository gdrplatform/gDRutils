## This global cache maintains a cache of identifiers, 
## headers, and their respective values. 

global_cache <- new.env(parent = emptyenv())
global_cache$identifiers_list <- list() 
global_cache$headers_list <- list()


#############
# Identifiers
#############

#' @keywords internal
.get_id <- function(k = NULL) {
  if (length(global_cache$identifiers_list) == 0L) {
    global_cache$identifiers_list <- IDENTIFIERS_LIST
  }
  
  valid_ids <- names(IDENTIFIERS_LIST)
  if (!is.null(k)) {
    checkmate::assert_string(k, null.ok = TRUE)
    checkmate::assert_choice(k, choices = valid_ids)
    
    return(global_cache$identifiers_list[[k]])
  } else {
    missing_ids <-
      setdiff(valid_ids, names(global_cache$identifiers_list))
    if (length(missing_ids) != 0L) {
      global_cache$identifiers_list[missing_ids] <-
        IDENTIFIERS_LIST[missing_ids]
    }
    return(global_cache$identifiers_list)
  }
}


#' @keywords internal
.set_id <- function(k, v) {
  valid_ids <- names(IDENTIFIERS_LIST)

  checkmate::assert_string(k, null.ok = TRUE)
  checkmate::assert_choice(k, choices = valid_ids)

  global_cache$identifiers_list[[k]] <- v
  invisible(NULL)
}


#' @keywords internal
.reset_ids <- function() {
  global_cache$identifiers_list <- list()
  invisible(NULL)
}


##########
# Headers
##########

#' @keywords internal
.get_header <- function(k = NULL) {
  ## The following .getHeadersList() call is inside the function .get_header
  ## to avoid cyclical dependencies and collation order problems.

  if (length(global_cache$headers_list) == 0L) {
    global_cache$headers_list <- .getHeadersList()
  }

  if (!is.null(k)) {
    valid_headers <- names(global_cache$headers_list) 

    checkmate::assert_string(k, null.ok = TRUE)
    checkmate::assert_choice(k, choices = valid_headers)

    return(global_cache$headers_list[[k]])
  } else {
    return(global_cache$headers_list)
  }
}

#' @keywords internal
.reset_headers <- function() {
  global_cache$headers_list <- list()
  invisible(NULL)
}
