#' Get and set metadata for parameters on a SummarizedExperiment object.
#'
#' Set metadata for the fitting parameters that define the Metrics assay in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to add fit parameter metadata.
#' @param value named list of metadata for fit parameters.
#' @param key_type string of a specific key type (i.e. 'nested_keys', etc.).
#' @param id_type string of a specific id type (i.e. 'duration', 'cellline_name', etc.).
#' @param simplify Boolean indicating whether output should be simplified.
#' @param append Boolean indicating whether to append the new metadata value to the existing entry. 
#'
#' @details
#' For \code{*et_SE_processing_metadata}, get/set metadata for the processing info that defines
#' the date_processed and packages versions in SummarizedExperiment object metadata.
#' For \code{*et_SE_fit_parameters}, get/set metadata for fit parameters
#' used to construct the Metrics assay in a SummarizedExperiment object.
#'
#' @return \code{se} with added metadata.
#' @name SE_metadata
#'
NULL


############
# Setters
############

#' @keywords SE_operators
#' @export
#' @rdname SE_metadata
#'
set_SE_fit_parameters <- function(se, value) {
  .set_SE_metadata(se, name = "fit_parameters", value)
}


#' @keywords SE_operators
#' @rdname SE_metadata
#' @export
#'
set_SE_processing_metadata <- function(se, value) {
  .set_SE_metadata(se, name = ".internal", value)
}


#' @keywords SE_operators
#' @rdname SE_metadata
#' @export
#'
set_SE_keys <- function(se, value) {
  .set_SE_metadata(se, name = "Keys", value)
}


#' @keywords SE_operators
#' @rdname SE_metadata
#' @export
#'
set_SE_experiment_metadata <- function(se, value, append = TRUE) {
  .set_SE_metadata(se, name = "experiment_metadata", value, append = append)
}

#' @keywords SE_operators
#' @rdname SE_metadata
#' @export
#'
set_SE_experiment_raw_data <- function(se, value) {
  .set_SE_metadata(se, name = "experiment_raw_data", value)
}



############
# Getters
############
#' @rdname SE_metadata
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs") 
#' se <- mae[[1]]
#' get_SE_fit_parameters(se)
#' 
#' @keywords SE_operators
#' @export
get_SE_fit_parameters <- function(se) {
  .get_SE_metadata(se, name = "fit_parameters")
}


#' @rdname SE_metadata
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs") 
#' se <- mae[[1]]
#' meta <- get_SE_processing_metadata(se)
#' 
#' @keywords SE_operators
#' @export
get_SE_processing_metadata <- function(se) {
  .get_SE_metadata(se, name = ".internal")
}

#' @rdname SE_metadata
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs") 
#' se <- mae[[1]]
#' get_SE_experiment_raw_data(se)
#' 
#' @keywords SE_operators
#' @export
get_SE_experiment_raw_data <- function(se) {
  .get_SE_metadata(se, name = "experiment_raw_data")
}

#' @rdname SE_metadata
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs") 
#' se <- mae[[1]]
#' get_SE_experiment_metadata(se)
#' 
#' @keywords SE_operators
#' @export
get_SE_experiment_metadata <- function(se) {
  .get_SE_metadata(se, name = "experiment_metadata")
}

#' @keywords SE_operators
#' @rdname SE_metadata
#' @export
#'
get_SE_keys <- function(se, key_type = NULL) {
  .get_SE_metadata(se, name = "Keys", subname = key_type)
}


##############
# Identifiers
##############

#' @rdname SE_metadata
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs") 
#' se <- mae[[1]]
#' get_SE_identifiers(se)
#' 
#' @keywords SE_operators
#' @export
get_SE_identifiers <- function(se, id_type = NULL, simplify = TRUE) {
  ## `strict = FALSE` is present for backwards compatibility.
  if (simplify) {
    if (length(id_type) > 1L) {
      stop("more than one identifier found, please use: simplify = FALSE")
    } else {
      out <- .get_SE_identifier(se, id_type)
    }
  } else {
    id_vector <- Vectorize(function(i) {
      unlist(unname(lapply(i, function(x) .get_SE_identifier(se, x))))
    },
    SIMPLIFY = FALSE)
    out <- id_vector(id_type)
  }
  
  out
}

#' get_MAE_identifiers
#'
#' get the identifiers of all SE's in the MAE
#' @param mae MultiAssayExperiment
#'
#' @return named list with identifiers for each SE
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' get_MAE_identifiers(mae)
#' 
#' @keywords SE_operators
#' @export
get_MAE_identifiers <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  
  MAEpply(mae, get_SE_identifiers)
}


#' @keywords internal
#' @noRd
.get_SE_identifier <- function(se, id_type) {
  checkmate::assert_choice(id_type, choices = names(IDENTIFIERS_LIST), null.ok = TRUE)
  out <-
    .get_SE_metadata(se,
                     name = "identifiers",
                     subname = id_type,
                     strict = FALSE)
  # For backwards compatibility, get identifiers from environment if not found on `se`.
  if (is.null(out)) {
    warning(sprintf("'se' was passed, but identifier '%s' not found on se's identifiers", id_type))
    out <- get_env_identifiers(id_type, simplify = TRUE)
  }
  out
}


#' @keywords SE_operators
#' @rdname SE_metadata
#' @export
#'
set_SE_identifiers <- function(se, value) {
  .set_SE_metadata(se, name = "identifiers", value)
}

###############
# Internals
###############

#' The primary purpose of this function is to allow other functions to create exposed getter functions.
#' @noRd
.get_SE_metadata <- function(se, name, subname = NULL, strict = TRUE) {
  v <- S4Vectors::metadata(se)[[name]]
  if (!is.null(subname)) {
    if (!subname %in% names(v) && strict) {
      stop(sprintf("'%s' is not one of valid subname(s): '%s'", subname, paste0(names(v), collapse = ", ")))
    }
    v <- v[[subname]]
  }
  v
}


#' The primary purpose of this function is to allow other functions to create exposed setter functions.
#' @noRd
.set_SE_metadata <- function(se, name, value, append = FALSE) {
  current_metadata <- .get_SE_metadata(se, name)
  
  if (is.null(current_metadata) || S4Vectors::isEmpty(current_metadata)) {
    current_metadata <- NULL
  }
  
  if (!is.null(current_metadata)) {
    if (append) {
      S4Vectors::metadata(se)[[name]] <- c(current_metadata, value)
    } else {
      warning(sprintf("overwriting existing metadata entry: '%s'", name))
      S4Vectors::metadata(se)[[name]] <- value
    }
  } else {
    S4Vectors::metadata(se)[[name]] <- value
  }
  
  se
}
