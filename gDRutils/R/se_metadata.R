#' Get and set metadata for parameters on a SummarizedExperiment object.
#'
#' Set metadata for the fitting parameters that define the Metrics assay in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to add fit parameter metadata.
#' @param value named list of metadata for fit parameters. 
#' @param key_type string of a specific key type (i.e. 'nested_keys', 'trt', 'masked_tag', etc.). 
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

#' @export
#' @rdname SE_metadata
#'
set_SE_fit_parameters <- function(se, value) {
  .set_SE_metadata(se, name = "fit_parameters", value)
}


#' @rdname SE_metadata
#' @export
#'
set_SE_processing_metadata <- function(se, value) {
  .set_SE_metadata(se, name = ".internal", value)
}


#' @rdname SE_metadata
#' @export
#'
set_SE_keys <- function(se, value) {
  .set_SE_metadata(se, name = "Keys", value)
}


#' @rdname SE_metadata
#' @export
#'
set_SE_experiment_metadata <- function(se, value) {
  .set_SE_metadata(se, name = "experiment_metadata", value)
}

############
# Getters
############

#' @rdname SE_metadata
#' @export
#'
get_SE_fit_parameters <- function(se) {
  .get_SE_metadata(se, name = "fit_parameters")
}


#' @rdname SE_metadata
#' @export
#'
get_SE_processing_metadata <- function(se) {
  .get_SE_metadata(se, name = ".internal")
}


#' @rdname SE_metadata
#' @export
#'
get_SE_experiment_metadata <- function(se) {
  .get_SE_metadata(se, name = "experiment_metadata")
}


#' @rdname SE_metadata
#' @export
#'
get_SE_keys <- function(se, key_type = NULL) {
  .get_SE_metadata(se, name = "Keys", subname = key_type)
}


##############
# Identifiers
##############

#' @rdname identifiers
#' @export
#'
get_SE_identifiers <- function(se, id_type = NULL) {
  ## `strict = FALSE` is present for backwards compatibility.
  out <-
    .get_SE_metadata(se,
                     name = "identifiers",
                     subname = id_type,
                     strict = FALSE)
  
  # throw error if an invalid identifier is provided with 'id_type'
  if (!is.null(id_type) && !id_type %in% names(get_identifier())) {
    stop(sprintf("Error: id_type:'%s' is an invalid identifier",
                 id_type))
  }
  
  # throw error if invalid identifier(s) is/are provided with metadata(se)
  invalid_idfs <- setdiff(names(out), names(get_identifier()))
  if (length(invalid_idfs) > 0) {
    stop(sprintf(
      "Error: metadata(se) contains invalid identifier(s): '%s'",
      toString(invalid_idfs)
    ))
  }
  
  out
}


#' @rdname identifiers
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
.set_SE_metadata <- function(se, name, value) {
  if (!is.null(.get_SE_metadata(se, name))) {
    warning(sprintf("overwriting existing metadata entry: '%s'", name))
  }
  S4Vectors::metadata(se)[[name]] <- value
  se
}
