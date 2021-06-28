#' Set metadata for parameters used in fitting in SummarizedExperiment object.
#'
#' Set metadata for the fitting parameters that define the Metrics assay in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to add fit parameter metadata.
#' @param value named list of metadata for fit parameters. 
#'
#' @return \code{se} with added metadata.
#'
#' @export
#'
set_SE_fit_parameters <- function(se, value) {
  .set_SE_metadata(se, name = "fit_parameters", value)
}

#' Set metadata for the processing info in SummarizedExperiment object.
#'
#' Set metadata for the processing info that defines 
#' the date_processed and packages versions in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which 
#' to add processing info metadata.
#' @param value named list of metadata for processing info. 
#'
#' @return \code{se} with added metadata.
#'
#' @export
#'
set_SE_processing_metadata <- function(se, value) {
  .set_SE_metadata(se, name = ".internal", value)
}


#' Set metadata for keys in SummarizedExperiment object.
#'
#' Set metadata for keys in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to add key metadata.
#' @param value named list of metadata for keys. 
#' Names of list should represent key types and list values should contain key type values.
#'
#' @return \code{se} with added metadata.
#'
#' @export
#'
set_SE_keys <- function(se, value) {
  .set_SE_metadata(se, name = "Keys", value)
}


#' Set experiment-level metadata for a SummarizedExperiment object.
#'
#' Set experiment-level metadata in the metadata slot of a \linkS4class{SummarizedExperiment} object.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to add experiment-level metadata.
#' @param value named list of metadata for the \code{se} object. 
#'
#' @return \code{se} with added metadata.
#'
#' @export
#'
set_SE_experiment_metadata <- function(se, value) {
  .set_SE_metadata(se, name = "experiment_metadata", value)
}


#' Get metadata for the raw data identifier mappings used to create a SummarizedExperiment object.
#'
#' Get metadata for the identifiers used to construct the metadata and assay data of a SummarizedExperiment.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to get identifiers.
#' @param id_type string of identifier type to retrieve.
#' Defaults to \code{NULL}.
#'
#' @return named list of identifiers used during \code{create_SE} operations.
#'
#' @export
#'
get_SE_identifiers <- function(se, id_type = NULL) {
  ## `strict = FALSE` is present for backwards compatibility.
  ## If the identifier does not exist on the se object, we fetch it from the environment.
  value <- .get_SE_metadata(se, name = "identifiers", subname = id_type, strict = FALSE)
  if (is.null(value)) {
    value <- get_identifier(id_type)
  }
  value
}


#' Get metadata for fit parameters metadata in SummarizedExperiment object.
#'
#' Get metadata for fit parameters used to construct the Metrics assay in a SummarizedExperiment object.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to get fit parameters.
#'
#' @return named list of fitting parameters.
#'
#' @export
#'
get_SE_fit_parameters <- function(se) {
  .get_SE_metadata(se, name = "fit_parameters")
}

#' Set metadata for the processing info in SummarizedExperiment object.
#'
#' Set metadata for the processing info that defines 
#' the date_processed and packages versions in SummarizedExperiment object metadata.
#' 
#' @param se a \linkS4class{SummarizedExperiment} object for which to get processing info.
#'
#' @return named list of processing info.
#'
#' @export
#'
get_SE_processing_metadata <- function(se) {
  .get_SE_metadata(se, name = ".internal")
}


#' Get metadata for experiment metadata in SummarizedExperiment object.
#'
#' Get metadata for experiment metadata in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to get experiment metadata.
#'
#' @return experiment metadata.
#'
#' @export
#'
get_SE_experiment_metadata <- function(se) {
  .get_SE_metadata(se, name = "experiment_metadata")
}


#' Get metadata for keys in SummarizedExperiment object.
#'
#' Get metadata for keys in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to get metadata relating to keys.
#' @param key_type string of a specific key type (i.e. 'nested_keys', 'trt', 'masked_tag', etc.). 
#'
#' @export
#'
get_SE_keys <- function(se, key_type = NULL) {
  .get_SE_metadata(se, name = "Keys", subname = key_type)
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
