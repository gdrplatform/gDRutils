#' Set metadata for parameters used in fitting in SummarizedExperiment object.
#'
#' Set metadata for the fitting parameters that define the Metrics assay in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to add fit parameter metadata.
#' @param value named list of metadata for fit parameters. 
#'
#' @return \linkS4class{se} with added metadata.
#'
#' @export
#'
set_SE_fit_parameters <- function(se, value) {
  .set_SE_metadata(se, name = "fit_parameters", value)
}


#' Set metadata for keys in SummarizedExperiment object.
#'
#' Set metadata for keys in SummarizedExperiment object metadata.
#'
#' @param se a \linkS4class{SummarizedExperiment} object for which to add key metadata.
#' @param value named list of metadata for keys. 
#' Names of list should represent key types and list values should contain key type values.
#'
#' @return \linkS4class{se} with added metadata.
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
#' @return \linkS4class{se} with added metadata.
#'
#' @export
#'
set_SE_experiment_metadata <- function(se, value) {
  .set_SE_metadata(se, name = "experiment_metadata", value)
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
#' @param key_type string of a specific key type (i.e. 'nested_keys', 'Trt', 'masked_tag', etc.). 
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
.get_SE_metadata <- function(se, name, subname = NULL) {
  v <- S4Vectors::metadata(se)[[name]]
  if (!is.null(subname)) {
    if (!subname %in% names(v)) {
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
