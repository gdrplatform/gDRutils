#' @export
.clean_key_inputs <- function(keys, cols) {
  dropped <- setdiff(keys, cols)
  if (length(dropped) != 0L) {
    warning(sprintf("ignoring input keys: '%s' which are not present in data.frame",
      paste0(dropped, collapse = ", ")))
  }
  intersect(keys, cols)
}


#' @noRd
#' @keywords internal
assert_equal_input_len <- function(outlier, ...) {
  first <- list(...)[[1]]
  h <- all(sapply(list(...), length) == length(first))
  if (!h) {
    stop("unequal length objects provided as input")
  }

  contains_length_one <- length(first) == 1L || length(outlier) == 1L
  if (length(first) != length(outlier) && !contains_length_one) {
    stop("unequal lengths detected, either the fit parameters must be length one, or the tested value")
  }

  invisible(NULL)
}

#' shorten normalization type
#' 
#' @param x string with normalization type
#' 
#' @export
shorten_normalization_type_name <- function(x) {
  checkmate::assert_choice(x, c("RelativeViability", "GRvalue"))
  dict <- c("RelativeViability" = "RV", "GRvalue" = "GR")
  dict[[x]]
}

#' extend abbreviated normalization type
#' 
#' @param x string with normalization type
#' 
#' @export
extend_normalization_type_name <- function(x) {
  checkmate::assert_choice(x, c("RV", "GR"))
  dict <- c("RV" = "RelativeViability", "GR" = "GRvalue")
  dict[[x]]
}

#' assert choices
#' 
#' @param x charvec expected subset
#' @param choices charvec reference set
#' @param ... Additional arguments to pass to \code{checkmate::test_choice}
#' @export
assert_choices <- function(x, choices, ...) {
  out <- vapply(x, function(y) {
    checkmate::test_choice(y, choices, ...)
  }, FUN.VALUE = logical(1))

  if (!all(out)) {
    msg <-
      sprintf(
        "Assertion on '%s' failed. Must be element(s) of {'%s'} set.",
        toString(x[!out]),
        toString(choices)
      )
    stop(msg)
  }
}

#' Lapply through all the experiments in MultiAssayExperiment object
#' 
#' @param mae MultiAssayExperiment object
#' @param FUN function that should be applied on each experiment of MultiAssayExperiment object
#' @param unify logical indicating if the output should be a unlisted object of unique
#' values across all the experiments 
#' @param ... Additional args to be passed to teh \code{FUN}.
#' @export
#' 
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
MAEpply <- function(mae, FUN, unify = FALSE, ...) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  checkmate::assert_function(FUN)
  checkmate::assert_flag(unify)
  experiments <- as.list(MultiAssayExperiment::experiments(mae))
  out <- lapply(experiments, FUN, ...)
  if (unify) {
    if (all(vapply(out, is.atomic, logical(1)))) {
      unlist(out, use.names = FALSE)
    } else {
      plyr::rbind.fill(lapply(out, data.frame))
    }
    
  } else {
    out
  }
}

#' is_mae_empty
#' 
#' check if all mae experiments are empty
#' @param mae MultiAssayExperiment object
#' @export
#' 
#' @return logical
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
is_mae_empty <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  
  all(MAEpply(mae, is_exp_empty, unify = TRUE))
}

#' is_any_exp_empty
#' 
#' check if any experiment is empty
#' @param mae MultiAssayExperiment object
#' @export
#' 
#' @return logical
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
is_any_exp_empty <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  
  any(MAEpply(mae, is_exp_empty, unify = TRUE))
}

#' is_exp_empty
#' 
#' check if experiment (SE) is empty
#' @param exp \linkS4class{SummarizedExperiment} object.
#' @export
#' 
#' @return logical
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
is_exp_empty <- function(exp) {
  checkmate::assert_class(exp, "SummarizedExperiment")
  
  names <- SummarizedExperiment::assayNames(exp)
  dt <- `if`(
    is.null(names),
    data.table::data.table(),
    convert_se_assay_to_dt(exp, names[[1]])
  )
  
  any(
    nrow(SummarizedExperiment::assay(exp)) == 0,
    nrow(dt) == 0
  )
}

#' get_non_empty_assays
#' 
#' get non empty assays
#' @param mae MultiAssayExperiment object
#' @export
#' 
#' @return charvec with non-empty experiments
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
get_non_empty_assays <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  
  ne_info <- MAEpply(mae, is_exp_empty) == FALSE
  names(ne_info[ne_info == TRUE])
}

#' mcolData
#'
#' get colData of all experiments
#' @param mae MultiAssayExperiment object
#' @export
#'
#' @return tibble with all-experiments colData
#'
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
mcolData <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  
  MAEpply(mae, SummarizedExperiment::colData, unify = TRUE)
}

#' mrowData
#'
#' get rowData of all experiments
#' @param mae MultiAssayExperiment object
#' @export
#'
#' @return tibble with all-experiments rowData
#'
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
mrowData <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  
  MAEpply(mae, SummarizedExperiment::rowData, unify = TRUE)
}
