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
#' @return shortened string representing the normalization type
#' @export
shorten_normalization_type_name <- function(x) {
  checkmate::assert_choice(x, c("RelativeViability", "GRvalue"))
  dict <- c("RelativeViability" = "RV", "GRvalue" = "GR")
  unname(dict[match(x, names(dict))])
}

#' extend abbreviated normalization type
#' 
#' @param x string with normalization type
#' 
#' @export
extend_normalization_type_name <- function(x) {
  checkmate::assert_choice(x, c("RV", "GR"))
  dict <- c("RV" = "RelativeViability", "GR" = "GRvalue")
  unname(dict[match(x, names(dict))])
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

#' Lapply or bplapply.
#'
#' @param x Vector (atomic or list) or an ‘expression’ object.
#' Other objects (including classed objects) will be coerced by
#' ‘base::as.list’.
#' @param FUN A user-defined function.
#' @param parallelize Logical indicating whether or not to parallelize the computation.
#' Defaults to \code{TRUE}.
#' @param ... optional argument passed to 
#' \link[BiocParallel]{bplapply} if \code{parallelize == TRUE},
#' else to \link[base]{lapply}.
#'
#' @return List containing output of \code{FUN} applied to every element in \code{x}.
#' @export
loop <- function(x, FUN, parallelize = TRUE, ...) {
  if (parallelize) {
    BiocParallel::bplapply(x, FUN, ...)
  } else {
    lapply(x, FUN, ...)
  }
}

#' Apply a function to every element of a bumpy matrix.
#'
#' Apply a user-specified function to every element of a bumpy matrix.
#'
#' @param se A \code{SummarizedExperiment} object with bumpy matrices.
#' @param FUN A function that will be applied to each element of the matrix in assay \code{req_assay_name}.
#' Output of the function must return a data.frame.
#' @param req_assay_name String of the assay name in the \code{se} that the \code{FUN} will act on.
#' @param out_assay_name String of the assay name that will contain the results of the applied function.
#' @param parallelize Logical indicating whether or not to parallelize the computation.
#'
#' @return The original \code{se} object with a new assay, \code{out_assay_name}.
#' @export
apply_bumpy_function <- function(se, FUN, req_assay_name, out_assay_name, parallelize = FALSE) {
  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(req_assay_name)
  checkmate::assert_string(out_assay_name)
  gDRutils::validate_se_assay_name(se, req_assay_name)

  asy <- SummarizedExperiment::assay(se, req_assay_name)
  checkmate::assert_class(asy, "BumpyDataFrameMatrix")
  df <- BumpyMatrix::unsplitAsDataFrame(asy, row.field = "row", column.field = "column")
  iterator <- unique(df[, c("column", "row")])
  out <- loop(seq_len(nrow(iterator)), FUN = function(elem) {
    x <- iterator[elem, ]
    i <- x[["row"]]
    j <- x[["column"]]
    elem_df <- asy[i, j][[1]]

    store <- FUN(elem_df)
    if (is(store, "data.frame") || is(store, "DFrame")) {
      if (nrow(store) != 0L) {
        store$row <- i
        store$column <- j
        store
      } else {
        NULL 
      }
    } else {
      stop("only data.frame objects supported as return values from FUN for now")
    }
  }, parallelize = parallelize)

  out <- S4Vectors::DataFrame(do.call("rbind", out))

  out_assay <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c("row", "column")],
    row = out$row,
    col = out$column)
  SummarizedExperiment::assays(se)[[out_assay_name]] <- out_assay
  se
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
