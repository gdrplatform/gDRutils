#' @keywords package_utils
#' @export
.clean_key_inputs <- function(keys, cols) {
  dropped <- setdiff(keys, cols)
  if (length(dropped) != 0L) {
    warning(sprintf("ignoring input keys: '%s' which are not present in data.table",
      paste0(dropped, collapse = ", ")))
  }
  intersect(keys, cols)
}


#' @noRd
#' @keywords internal
assert_equal_input_len <- function(outlier, ...) {
  first <- list(...)[[1]]
  h <- all(vapply(list(...), length, FUN.VALUE = 1) == length(first))
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
#' 
#' @examples 
#' shorten_normalization_type_name("GRvalue")
#' 
#' @keywords package_utils
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
#' @return string
#' 
#' @examples 
#' extend_normalization_type_name("GR")
#' 
#' @keywords package_utils
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
#' 
#' @return \code{NULL}
#' 
#' @examples 
#' assert_choices("x", c("x","y"))
#' 
#' @keywords package_utils
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
#' @keywords package_utils
#' @export
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
#' 
#' @return list or vector depends on unify param
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' MAEpply(mae, SummarizedExperiment::assayNames)
#' 
#' @keywords package_utils
#' @export
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
      data.table::rbindlist(lapply(out, data.table::as.data.table), fill = TRUE)
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
#' 
#' @examples 
#' loop(list(c(1,2), c(2,3)), sum, parallelize = FALSE)
#' 
#' @keywords package_utils
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
#' Output of the function must return a data.table.
#' @param req_assay_name String of the assay name in the \code{se} that the \code{FUN} will act on.
#' @param out_assay_name String of the assay name that will contain the results of the applied function.
#' @param parallelize Logical indicating whether or not to parallelize the computation.
#' @param ... Additional args to be passed to teh \code{FUN}.
#' @return The original \code{se} object with a new assay, \code{out_assay_name}.
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' se <- mae[[1]]
#' FUN <- function(x) {
#'   data.table::data.table(Concentration = x$Concentration, CorrectedReadout = x$CorrectedReadout)
#' } 
#' apply_bumpy_function(
#'   se, 
#'   FUN = FUN, 
#'   req_assay_name = "RawTreated", 
#'   out_assay_name = "CorrectedReadout"
#' )
#' 
#' @keywords package_utils
#' @export
apply_bumpy_function <- function(se,
                                 FUN,
                                 req_assay_name,
                                 out_assay_name,
                                 parallelize = FALSE,
                                 ...) {
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
    store <- FUN(elem_df, ...)
    if (is(store, "data.table") || is(store, "DFrame")) {
      if (nrow(store) != 0L) {
        store$row <- i
        store$column <- j
        store
      } else {
        NULL
      }
    } else {
      stop("only data.table objects supported as return values from FUN for now")
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
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
#' 
#' @return logical
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' is_mae_empty(mae)
#' 
#' @keywords package_utils
#' @export
is_mae_empty <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")

  all(MAEpply(mae, is_exp_empty, unify = TRUE))
}

#' is_any_exp_empty
#'
#' check if any experiment is empty
#' @param mae MultiAssayExperiment object
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
#' 
#' @return logical
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' is_any_exp_empty(mae)
#' 
#' @keywords package_utils
#' @export
is_any_exp_empty <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")

  any(MAEpply(mae, is_exp_empty, unify = TRUE))
}

#' is_exp_empty
#'
#' check if experiment (SE) is empty
#' @param exp \linkS4class{SummarizedExperiment} object.
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
#' 
#' @return logical
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' se <- mae[[1]]
#' is_exp_empty(se)
#' 
#' @keywords package_utils
#' @export
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
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
#' 
#' @return charvec with non-empty experiments
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' get_non_empty_assays(mae)
#' 
#' @keywords package_utils
#' @export
get_non_empty_assays <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")

  ne_info <- MAEpply(mae, is_exp_empty) == FALSE
  names(ne_info[ne_info == TRUE])
}

#' mcolData
#'
#' get colData of all experiments
#' @param mae MultiAssayExperiment object
#' 
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
#' 
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs")
#' mcolData(mae)
#'
#' @return data.table with all-experiments colData
#'
#' @keywords package_utils
#' @export
mcolData <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")

  MAEpply(mae, SummarizedExperiment::colData, unify = TRUE)
}

#' mrowData
#'
#' get rowData of all experiments
#' @param mae MultiAssayExperiment object
#' @keywords package_utils
#' @export
#'
#' @return data.table with all-experiments rowData
#'
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small.qs") 
#' mrowData(mae)
#'
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
mrowData <- function(mae) {
  checkmate::assert_class(mae, "MultiAssayExperiment")

  MAEpply(mae, SummarizedExperiment::rowData, unify = TRUE)
}

#' Get synthetic data from gDRtestData package
#'
#' @param qs qs filename
#' 
#' @keywords package_utils
#' @export
#' 
#' @examples 
#' get_synthetic_data("finalMAE_small.qs") 
#'
#' @return loaded data
#' 
get_synthetic_data <- function(qs) {
  # check if prefix exist, if not add one
  if (!grepl("finalMAE", qs)) {
    qs <- paste("finalMAE", qs, sep = "_")
  }
  # check if extension exist and is supported -`qs`, if not add one or replace
  if (!grepl(".qs$", qs)) {
    if (grepl(".RDS$", qs)) {
      qs <- gsub(".RDS", ".qs", qs)
    } else {
      qs <- paste0(qs, ".qs")
    }
  }
  qs::qread(system.file("testdata", qs, package = "gDRtestData"))
}


#' Geometric mean
#' 
#' Auxiliary function for calculating geometric mean with possibility to handle -Inf
#' 
#' @param x numeric vector
#' @param fixed flag should be add fix for -Inf 
#' @param maxlog10Concentration numeric value needed to calculate minimal value
#' 
#' @return numeric vector
#' 
#' @examples 
#' geometric_mean(c(2, 8))
#' 
#' @keywords package_utils
#' @export
#' 
#' @keywords internal
geometric_mean <- function(x, fixed = TRUE, maxlog10Concentration = 1) {
  checkmate::assert_numeric(x)
  checkmate::assert_flag(fixed)
  checkmate::assert_numeric(maxlog10Concentration)
  
  if (fixed) {
    x <- pmax(
      10 ^ maxlog10Concentration / 1e6,
      pmin(5 * 10 ^ maxlog10Concentration, x)
    )
  }
  exp(mean(log(x)))
}

#' Average biological replicates.
#'
#' Average biological replicates on the data table side. 
#'
#' @param dt data.table with Metric data
#' @param var String representing additional metadata of replicates
#' @param pidfs list of prettified identifiers
#' @param fixed Flag should be add fix for -Inf in geometric mean.
#' @param geometric_average_fields Character vector of column names in \code{dt} 
#' to take the geometric average of.
#' 
#' @examples
#' dt <- data.table::data.table(a = c(1:10, 1),
#' b = c(rep("drugA", 10), rep("drugB", 1)))
#' average_biological_replicates_dt(dt, var = "a")
#' 
#' @return data.table without replicates
#' @keywords package_utils
#' @export
average_biological_replicates_dt <- function(
    dt,
    var,
    pidfs = get_prettified_identifiers(),
    fixed = TRUE,
    geometric_average_fields = get_header("metric_average_fields")$geometric_mean) {
  data <- data.table::copy(dt)
  average_fields <- setdiff(names(Filter(is.numeric, data)), c(unlist(pidfs),
                                                               var,
                                                               prettify_flat_metrics(get_header("iso_position"),
                                                                                     human_readable = TRUE)))
  geometric_average_fields <- intersect(geometric_average_fields, names(dt))
  group_by <- setdiff(names(data),
                      c(average_fields, var, prettify_flat_metrics(get_header("id"), human_readable = TRUE)))
  group_by <- grep("Fit Type", group_by, invert = TRUE, value = TRUE)
  data <- data[, (var) := NULL][, 
                                (average_fields) := lapply(.SD, mean, na.rm = TRUE), 
                                .SDcols = average_fields, 
                                by = group_by][,
                                               (geometric_average_fields) := lapply(.SD, FUN = function(x) {
                                                 geometric_mean(x, fixed = fixed)
                                               }), 
                                               .SDcols = geometric_average_fields, 
                                               by = group_by]
  unique(data, by = group_by)
}

#' Helper function to find duplicated rows
#'
#' @param x data frame
#' @param col_names character vector, columns in which duplication are searched for
#' @return integer vector
#' @examples
#' dt <- data.table::data.table(a = c(1, 2, 3), b = c(3, 2, 2))
#' get_duplicated_rows(dt, "b")
#' @keywords package_utils
#' @export
get_duplicated_rows <- function(x, col_names = NULL) {
  checkmate::assertMultiClass(x, c("data.table", "DataFrame"))
  checkmate::assert_true(all(col_names %in% colnames(x)))
  
  if (!is.null(col_names)) {
    x <- subset(x, select = col_names)
  }
  which(duplicated(x) | duplicated(x, fromLast = TRUE))
}
