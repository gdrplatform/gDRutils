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
  validate_se_assay_name(se, req_assay_name)

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

  out <- S4Vectors::DataFrame(do.call(rbind, out))

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

#' Average biological replicates on the data table side. 
#'
#' @param dt data.table with Metric data
#' @param var String representing additional metadata of replicates
#' @param prettified Flag indicating if the provided identifiers in the dt are prettified
#' @param fixed Flag indicating whether to add a fix for -Inf in the geometric mean.
#' @param geometric_average_fields Character vector of column names in \code{dt} 
#' to take the geometric average of.
#' @param fit_type_average_fields Character vector of column names in \code{dt} 
#' that should be treated as a column with fit type data
#' @param blacklisted_fields Character vector of column names in \code{dt} 
#' that should be skipped in averaging
#' @param add_sd Flag indicating whether to add standard deviation and count columns.
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
    prettified = FALSE,
    fixed = TRUE,
    geometric_average_fields = get_header("metric_average_fields")$geometric_mean,
    fit_type_average_fields = get_header("metric_average_fields")$fit_type,
    blacklisted_fields = get_header("metric_average_fields")$blacklisted,
    add_sd = FALSE) {
  
  checkmate::assert_data_table(dt)
  checkmate::assert_string(var)
  checkmate::assert_flag(prettified)
  checkmate::assert_character(geometric_average_fields)
  checkmate::assert_character(fit_type_average_fields)
  checkmate::assert_flag(add_sd)
  
  data <- data.table::copy(dt)
  
  if (prettified) {
    pidfs <- get_prettified_identifiers()
    iso_cols <- prettify_flat_metrics(get_header("iso_position"), human_readable = TRUE)
    id_cols <- prettify_flat_metrics(get_header("id"), human_readable = TRUE)
  } else {
    pidfs <- get_env_identifiers()
    iso_cols <- get_header("iso_position")
    id_cols <- prettify_flat_metrics(get_header("id"))
  }
  
  average_fields <- setdiff(names(Filter(is.numeric, data)), c(unlist(pidfs), var, iso_cols))
  geometric_average_fields <- intersect(geometric_average_fields, names(dt))
  fit_type_average_fields <- intersect(fit_type_average_fields, names(dt))
  blacklisted_fields <- intersect(blacklisted_fields, names(dt))
  group_by <- setdiff(names(data), c(average_fields, var, id_cols, fit_type_average_fields, blacklisted_fields))
  
  if (add_sd) {
    # Calculate standard deviation for both average_fields and geometric_average_fields
    sd_fields <- paste0(average_fields, "_sd")
    geom_sd_fields <- paste0(geometric_average_fields, "_sd")
    
    data <- data[, (sd_fields) := lapply(.SD, calc_sd),
                 .SDcols = average_fields, by = group_by]
    data <- data[, (geom_sd_fields) := lapply(.SD, calc_sd),
                 .SDcols = geometric_average_fields, by = group_by]
    
    # Calculate count and add as a single column
    data <- data[, count := .N, by = group_by]
  }
  
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

#' Checks if \code{se} is combo dataset.
#'
#' @param se SummarizedExperiment
#' 
#' @examples
#' se <- get_synthetic_data("combo_matrix")[[1]]
#' is_combo_data(se)
#' se <- get_synthetic_data("combo_matrix")[[2]]
#' is_combo_data(se)
#' se <- get_synthetic_data("small")[[1]]
#' is_combo_data(se)
#'
#' @return logical
#' @keywords combination_data
#' 
#' @export
is_combo_data <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  
  all(get_combo_assay_names() %in% SummarizedExperiment::assayNames(se))
}

#' Has Single Codrug Data
#'
#' @param cols character vector with the columns of the input data
#' @param prettify_identifiers logical flag specifying if identifiers are expected to be prettified
#' @param codrug_identifiers character vector with identifiers for the codrug columns
#' 
#' @examples
#' has_single_codrug_data("Drug Name")
#' has_single_codrug_data(c("Drug Name", "Cell Lines"))
#' has_single_codrug_data(c("Drug Name 2", "Concentration 2"))
#' has_single_codrug_data(
#'   get_prettified_identifiers(
#'     c("concentration2", "drug_name2"), 
#'     simplify = FALSE
#'   )
#' )
#'
#' @keywords combination_data
#' @return logical flag
#' 
#' @export
has_single_codrug_data <-
  function(cols,
           prettify_identifiers = TRUE,
           codrug_identifiers = c("drug_name2", "concentration2")) {
    
    checkmate::assert_true(all(codrug_identifiers %in% names(get_env_identifiers(simplify = TRUE))))
    checkmate::assert_flag(prettify_identifiers)
    
    codrug_colnames <- if (prettify_identifiers) {
      get_prettified_identifiers(codrug_identifiers, simplify = FALSE)
    } else {
      unname(unlist(get_env_identifiers(codrug_identifiers, simplify = FALSE)))
    }
    checkmate::assert_character(cols, any.missing = FALSE)
    checkmate::assert_character(codrug_colnames, any.missing = FALSE)
    
    all(codrug_colnames %in% cols)
  }


#' Has Valid Codrug Data
#'
#' @param data data.table with input data
#' @param prettify_identifiers logical flag specifying if identifiers are expected to be prettified
#' @param codrug_name_identifier string with the identifiers for the codrug drug_name column
#' @param codrug_conc_identifier string with the identifiers for the codrug concentration column
#' 
#' @examples
#' dt <-
#'   data.table::data.table(
#'     "Drug Name" = letters[seq_len(3)],
#'     "Concentration" = seq_len(3),
#'     "Drug Name 2" = letters[4:6],
#'     "Concentration 2" = 4:6
#'   )
#' has_valid_codrug_data(dt)
#' 
#' dt$`Concentration 2` <- NULL
#' has_valid_codrug_data(dt)
#'
#' @keywords combination_data
#' @return logical flag
#' 
#' @export
has_valid_codrug_data <-
  function(data,
           prettify_identifiers = TRUE,
           codrug_name_identifier = "drug_name2",
           codrug_conc_identifier = "concentration2") {
    checkmate::assert_data_table(data)
    dcols <- colnames(data)
    checkmate::assert_flag(prettify_identifiers)
    checkmate::assert_string(codrug_name_identifier)
    checkmate::assert_string(codrug_conc_identifier)
    
    idfs <- if (prettify_identifiers) {
      get_prettified_identifiers(simplify = TRUE)
    } else {
      get_env_identifiers()
    }
    
    codrug_v <- c(codrug_name_identifier, codrug_conc_identifier)
    
    status <-
      # codrug data not present for drug_name and/or concentration data
      if (!has_single_codrug_data(dcols, prettify_identifiers, codrug_v)) {
        FALSE
      }  else {
        codrug_cols <- as.character(idfs[codrug_v])
        
        # codrug data not valid (for drug names and/or concentration data)
        if (all(data[[codrug_cols[1]]] %in% idfs[["untreated_tag"]]) ||
            all(is.na(data[[codrug_cols[2]]]))) {
          FALSE
        } else {
          TRUE
        }
      }
    status
  }

#' Remove Codrug Data
#'
#' @param data data.table with input data
#' @param prettify_identifiers logical flag specifying if identifiers are expected to be prettified
#' @param codrug_identifiers character vector with identifiers for the codrug columns
#' 
#' @examples
#' 
#' dt <-
#'   data.table::data.table(
#'     "Drug Name" = letters[seq_len(3)],
#'     "Concentration" = seq_len(3),
#'     "Drug Name 2" = letters[4:6],
#'     "Concentration 2" = 4:6
#'   )
#' dt
#' remove_codrug_data(dt)
#'
#' @keywords combination_data
#' @return data.table without combination columns
#' 
#' @export
remove_codrug_data <-
  function(data,
           prettify_identifiers = TRUE,
           codrug_identifiers = c("drug_name2", "concentration2")) {
    
    checkmate::assert_true(all(codrug_identifiers %in% names(get_env_identifiers())))
    checkmate::assert_data_table(data)
    checkmate::assert_flag(prettify_identifiers)
    
    codrug_colnames <- if (prettify_identifiers) {
      vapply(codrug_identifiers, function(x) get_prettified_identifiers(x), character(1))
    } else {
      vapply(codrug_identifiers, function(x) get_env_identifiers(x), character(1))
    }
    checkmate::assert_character(codrug_colnames, any.missing = FALSE)
    
    idx <- which(!colnames(data) %in% codrug_colnames)
    
    # support both: data.table and data.frame
    subset(data, select = idx)
  }

#' Identify and return additional variables in list of dt
#'
#' @param dt_list list of data.table or data.table containing additional variables
#' @param unique logical flag indicating if all variables should be returned 
#' or only those containing more than one unique value
#' @param prettified Flag indicating if the provided identifiers in the dt are prettified
#' 
#' @examples
#' dt <- data.table::data.table(
#'   Gnumber = seq_len(10), 
#'   Concentration = runif(10), 
#'   Ligand = c(rep(0.5, 5), rep(0, 5))
#' )
#' get_additional_variables(dt)
#'
#' @return vector of variable names with additional variables
#' 
#' @keywords combination_data
#' @export
get_additional_variables <- function(dt_list,
                                     unique = FALSE,
                                     prettified = FALSE) {
  
  
  if (data.table::is.data.table(dt_list)) {
    dt_list <- list(dt_list)
  }
  checkmate::assert_flag(unique)
  checkmate::assert_flag(prettified)
  
  if (prettified) {
    headers <- prettify_flat_metrics(unlist(get_header()), human_readable = TRUE)
    pidfs <- get_prettified_identifiers()
    idf2keep <- pidfs[c("drug3", "concentration3", "duration")]
    idfs <- setdiff(unique(c(headers, pidfs)), idf2keep)
  } else {
    headers <- unlist(get_header())
    pidfs <- get_env_identifiers()
    idf2keep <- pidfs[c("drug3", "concentration3", "duration")]
    idfs <- setdiff(unique(c(headers, pidfs)), idf2keep)
  }
  

  additional_perturbations <- unique(unlist(lapply(dt_list, function(x) {
    setdiff(sub(" \\(.*\\)$", "", names(x)), idfs)
  })))
  
  if (unique) {
    additional_perturbations
  } else {
    unlist(lapply(additional_perturbations, function(x) {
      if (any(vapply(dt_list, function(y) {
        length(unique(y[[x]])) > 1
      }, FUN.VALUE = logical(1)))) {
        x
      } else {
        NULL
      }
    }))
  }
}

#' Calculate Standard Deviation or Return Zero
#'
#' This function calculates the standard deviation of a numeric vector.
#' If the vector has a length of 1 and it is numeric, it returns 0.
#'
#' @param x A numeric vector.
#' @return The standard deviation of the vector if its length is greater than 1 or it is not numeric, otherwise 0.
#' @examples
#' calc_sd(c(1, 2, 3, 4, 5)) # Should return the standard deviation
#' calc_sd(c(1)) # Should return 0
#' calc_sd(numeric(0)) # Should return NA
#' calc_sd(c("a", "b", "c")) # Should return NA
#' @keywords package_utils
#' @export
calc_sd <- function(x) {
  if (length(x) == 1 && is.numeric(x) && !is.na(x)) {
    return(0)
  } else {
    return(stats::sd(x, na.rm = TRUE))
  }
}



#' safe wrapper of Sys.getenv()
#' 
#' So far the helper is needed to handle env vars containing `:` 
#' for which the backslash  is automatically added in some contexts
#' and R could not get the original value for these env vars.
#' 
#' @param x string with the name of the environemntal variable
#' @param ... additional params for Sys.getenev
#' @keywords package_utils
#' 
#' @examples 
#' get_env_var("HOME")
#
#' @export 
#' @return sanitized value of the env variable
get_env_var <- function(x, ...) {
  gsub("\\\\", "", Sys.getenv(x, ...))
}

#' Remove batch substring from drug id
#'
#' Gnumber, i.e. "G12345678" is currently the default format of drug_id. It's also used as a drug name in some cases.
#'
#' By default, Gnumber(s) followed by any character (except for underscore and any digit) 
#' and any batch substring are cleaned:
#'  * G00060245.18 => G00060245
#'  * G00060245.1-8 => G00060245
#'  * G02948263.1-1.DMA => G02948263
#'  * Gnumber followed by the codrug
#'    * G03252046.1-2;G00376771 => G03252046
#'  * Gnumber followed by the two codrugs
#'    * G03256376.1-2;G00376771.1-19;G02557755 => G03256376
#'  * Gnumber followed by the drug name
#'    * G00018838, Cisplatin => G00018838
#'
#' By default, Gnumber(s) followed by the "_" or digit (regardless the batch substring) are not cleaned:
#'  *  Gnumber with suffix added to prevent duplicated ids
#'     * G00060245_(G00060245.1-8)
#'  *  too long Gnumber 
#'     * G123456789.1-12
#' 
#' @param drug_vec character vector with drug id(s)
#' @param drug_p string with regex pattern for drug id. Set to Gnumber format by default: "G\[0-9\]\{8\}".
#' @param sep_p string with regex pattern for separator. Set to any character except for digit and space
#' @param batch_p string with regex pattern for batch substring. 
#'        By default set to any character(s): ".+"
#'
#' @examples
#' remove_drug_batch("G00060245.18")
#' remove_drug_batch("G00060245.1-8")
#' remove_drug_batch("G00060245.1-1.DMA")
#'
#' remove_drug_batch("G03252046.1-2;G00376771")
#' remove_drug_batch("G00018838, Cisplatin")
#' remove_drug_batch("G03256376.1-2;G00376771.1-19;G02557755")
#' remove_drug_batch("G00060245_(G00060245.1-8)")
#' remove_drug_batch(c("G00060245.18", "G00060245.1-8", "G00060245.1-1.DMA"))
#' 
#' remove_drug_batch("DRUG_01.123", drug_p = "DRUG_[0-9]+")
#' remove_drug_batch("G00001234:22-1", sep_p = ":")
#' remove_drug_batch("G00001234.28", batch_p = "[0-9]+")
#'
#' @return charvec with Gnumber(s)
#' @export
#' @keywords package_utils
remove_drug_batch <- function(drug_vec,
                              drug_p = "^G[0-9]{8}",
                              sep_p = "[^0-9|^_]",
                              batch_p = ".+") {
  checkmate::assert_character(drug_vec)
  checkmate::assert_string(drug_p)
  checkmate::assert_string(sep_p)
  checkmate::assert_string(batch_p)
  
  p <- paste0("(", drug_p, ")", sep_p, batch_p, "$")
  r <- "\\1"
  sub(p, r, drug_vec)
}


#' Cap infinity values (Inf, -Inf) in the assay data
#'
#' @param conc_assay_dt assay data in data.table format with Concentration data
#' @param assay_dt assay data in data.table format with infinity values to be capped
#' @param experiment_name string with the name of the experiment
#' @param col string with column name to be capped in assay_dt ("xc50" by default)
#' @param scaling_factor number with scaling factor for min and max concentration values
#'                       final formulas are min / scaling_factor and max * scaling_factor
#'        
#' @examples
#'  # single-agent data 
#'  sdata <- get_synthetic_data("finalMAE_small")
#'  smetrics_data <- convert_se_assay_to_dt(sdata[[get_supported_experiments("sa")]], "Metrics")
#'  saveraged_data <- convert_se_assay_to_dt(sdata[[get_supported_experiments("sa")]], "Averaged")
#'  smetrics_data_capped <- cap_assay_infinities(saveraged_data, smetrics_data, experiment_name = "single-agent")
#'
#'  # combination data 
#'  cdata <- get_synthetic_data("finalMAE_combo_matrix_small")
#'  scaveraged_data <- convert_se_assay_to_dt(cdata[[get_supported_experiments("combo")]], "Averaged")
#'  scmetrics_data <- convert_se_assay_to_dt(cdata[[get_supported_experiments("combo")]], "Metrics")
#'  scmetrics_data_capped <- cap_assay_infinities(scaveraged_data, scmetrics_data, experiment_name = "combination")
#' 
#' @return data.table without replicates
#' @keywords package_utils
#' @export
cap_assay_infinities <- function(conc_assay_dt,
                               assay_dt,
                               experiment_name,
                               col = "xc50",
                               scaling_factor = 1) {
  
  checkmate::assert_data_table(conc_assay_dt)
  checkmate::assert_data_table(assay_dt)
  checkmate::assert_string(experiment_name)
  checkmate::assert_choice(experiment_name, get_supported_experiments())
  checkmate::assert_string(col)
  checkmate::assert_choice(col, colnames(assay_dt))
  checkmate::assert_number(scaling_factor, lower = 1)
  
  conc_col <- if (experiment_name %in% c(get_experiment_groups("single-agent")[["single-agent"]],
                                         get_experiment_groups("combination"))) {
    
    # in combination experiments there is a matrix of drug1 X drug2 concentrations
    # as the matrix is symmetric the values of concentrations for drug1 ("Concentration")
    # and drug2 ("Concentration_2") are identical
    # thus the logic for single-agent and combination experiment is identical in this case
    get_env_identifiers("concentration")
  } else {
    stop(sprintf("unsupported experiment:'%s'", experiment_name))
  }
  
  # remove records for 0 concentrations 
  conc_assay_dt <- conc_assay_dt[conc_assay_dt[[conc_col]] != 0, ]
  
  group_cols <- if (experiment_name == get_experiment_groups("single-agent")[["single-agent"]]) {
    as.character(gDRutils::get_env_identifiers(c("drug", "cellline"), simplify = FALSE))
  } else if (experiment_name == get_experiment_groups("combination")) {
    as.character(gDRutils::get_env_identifiers(c("drug", "drug2", "cellline"), simplify = FALSE))
  } else {
    sprintf("unsupported experiment:'%s'", experiment_name)
  }
  
  out_dt <- if (any(assay_dt[[col]] %in% c(Inf, -Inf))) {
    min_max_conc <- conc_assay_dt[, .(min = min(get(conc_col)), max = max(get(conc_col))), by = group_cols]
    mt <- merge(assay_dt, min_max_conc, by = group_cols)
    mt[mt[[col]] == -Inf, col] <- mt[mt[[col]] == -Inf, "min"] / scaling_factor
    mt[mt[[col]] == Inf, col] <- mt[mt[[col]] == Inf, "max"] * scaling_factor
    data.table::setkey(mt, NULL)
    mt[, -c("min", "max")]
    
  } else {
    assay_dt
  }
  out_dt
  
}
