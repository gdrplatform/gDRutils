#' Convert a SummarizedExperiment assay to a long data.table
#'
#' Convert an assay within a \linkS4class{SummarizedExperiment} object to a long data.table.
#'
#' @details NOTE: to extract information about 'Control' data, simply call the
#' function with the name of the assay holding data on controls.
#' To extract the reference data in to same format as 'Averaged' use \code{convert_se_ref_assay_to_dt}.
#'
#' @param se A \linkS4class{SummarizedExperiment} object holding raw and/or processed dose-response data in its assays.
#' @param assay_name String of name of the assay to transform within the \code{se}.
#' @param include_metadata Boolean indicating whether or not to include \code{rowData(se)}
#' and \code{colData(se)} in the returned data.table.
#' Defaults to \code{TRUE}.
#' @param retain_nested_rownames Boolean indicating whether or not to retain the rownames 
#' nested within a \code{BumpyMatrix} assay.
#' Defaults to \code{FALSE}.
#' If the \code{assay_name} is not of the \code{BumpyMatrix} class, this argument's value is ignored.
#' If \code{TRUE}, the resulting column in the data.table will be named as \code{"<assay_name>_rownames"}.
#' @param wide_structure Boolean indicating whether or not to transform data.table into wide format.
#' `wide_structure = TRUE` requires `retain_nested_rownames = TRUE`.
#' @param unify_metadata Boolean indicating whether to unify DrugName and CellLineName in cases where DrugNames
#' and CellLineNames are shared by more than one Gnumber and/or clid within the experiment.
#' @param drop_masked Boolean indicating whether to drop masked values; TRUE by default.
#' @param merge_additional_variables Boolean indicating whether to merge additional variables identified by
#' \code{get_additional_variables} into the \code{DrugName} column. Defaults to \code{FALSE}.
#' @keywords convert
#'
#' @return data.table representation of the data in \code{assay_name}.
#'
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' convert_se_assay_to_dt(se, "Metrics")
#' 
#' @seealso flatten
#' @export
convert_se_assay_to_dt <- function(se,
                                   assay_name,
                                   include_metadata = TRUE,
                                   retain_nested_rownames = FALSE,
                                   wide_structure = FALSE,
                                   unify_metadata = FALSE,
                                   drop_masked = TRUE,
                                   merge_additional_variables = FALSE) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_flag(include_metadata)
  checkmate::assert_flag(retain_nested_rownames)
  checkmate::assert_flag(wide_structure)
  checkmate::assert_flag(unify_metadata)
  checkmate::assert_flag(merge_additional_variables)
  validate_se_assay_name(se, assay_name)
  if (wide_structure) {
    # wide_structure works only with `normalization_type` column in the assay 
    # and only for assays class "BumpyMatrix"
    if (!inherits(SummarizedExperiment::assay(se, assay_name), "BumpyDataFrameMatrix")) {
      warning("assay is not class `BumpyMatrix`, wide_structure=TRUE ignored")
      wide_structure <- FALSE
    } else if ("normalization_type" %in%
               BumpyMatrix::commonColnames(SummarizedExperiment::assay(se, assay_name))) {
      retain_nested_rownames <- TRUE
    } else {
      warning("'normalization_type' not found in assay, wide_structure=TRUE ignored")
      wide_structure <- FALSE
    }
  }
  dt <- .convert_se_assay_to_dt(se, assay_name, retain_nested_rownames = retain_nested_rownames)
  if (nrow(dt) == 0L) {
    return(dt)
  }
  if (drop_masked) {
    conc <- get_env_identifiers("concentration")
    masked_tag <- get_env_identifiers("masked_tag")
    if ("normalization_type" %in% names(dt)) {
      dt <- dt[!is.na(normalization_type)]
    }
    if (conc %in% names(dt)) {
      dt <- dt[!is.na(get(conc))]
    }
    if (masked_tag %in% names(dt)) {
      dt <- dt[get(masked_tag) == FALSE]
    }
  }
  if (include_metadata) {
    dt <- .extract_and_merge_metadata(se, data.table::copy(dt))
    
    if (merge_additional_variables) {
      additional_vars <- get_additional_variables(list(dt)) 

      if (!is.null(additional_vars) && length(additional_vars) > 0) {
        dt <- update_drug_name(dt, additional_vars)
      }
    }
  }
  if (wide_structure) {
    id_col <- paste0(assay_name, "_rownames")
    dt$id <- gsub("_.*", "", dt[[id_col]])
    dt[[id_col]] <- NULL
    normalization_cols <- unique(c(grep("^x$|x_+", names(dt), value = TRUE),
                                   intersect(unlist(get_header()[c("excess", "scores", "response_metrics")]),
                                             names(dt))))
    rest_cols <- setdiff(colnames(dt), c(normalization_cols, "normalization_type"))
    dcast_formula <- paste0(paste0(rest_cols, collapse = " + "), " ~  normalization_type")
    new_cols <- as.vector(outer(normalization_cols, unique(dt$normalization_type),
                                paste, sep = "_"))
    new_cols_rename <- unlist(lapply(strsplit(new_cols, "_"), function(x) {
      x[length(x)] <- extend_normalization_type_name(x[length(x)])
      if (grepl("^x$|x_+", x[1])) {
        paste(x[-1], collapse = "_")
      } else {
        paste(x, collapse = "_")
      }
    }))
    dt <- data.table::dcast(dt, dcast_formula, value.var = normalization_cols)
    dt$id <- NULL
    if (!all(new_cols %in% names(dt))) {
      new_cols <- gsub("x_", "", new_cols)
    }
    if (!any(duplicated(new_cols_rename))) {
      data.table::setnames(dt, new_cols, new_cols_rename, skip_absent = TRUE)
    }
  }
  if (unify_metadata) {
    dt <- gDRutils::set_unique_drug_names_dt(dt)
    dt <- gDRutils::set_unique_cl_names_dt(dt)
  }
  dt
}

#' @keywords internal
#' @return data.table containing merged assay data and metadata.
#' @noRd
.extract_and_merge_metadata <- function(se, dt) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_data_table(dt)
  
  rData <- data.table::as.data.table(rowData(se))
  rData[, rId := rownames(se)]
  cData <- data.table::as.data.table(colData(se))
  cData[, cId := colnames(se)]
  
  ids <- data.table::CJ(cData$cId, rData$rId)
  data.table::setnames(ids, c("cId", "rId"))
  ids[, names(ids) := lapply(.SD, as.character), .SDcols = names(ids)]
  annotations <- cData[rData[ids, on = "rId"], on = "cId"][dt, on = c("rId", "cId")]
  data.table::setcolorder(annotations, Reduce(union, list(names(dt), names(rData), names(cData))))
  annotations
}

#' Convert assay data into data.table.
#' @return data.table of assay data.
#' @keywords internal
#' @noRd
#'
.convert_se_assay_to_dt <- function(se, assay_name, retain_nested_rownames) {
  
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  
  object <- assays(se)[[assay_name]]
  checkmate::assert_true(inherits(object, "BumpyDataFrameMatrix") || inherits(object, "matrix"))
  
  rowfield <- "rId"
  colfield <- "cId"
  
  if (methods::is(object, "BumpyDataFrameMatrix")) {
    as_df <- BumpyMatrix::unsplitAsDataFrame(object, row.field = rowfield, column.field = colfield)
    # Retain nested rownames.
    if (retain_nested_rownames) {
      if (is.character(rownames(as_df))) {
        as_df[[paste0(assay_name, "_rownames")]] <- rownames(as_df)
      }
    }
    as_dt <- data.table::as.data.table(as_df)
    
  } else if (methods::is(object, "matrix")) {
    first <- object[1, 1][[1]]
    if (is.numeric(first)) {
      as_dt <-
        data.table::melt(data.table::as.data.table(object, keep.rownames = TRUE),
                         measure.vars = colnames(object))
      data.table::setnames(as_dt, c(rowfield, colfield, assay_name))
    } else {
      stop(sprintf("matrix with nested objects of class '%s' is not supported", class(first)))
    }
  } else {
    stop(sprintf("assay of class '%s' is not supported", class(object)))
  }
  as_dt
}


########################
# Convenience functions
########################

#' Convert a MultiAssayExperiment assay to a long data.table
#'
#' Convert an assay within a \linkS4class{SummarizedExperiment} object in a MultiAssayExperiment
#' to a long data.table.
#'
#' @details NOTE: to extract information about 'Control' data, simply call the
#' function with the name of the assay holding data on controls.
#'
#' @param mae A \linkS4class{MultiAssayExperiment} object holding experiments with 
#' raw and/or processed dose-response data in its assays.
#' @param assay_name String of name of the assay to transform within an experiment of the \code{mae}.
#' @param experiment_name String of name of the experiment in \code{mae} whose \code{assay_name} should be converted.
#' Defaults to \code{NULL} to indicate to convert assay in all experiments into one data.table object.
#' @param include_metadata Boolean indicating whether or not to include \code{rowData()}
#' and \code{colData()} in the returned data.table.
#' Defaults to \code{TRUE}.
#' @param retain_nested_rownames Boolean indicating whether or not to retain the rownames 
#' nested within a \code{BumpyMatrix} assay.
#' Defaults to \code{FALSE}.
#' If the \code{assay_name} is not of the \code{BumpyMatrix} class, this argument's value is ignored.
#' If \code{TRUE}, the resulting column in the data.table will be named as \code{"<assay_name>_rownames"}.
#' @param wide_structure Boolean indicating whether or not to transform data.table into wide format.
#' `wide_structure = TRUE` requires `retain_nested_rownames = TRUE` however that will be validated 
#' in `convert_se_assay_to_dt` function
#' @param drop_masked Boolean indicating whether to drop masked values; TRUE by default.
#' @param merge_additional_variables Boolean indicating whether to merge additional variables identified by
#' \code{get_additional_variables} into the \code{DrugName} column. Defaults to \code{FALSE}.
#' @keywords convert
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
#' 
#' @return data.table representation of the data in \code{assay_name}.
#'
#' @seealso flatten convert_se_assay_to_dt
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small")
#' convert_mae_assay_to_dt(mae, "Metrics")
#' 
#' @export
convert_mae_assay_to_dt <- function(mae,
                                    assay_name,
                                    experiment_name = NULL,
                                    include_metadata = TRUE,
                                    retain_nested_rownames = FALSE,
                                    wide_structure = FALSE,
                                    drop_masked = TRUE,
                                    merge_additional_variables = FALSE) {
  
  # Assertions.
  checkmate::assert_class(mae, "MultiAssayExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_choice(experiment_name, names(mae), null.ok = TRUE)
  checkmate::assert_flag(include_metadata)
  checkmate::assert_flag(retain_nested_rownames)
  checkmate::assert_flag(wide_structure)
  checkmate::assert_flag(merge_additional_variables)
  
  if (is.null(experiment_name)) {
    experiment_name <- names(mae)
  }
  
  dtList <- lapply(experiment_name, function(x) {
    if (!assay_name %in% assayNames(mae[[x]])) {
      return()
    }
    convert_se_assay_to_dt(mae[[x]],
                           assay_name = assay_name,
                           include_metadata = include_metadata,
                           retain_nested_rownames = retain_nested_rownames,
                           wide_structure = wide_structure,
                           drop_masked = drop_masked,
                           merge_additional_variables = merge_additional_variables)
  })
  if (all(vapply(dtList, is.null, logical(1)))) {
    warning(sprintf("assay '%s' was not found in any of the following experiments: '%s'",
                    assay_name,
                    paste(experiment_name, collapse = ", ")))
  }
  data.table::rbindlist(dtList, fill = TRUE, use.names = TRUE)
}


#' Convert a SummarizedExperiment assay to a long data.table and conduct some post processing steps
#'
#' Convert an assay within a SummarizedExperiment object to a long data.table. Then
#' conduct some post processing steps.
#'
#' Current strategy is per-assay specific.
#' 1. combo assays: conversion to data.table only (with `wide_structure` = FALSE)
#' 2. 'Metrics' assay can be converted to three types of outputs:
#'   - Metrics_initial (conversion to data.table only, with `wide_structure` = FALSE)
#'   - Metrics_raw: same as Metrics_initial followed by:
#'     * fix for 'EC50' and 'Metrics_rownames'
#'     * flatten
#'     * prettifying and dropping excess variables
#'   - Metrics (same as Metrics_raw + cap_values if `cap_values = TRUE`)
#' 3. 'Normalization' and 'Averaged' assay:
#'   - conversion to data.table (with `wide_structure` = TRUE)
#'   - prettifying and dropping excess variables
#'
#' @details NOTE: to extract information about 'Control' data, simply call the
#' function with the name of the assay holding data on controls.
#' To extract the reference data in the same format as 'Averaged' use \code{convert_se_ref_assay_to_dt}.
#'
#' @param se A SummarizedExperiment object holding raw and/or processed dose-response data in its assays.
#' @param assay_name String of name of the assay to transform within the \code{se}.
#' @param output_table String of type name of the output data.table.
#' @param cap_values Logical indicating whether to apply capping (via `capVals`) for "Metrics" output. Default is FALSE.
#' @return data.table representation of the data in \code{assay_name} with added information from \code{colData}.
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' convert_se_assay_to_custom_dt(se, "Metrics")
#' convert_se_assay_to_custom_dt(se, "Metrics", output_table = "Metrics_raw")
#' convert_se_assay_to_custom_dt(se, "Metrics", output_table = "Metrics_initial")
#' convert_se_assay_to_custom_dt(se, "Averaged")
#' convert_se_assay_to_custom_dt(se, "Metrics", cap_values = TRUE)
#'
#' @seealso convert_se_assay_to_dt
#' @keywords convert
#' @export
convert_se_assay_to_custom_dt <- function(se,
                                          assay_name,
                                          output_table = NULL,
                                          cap_values = FALSE) {
  
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_string(output_table, null.ok = TRUE)
  checkmate::assert_choice(tolower(assay_name), as.character(tolower(get_assay_names())))
  checkmate::assert_choice(output_table,
                           c(get_assay_names(), "Metrics_initial", "Metrics_raw"),
                           null.ok = TRUE)
  checkmate::assert_flag(cap_values)
  
  if (is.null(output_table)) {
    output_table <- assay_name
  }
  if (output_table %in% c("Metrics_initial", "Metrics_raw")) {
    stopifnot(assay_name == "Metrics")
  }
  
  wide_structure <- assay_name %in% c("Normalized", "Averaged")
  dt <- convert_se_assay_to_dt(se,
                               assay_name,
                               include_metadata = TRUE,
                               wide_structure = wide_structure)
  
  if (output_table %in% c(get_combo_assay_names(), "Metrics_initial")) {
    return(dt)
  }
  
  if (output_table %in% c("Metrics", "Metrics_raw")) {
    # SE*.qs files contain 'c50' column instead of the 'ec50' in the metrics
    # this is a temporary fix that should be removed once data (qs files) is reprocessed
    data.table::setnames(dt, "c50", "ec50", skip_absent = TRUE)
    
    groups <- c("normalization_type", "fit_source")
    if (all(groups %in% names(dt))) {
      dt <- flatten(
        tbl = dt,
        groups = groups,
        wide_cols = get_header("response_metrics")
      )
    }
  }
  
  if (output_table %in% c("Metrics_raw", "Metrics", "Normalized", "Averaged")) {
    # TODO GDR-2513 # nolint start
    # pidfs <- get_SE_identifiers(se)
    # udrugs <- unique(as.character(dt[[pidfs[["drug_name"]]]]))
    # if (all(udrugs %in% get_SE_identifiers(se, "untreated_tag"))) {
    #   futile.logger::flog.trace("Main App: \t dropping excess variables", name = "trace.logger")
    #   vars <-
    #     as.character(pidfs[c("drug_name2", "concentration2", "drug2", "drug_moa2")])
    #   dt[vars] <- NULL
    # }
    # nolint end
    
    # add identifiers specific for given SE
    colnames(dt) <- prettify_flat_metrics(colnames(dt), human_readable = TRUE)
  }
  
  if (output_table == "Metrics" && cap_values) {
    dt <- capVals(dt)
  }
  
  dt
}

#' Cap metric values
#'
#' Convenience function to apply caps to outlying metric values.
#'
#' The following metrics are capped at the respective values:
#' \itemize{
#'   \item{\code{E max}: 0  - 1.1}
#'   \item{\code{GR max}: -1  - 1.1}
#'   \item{\code{RV AOC within set range}: over -0.1}
#'   \item{\code{GR AOC within set range}: over of -0.1}
#'   \item{\code{GR50}: 1e-4 to 30}
#'   \item{\code{IC50}: 1e-4 to 30}
#'   \item{\code{EC50}: 1e-4 to 30} (change 0 to NA beforehand)
#' }
#'
#' @param x \code{data.table} containing growth metrics extracted from a \code{SummarizedExperiment}
#' 
#' @examples
#' dt <- data.table::data.table(
#'   `E Max` = c(-0.1, 0, 0.5, 1.2),
#'   `GR Max` = c(-1.1, -1, 0.5, 1.2),
#'   `RV AOC within set range` = c(-0.2, -0.1, 0, 3),
#'   `GR AOC within set range` = c(-0.2, -0.1, 0, 3), 
#'   `GR50` = c(0, 1e-7, 10, 34),
#'   `IC50` = c(0, 1e-7, 10, 34),
#'   `EC50` = c(0, 1e-7, 10, 34),
#'   check.names = FALSE
#' )
#' dt
#' dt1 <- capVals(dt)
#' dt1
#' 
#' @return A data table with capped values.
#' @keywords internal
#'
#' @seealso \code{convert_se_assay_to_dt}, \code{\link[scales]{oob}}
#'
#' @export
capVals <- function(x) {
  
  checkmate::assert_data_table(x)
  
  json_path <- system.file(package = "gDRutils", "settings.json")
  s <- get_settings_from_json("capVals", json_path)
  # fifty_lower_limit numeric value of the lower limit to cap all x50 metrics
  checkmate::assert_number(s$fifty_lower_limit)
  # fifty_upper_limit numeric value of the upper limit to cap all x50 metrics
  checkmate::assert_number(s$fifty_upper_limit)
  # max_upper_limit numeric value of the upper limit to cap all xmax metrics
  checkmate::assert_number(s$max_upper_limit)
  # range_lower_limit numeric value of the lower limit to cap all xrange metrics
  checkmate::assert_number(s$range_lower_limit)
  
  s_col <- get_settings_from_json("CAP_VALS_COLS", json_path)
  if (!NROW(intersect(s_col, names(x)))) return(x) # no columns to capped
  
  X <- data.table::copy(x)
  if ("E Max" %in% names(X)) {
    X[, `E Max` := scales::oob_squish_any(`E Max`, range = c(0, s$max_upper_limit))]
  }
  if ("GR Max" %in% names(X)) {
    X[, `GR Max` := scales::oob_squish_any(`GR Max`, range = c(-1, s$max_upper_limit))]
  }
  if ("RV AOC within set range" %in% names(X)) {
    X[, `RV AOC within set range` := scales::oob_squish_any(`RV AOC within set range`,
                                                            range = c(s$range_lower_limit, NA))]
  }
  if ("GR AOC within set range" %in% names(X)) {
    X[, `GR AOC within set range` := scales::oob_squish_any(`GR AOC within set range`,
                                                            range = c(s$range_lower_limit, NA))]
  }
  if ("GR50" %in% names(X)) {
    X[, GR50 := scales::oob_squish_any(GR50, range = c(s$fifty_lower_limit, s$fifty_upper_limit))]
  }
  if ("IC50" %in% names(X)) {
    X[, IC50 := scales::oob_squish_any(IC50, range = c(s$fifty_lower_limit, s$fifty_upper_limit))]
  }
  if ("EC50" %in% names(X)) {
    ## An EC50 value of 0 indicates a flat fit.
    ## As such, it is appropriate to set the value to NA,
    ## and blank entries in the metric clustering is to be expected.
    X[EC50 == 0, EC50 := NA]
    X[, EC50 := scales::oob_squish_any(EC50, range = c(s$fifty_lower_limit, s$fifty_upper_limit))]
  }
  return(X)
}

#' Update drug name with additional variables
#'
#' Concatenates the values of specified additional variables to the existing
#' drug identifier columns in a data.table, using the variables defined in
#' \code{get_env_identifiers}.
#'
#' @param dt A data.table containing drug-response information, including drug
#' identifier columns (e.g., \code{DrugName}, \code{Gnumber}) and the \code{additional_vars}.
#' @param additional_vars Character vector of column names (variables) to merge
#' into the drug identifier columns.
#'
#' @return A copy of the input data.table \code{dt} with the relevant drug
#' identifier columns updated to include the additional variable information in the format:
#' \code{Identifier (variable = value)}.
#'
#' @examples
#' # Assuming get_env_identifiers() returns c("DrugName", "Gnumber") for drug identifiers
#' dt <- data.table::data.table(
#'   DrugName = c("DrugA", "DrugA", "DrugB"),
#'   Gnumber = c("G1", "G1", "G2"),
#'   Var1 = c(NA, "X", NA),
#'   Var2 = c(NA, "Y", "Z")
#' )
#' additional_vars <- c("Var1", "Var2")
#' # update_drug_name(dt, additional_vars) # Would update DrugName and Gnumber
#'
#' @keywords internal
#' @export
update_drug_name <- function(dt, additional_vars) {
  checkmate::assert_data_table(dt)
  checkmate::assert_character(additional_vars)
  
  dt <- data.table::copy(dt)
  
  # Identify the columns to merge the additional info into
  cols_to_merge <- unlist(get_env_identifiers(c("drug", "drug_name"), simplify = FALSE))
  
  for (var in additional_vars) {
    if (!var %in% names(dt)) {
      warning(sprintf("Additional variable '%s' not found in data.table. Skipping merge for this variable.", var))
      next
    }
    
    # Iterate over all drug identifier columns
    for (col in cols_to_merge) {
      if (!col %in% names(dt)) {
        warning(sprintf("Drug identifier column '%s' not found in data.table. Skipping update for this column.", col))
        next
      }

      dt[, (col) := ifelse(
        is.na(dt[[var]]),
        get(col),
        paste0(get(col), " (", var, " = ", get(var), ")")
      )]
    }
  }
  return(dt)
}