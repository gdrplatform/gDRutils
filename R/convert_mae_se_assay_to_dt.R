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
                                   wide_structure = FALSE) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_flag(include_metadata)
  checkmate::assert_flag(retain_nested_rownames)
  validate_se_assay_name(se, assay_name)
  if (wide_structure) {
    # wide_structure works only with `normalization_type` column in the assay
    if ("normalization_type" %in%
        BumpyMatrix::commonColnames(SummarizedExperiment::assay(se, assay_name))) {
      retain_nested_rownames <- TRUE 
    } else {
      wide_structure <- FALSE
    }
  }
  dt <- .convert_se_assay_to_dt(se, assay_name, retain_nested_rownames = retain_nested_rownames)
  if (nrow(dt) == 0L) {
    return(dt) # TODO: Should this return something else?
  }
  if (include_metadata) {
    dt <- .extract__and_merge_metadata(se, data.table::copy(dt))
  }
  if (wide_structure) {
    normalization_cols <- grep("^x$|x_+", names(dt), value = TRUE)
    id_col <- paste0(assay_name, "_rownames")
    dt$id <- gsub("_.*", "", dt[[id_col]])
    dt[[id_col]] <- NULL
    rest_cols <- setdiff(colnames(dt), c(normalization_cols, "normalization_type"))
    dcast_formula <- paste0(paste0(rest_cols, collapse = " + "), " ~  normalization_type")
    new_cols <- as.vector(outer(normalization_cols, unique(dt$normalization_type),
                                paste, sep = "_"))
    new_cols_rename <- unlist(lapply(strsplit(new_cols, "_"), function(x) {
      x[length(x)] <- gDRutils::extend_normalization_type_name(x[length(x)])
      paste(x[-1], collapse = "_")
      }))
    dt <- data.table::dcast(dt, dcast_formula, value.var = normalization_cols)
    dt$id <- NULL 
    if (!all(new_cols %in% names(dt))) {
      new_cols <- gsub("x_", "", new_cols)
    }
    data.table::setnames(dt, new_cols, new_cols_rename, skip_absent = TRUE)
  }
  dt
}

#' @keywords internal
.extract__and_merge_metadata <- function(se, dt) {
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

  if (methods::is(object, "BumpyDataFrameMatrix")) {
    as_df <- BumpyMatrix::unsplitAsDataFrame(object, row.field = "rId", column.field = "cId")
    # Retain nested rownames.
    if (retain_nested_rownames) {
      checkmate::assert_string(rownames(as_df))
      as_df[[paste0(assay_name, "_rownames")]] <- rownames(as_df)
    }
    as_dt <- data.table::as.data.table(as_df)

  } else if (methods::is(object, "matrix")) {
    first <- object[1, 1][[1]]
    if (is.numeric(first)) {
      as_dt <-
        data.table::melt(data.table::as.data.table(object, keep.rownames = TRUE),
                         measure.vars = colnames(object))
      data.table::setnames(as_dt, c("rId", "cId", assay_name))
    } else {
      stop(sprintf("matrix with nested objects of class '%s' is not supported", class(first)))
    }
    as_dt
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
#' To extract the reference data in to same format as 'Averaged' use \code{convert_mae_ref_assay_to_dt}.
#'
#' @param mae A \linkS4class{MultiAssayExperiment} object holding experiments with 
#' raw and/or processed dose-response data in its assays.
#' @param assay_name String of name of the assay to transform within the \code{se}.
#' @param experiment_name String of name of the experiment in `mae` whose `assay_name` should be converted.
#' Default to `NULL` that all the experiment should be converted into one data.table object.
#' @param include_metadata Boolean indicating whether or not to include \code{rowData(se)}
#' and \code{colData(se)} in the returned data.table.
#' Defaults to \code{TRUE}.
#' @param retain_nested_rownames Boolean indicating whether or not to retain the rownames 
#' nested within a \code{BumpyMatrix} assay.
#' Defaults to \code{FALSE}.
#' If the \code{assay_name} is not of the \code{BumpyMatrix} class, this argument's value is ignored.
#' If \code{TRUE}, the resulting column in the data.table will be named as \code{"<assay_name>_rownames"}.
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
#' 
#' @return data.table representation of the data in \code{assay_name}.
#'
#' @seealso flatten
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
                                    retain_nested_rownames = FALSE) {
  
  # Assertions.
  checkmate::assert_class(mae, "MultiAssayExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_choice(experiment_name, names(mae), null.ok = TRUE)
  checkmate::assert_flag(include_metadata)
  checkmate::assert_flag(retain_nested_rownames)
  
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
                           retain_nested_rownames = retain_nested_rownames)
  })
  data.table::rbindlist(dtList, fill = TRUE, use.names = TRUE)
}
