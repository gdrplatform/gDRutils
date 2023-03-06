#' Convert a SummarizedExperiment assay to a long data.table
#'
#' Convert an assay within a \linkS4class{SummarizedExperiment} object to a long data.table.
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
#'
#' @return data.table representation of the data in \code{assay_name}.
#'
#' @details NOTE: to extract information about 'Control' data, simply call the
#' function with the name of the assay holding data on controls.
#' To extract the reference data in to same format as 'Averaged' use \code{convert_se_ref_assay_to_dt}.
#'
#' @seealso flatten
#' @export
#'
convert_se_assay_to_dt <- function(se,
                                   assay_name,
                                   include_metadata = TRUE,
                                   retain_nested_rownames = FALSE) {

  # Assertions.
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_flag(include_metadata)
  checkmate::assert_flag(retain_nested_rownames)
  
  validate_se_assay_name(se, assay_name)

  dt <- .convert_se_assay_to_dt(se, assay_name, retain_nested_rownames = retain_nested_rownames)

  if (nrow(dt) == 0L) {
    return(dt) # TODO: Should this return something else?
  }

  if (include_metadata) {
    rData <- rowData(se)
    rData$rId <- rownames(rData)

    cData <- colData(se)
    cData$cId <- rownames(cData)

    ids <- expand.grid(rData$rId, cData$cId)
    colnames(ids) <- c("rId", "cId")
    ids[] <- lapply(ids, as.character)

    annotations <- merge(ids, rData, by = "rId", all.x = TRUE)
    annotations <- merge(annotations, cData, by = "cId", all.x = TRUE)

    dt <- merge(dt, annotations, by = c("rId", "cId"), all.x = TRUE)
  }

  data.table::as.data.table(dt)
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
      as_df[[paste0(assay_name, "_rownames")]] <- rownames(as_df)
    }

  } else if (methods::is(object, "matrix")) {
    first <- object[1, 1][[1]]
    if (is.numeric(first)) {
      as_df <- reshape2::melt(object, varnames = c("rId", "cId"), value.name = assay_name)
    } else {
      stop(sprintf("matrix with nested objects of class '%s' is not supported", class(first)))
    }
    as_df
  } else {
    stop(sprintf("assay of class '%s' is not supported", class(object)))
  }
  data.table::as.data.table(as_df)
}


########################
# Convenience functions
########################

#' Convert a MultiAssayExperiment assay to a long data.table
#'
#' Convert an assay within a \linkS4class{SummarizedExperiment} object in a MultiAssayExperiment
#' to a long data.table.
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
#' @return data.table representation of the data in \code{assay_name}.
#'
#' @details NOTE: to extract information about 'Control' data, simply call the
#' function with the name of the assay holding data on controls.
#' To extract the reference data in to same format as 'Averaged' use \code{convert_mae_ref_assay_to_dt}.
#'
#' @seealso flatten
#' @export
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
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
  dt <- plyr::rbind.fill(dtList)
  dt
}
