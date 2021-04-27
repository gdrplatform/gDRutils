#' Convert a SummarizedExperiment assay to a long data.table
#'
#' Convert an assay within a \linkS4class{SummarizedExperiment} object to a long data.table.
#'
#' @param se A \linkS4class{SummarizedExperiment} object holding raw and/or processed dose-response data in its assays.
#' @param assay_name String of name of the assay to transform within the \code{se}.
#' @param include_metadata Boolean indicating whether or not to include \code{rowData(se)}
#' and \code{colData(se)} in the returned data.table.
#' Defaults to \code{TRUE}.
#'
#' @return data.table representation of the data in \code{assay_name}.
#'
#' @details NOTE: to extract information about 'Control' data, simply call the
#' function with the name of the assay holding data on controls.
#' To extract the reference data in to same format as 'Averaged' use \code{convert_se_ref_assay_to_dt}.
#'
#' @export
#'
convert_se_assay_to_dt <- function(se,
                                   assay_name,
                                   include_metadata = TRUE) {

  # Assertions.
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::test_string(assay_name)
  checkmate::assert_flag(include_metadata)
  checkmate::assert_character(sel_metric_type)	
  
  validate_se_assay_name(se, assay_name)

  dt <- .convert_se_assay_to_dt(se, assay_name)

  if (nrow(dt) == 0L) {
    return(dt) # TODO: Should this return something else?
  }

  if (include_metadata) {
    rData <- SummarizedExperiment::rowData(se)
    rData$rId <- rownames(rData)

    cData <- SummarizedExperiment::colData(se)
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
.convert_se_assay_to_dt <- function(se, assay_name) {
  object <- SummarizedExperiment::assays(se)[[assay_name]]
  checkmate::assert_true(inherits(object, "BumpyDataFrameMatrix") || inherits(object, "matrix"))

  if (methods::is(object, "BumpyDataFrameMatrix")) {
    as_df <- BumpyMatrix::unsplitAsDataFrame(object, row.field = "rId", column.field = "cId")
    # Retain nested rownames.
    as_df[[paste0(assay_name, "_rownames")]] <- rownames(as_df)

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

#' Convert the reference values from a SummarizedExperiment assay to a long data.table
#'
#' Transform the Ref[RelativeViability/GRvalue] within a \linkS4class{SummarizedExperiment} object to a long data.table.
#' Clean up the column names and add columns to match the format of the data.table from the \code{'Averaged'} assay.
#'
#' @param se A \linkS4class{SummarizedExperiment} object holding reference data in its assays.
#' @param ref_relative_viability_assay String of the name of the assay in the \code{se} 
#' holding the reference relative viability data.
#' @param ref_gr_value_assay String of the name of the assay in the \code{se} holding the reference GR value data.
#'
#' @return data.table representation of the reference data.
#'
#' @details This is a convenience function to massage the reference data into the same format as the \code{"Averaged"}
#' assay data.
#'
#' @export
#'
convert_se_ref_assay_to_dt <- function(se,
                                       ref_relative_viability_assay = "RefRelativeViability",
                                       ref_gr_value_assay = "RefGRvalue") {
  rv <- convert_se_assay_to_dt(se, ref_relative_viability_assay, include_metadata = TRUE)
  colnames(rv)[colnames(rv) == ref_relative_viability_assay] <- "RelativeViability"
  rv$std_RelativeViability <- NA

  gr <- convert_se_assay_to_dt(se, ref_gr_value_assay)
  colnames(gr)[colnames(gr) == ref_gr_value_assay] <- "GRvalue"
  gr$std_GRvalue <- NA

  dt <- merge(rv, gr, all = TRUE)

  # Fill primary drug with 'untreated_tag'.
  dt$Concentration <- 0
  untreated <- get_SE_identifiers(se, "untreated_tag")[1]
  dt[, get_SE_identifiers(se, "drug")] <- untreated
  dt[, get_SE_identifiers(se, "drugname")] <- untreated
  dt[, get_SE_identifiers(se, "drug_moa")] <- untreated

  data.table::as.data.table(dt)
}


#' Flatten a data.table
#'
#' Flatten a stacked data.table.
#'
#' @param columns character vector of column names representing uniquifying groups in expansion.
#' @param flatten character vector of column names to flatten.
#' @param sep 
#' Defaults to \code{"_"}.
#' @return data.table with assay data
#' @details flattened columns will be named according to the columns.
#' This is useful in trying to get a flattened version of the \code{"Metrics"} assay.
#' @export
#'
flatten_stacked <- function(dt, columns, flatten, sep = "_") {
  if (!all(columns %in% colnames(dt))) {	
    stop(sprintf("missing expected uniquifying columns: '%s'",
      paste0(setdiff(colnames(dt), columns), collapse = ", ")))
  }	
  
  idx <- which(columns %in% colnames(dt))
  uniquifying <- dt[, idx, drop = FALSE]
  uniquifying <- unique(uniquifying)

  out <- split(dt[, -idx], dt[, idx, drop = FALSE], sep = sep)
  missing <- setdiff(flatten, colnames(dt))
  if (length(missing) != 0L) {
    warning(sprintf("missing listed flatten columns: '%s'", paste0(missing, collapse = ", ")))
  }

  rename <- colnames(out[[1]]) %in% flatten 
  for (grp in names(out)) {
    group <- out[[grp]]
    colnames(group)[rename] <- paste0(grp, sep, colnames(group)[rename])
    out[[grp]] <- group
  }

  Reduce(merge, out)
}
