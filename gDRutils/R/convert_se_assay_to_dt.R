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
#' To extract the reference data in to same format as 'Averaged' use \code{convert_se_ref_assay_to_dt}
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

  is_valid_se_assay_name(se, assay_name)

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


#' Convert the reference values from a SummarizedExperiment assay to a long data.table
#'
#' Transform the Ref[RelativeViability/GRvalue] within a \linkS4class{SummarizedExperiment} object to a long data.table.
#' Clean up the column names and add columns to match the format of the data.table from the \code{'Averaged'} assay.
#'
#' @param se A \linkS4class{SummarizedExperiment} object holding reference data in its assays.
#' @param ref_relative_viability_assay String of the name of the assay in the \code{se} holding the reference relative viability data.
#' @param ref_gr_value_assay String of the name of the assay in the \code{se} holding the reference GR value data.
#'
#' @return data.table representation of the reference data.
#'
#' @details This is a convenience function to massage the reference data into the same format as the code{"Averaged"} assay data.
#'
#' @export
#'
convert_se_ref_assay_to_dt <- function(se,
                                       ref_relative_viability_assay = "RefRelativeViability",
                                       ref_gr_value_assay = "RefGRvalue") {
  rv <- convert_se_assay_to_dt(se, ref_relative_viability_assay, include_metadata=TRUE)
  colnames(rv)[colnames(rv) == ref_relative_viability_assay] <- "RelativeViability"
  rv$std_RelativeViability <- NA

  gr <- convert_se_assay_to_dt(se, ref_gr_value_assay)
  colnames(gr)[colnames(gr) == ref_gr_value_assay] <- "GRvalue"
  gr$std_GRvalue <- NA

  dt <- merge(rv, gr, all = TRUE)

  # Fill primary drug with 'untreated_tag'.
  dt$Concentration <- 0
  # TODO: See if these identifiers are already on the SE object.
  dt[, gDRutils::get_identifier("drug")] <- get_identifier("untreated_tag")[1]
  dt[, gDRutils::get_identifier("drugname")] <- get_identifier("untreated_tag")[1]
  dt[, gDRutils::get_identifier("drug_moa")] <- get_identifier("untreated_tag")[1]

  data.table::as.data.table(dt)
}


#' Convert assay data into data.table
#' @return data.table with assay data
#' @keywords internal
#' @noRd
#'
.convert_se_assay_to_dt <- function(se, assay_name) {
  object <- SummarizedExperiment::assays(se)[[assay_name]]
  checkmate::assert_true(inherits(object, "BumpyDataFrameMatrix") || inherits(object, "matrix"))

  if (is(object, "BumpyDataFrameMatrix")) {
    as_df <- BumpyMatrix::unsplitAsDataFrame(object, row.field = "rId", column.field = "cId")

  } else if (methods::is(object, "matrix")) {
    first <- object[1, 1][[1]]
    if (is.numeric(first)) {
      as_df <- reshape2::melt(object, varnames = c("rId", "cId"), value.name = assay_name)
    } else if (methods::is(first, "DFrame") || methods::is(first, "data.frame")) {

      # TODO: Deprecate me.
      .Deprecated(msg = paste("support for nested DataFrames of class `matrix`",
                              "will be deprecated in the next release cycle.",
                              "See `BumpyDataFrameMatrix` instead"))
      asL <-
        lapply(seq_len(ncol(object)), function(x) {
          myL <- object[, x, drop = FALSE]
          names(myL) <- rownames(object)

      # in some datasets there might be no data for given drug/cell_line combination
      # under such circumstances DataFrame will be empty
      # and should be filtered out
      # example: testdata 7 - SummarizedExperiment::assay(seL[[7]],"Normalized")[5:8,1]
          myL <- myL[vapply(myL, nrow, integer(1)) > 0]

          myV <- vapply(myL, nrow, integer(1))
          rCol <- rep(names(myV), as.numeric(myV))
      # there might be DataFrames with different number of columns
      # let's fill with NAs where necessary
      if (length(unique(sapply(myL, ncol))) > 1) {
        df <-
          do.call(plyr::rbind.fill,
            lapply(myL, data.table::as.data.table))
      } else {
        df <-
          data.table::rbindlist(lapply(myL, data.table::as.data.table), fill = TRUE)
      }
      if (nrow(df) == 0)
        return()
      df$rId <- rCol
      df$cId <- colnames(object)[x]
      df
    })
      as_df <- data.table::rbindlist(asL)
    }
  as_df
  }
}


# TODO: Deprecate me.
#' Transform a SummarizedExperiment assay to a long data.table
#'
#' Transform a SummarizedExperiment assay to a long data.table with a single entry for each row and column combination.
#'
#' @param se \linkS4class{SummarizedExperiment} object with dose-response data.
#' @param assay_name String of name of the assay or index of the assay in the \code{se}.
#' @param merge_metrics Logical indicating whether the metrics should be merged.
#' Defaults to \code{FALSE}.
#' @param include_metadata Boolean indicating whether to include the metadata on the SummarizedExperiment.
#' Defaults to \code{TRUE}.
#'
#' @return data.table with dose-response data
#'
#' @export
#'
assay_to_dt <- function(se,
                        assay_name,
                        merge_metrics = FALSE,
                        include_metadata = TRUE) {

  # Assertions.
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assertTRUE(checkmate::test_count(assay_name) ||
                          checkmate::test_string(assay_name))
  checkmate::assert_flag(merge_metrics)

  .Deprecated(new = convert_se_assay_to_dt,
    msg = "support for 'assay_to_dt' will be dropped next release cycle. See 'convert_se_assay_to_dt' instead")

  if (is.integer(assay_name)) {
    assay_name <- SummarizedExperiment::assayNames(se)[assay_name]
  }

  as_dt <- convert_se_assay_to_dt(se, assay_name, include_metadata = include_metadata)
  if (assay_name == "Metrics") {
    ## TODO: Put in issue to BumpyMatrix::unsplitAsBumpyMatrix to also return nested rownames.
    ## Then can remove all hard-coded logic below regarding metrics.
    as_dt$dr_metric <- rep_len(c("RV", "GR"), nrow(as_dt))

    ## NOTE: assay_to_dt function is deprecated for convert_se_assay_to_dt function,
    ## which noteably does NOT support the merge_metrics argument.
    if (merge_metrics) {
      metric_headers <- get_header("response_metrics")
      if (!all(metric_headers %in% colnames(as_dt))) {
        stop(sprintf("missing expected metric headers: '%s'",
          paste0(setdiff(metric_headers, colnames(as_dt)), collapse = ", ")))
      }

      id_headers <- c("rId", "cId")
      headers <- c(id_headers, metric_headers)

      metric_headers <- colnames(as_dt)[colnames(as_dt) %in% metric_headers]

      Df_RV <- as_dt[dr_metric == "RV", ..headers]
      rv_map <- get_header("RV_metrics")

      Df_GR <- as_dt[dr_metric == "GR", - "dr_metric"]
      gr_map <- get_header("GR_metrics")

      data.table::setnames(Df_RV,
               old = metric_headers,
               new = unname(rv_map[metric_headers]))

      data.table::setnames(Df_GR,
               old = metric_headers,
               new = unname(gr_map[metric_headers]))

      as_dt <- merge(Df_RV,
                     Df_GR,
                     by = id_headers,
                     all = TRUE)
    }
  }
  return(as_dt)
}
