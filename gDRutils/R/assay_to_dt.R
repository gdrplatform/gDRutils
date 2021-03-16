#' Transform a SummarizedExperiment assay to a long data.table
#'
#' Transform a SummarizedExperiment assay to a concatenated data.table.
#'
#' @param se A \linkS4class{SummarizedExperiment} object with dose-response data in its assays.
#' @param assay_name String of name of the assay to transform within the \code{se}.
#' @param include_metadata Boolean indicating whether or not to include \code{rowData(se)}
#' and \code{colData(se)} in the response.
#' Defaults to \code{TRUE}.
#'
#' @return data.table of all data specified by \code{assay_name}.
#'
#' @details This is a base function that can be repurposed for more complex logic. 
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

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("'%s' is not on of the available assays: '%s'", 
      assay_name, paste0(SummarizedExperiment::assayNames(se))))
  }
 
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
    ## Then can remove this hard coding and assumptions. 
    as_dt$dr_metric <-  rep_len(c("RV", "GR"), nrow(as_dt))
    if (merge_metrics) {
      colnames_RV <- get_header("RV_metrics")
      diff_RV_columns <- setdiff(names(colnames_RV), colnames(as_dt))
      if (length(diff_RV_columns) > 0) {
        warning(paste(
          "missing column(s) in SE:",
          paste(diff_RV_columns, collapse = ", ")
        ))
        colnames_RV <- colnames_RV[!names(colnames_RV) %in% diff_RV_columns]
      }
      colnames_GR <- get_header("GR_metrics")
      diff_GR_columns <- setdiff(names(colnames_GR), colnames(as_dt))
      if (length(diff_GR_columns) > 0) {
        warning(paste(
          "missing column(s) in SE:",
          paste(diff_GR_columns, collapse = ", ")
        ))
        colnames_GR <- colnames_GR[!names(colnames_GR) %in% diff_GR_columns]
      }
      
      vars <- c("rId", "cId", names(colnames_RV))
      Df_RV <- as_dt[dr_metric == "RV", ..vars]
      Df_GR <- as_dt[dr_metric == "GR", - "dr_metric"]
      
      data.table::setnames(Df_RV,
                           old = names(colnames_RV),
                           new = unname(colnames_RV))
      data.table::setnames(Df_GR,
                           old = names(colnames_GR),
                           new = unname(colnames_GR))
      as_dt <- merge(Df_RV,
		     Df_GR,
		     by = c("rId", "cId"),
		     all = TRUE)
    }
  as_dt
  }
}


#' Convert assay data into data.table
#'
#' @param object An object comprising assay in SummarizedExperiment
#'
#' @return data.table with assay data
#'
.convert_se_assay_to_dt <- function(se, assay_name) {
  object <- SummarizedExperiment::assays(se)[[assay_name]]
  if (is(object, "BumpyDataFrameMatrix")) {
    as_df <- BumpyMatrix::unsplitAsDataFrame(object, row.field = "rId", column.field = "cId")
  } else if (is(object, "matrix")) {
    first <- object[1, 1]
    if (is.numeric(first)) {
      as_df <- reshape2::melt(object, varnames = c("rId", "cId"), value.name = assay_name)
    } else if (is(first, "DFrame")) {
      # TODO: Deprecate me. 
      .Deprecated("Support for nested DataFrames not of class `BumpyDataFrameMatrix` will be deprecated")
      asL <-
	lapply(seq_len(ncol(object)), function(x) {
	  myL <- object[, x, drop = FALSE]
	  
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

    } else {
      stop(
	paste(
	  "Unable to find an inherited method for function 'convert_assay_data_to_dt' for signature",
	  shQuote(class(object))
	)
      )
    }
  }
}
