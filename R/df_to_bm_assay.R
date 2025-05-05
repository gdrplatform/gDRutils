

#' df_to_bm_assay
#'
#' Convert data.table with dose-response data into a BumpyMatrix assay.
#'
#' The 'assay' is simply a BumpyMatrix object with rownames being the treatment ids,
#' colnames being the ids of the cell lines and values with
#' dose-response data for given cell lines under given conditions.
#'
#' @param data data.table with drug-response data
#' @param discard_keys a vector of keys that should be discarded
#' @keywords convert
#'
#' @examples
#' df_to_bm_assay(data.table::data.table(Gnumber = 2, clid = "A"))
#'
#' @return BumpyMatrix object
#'
#' @export
df_to_bm_assay <-
  function(data,
           discard_keys = NULL) {
    stopifnot(any(inherits(data, "data.table"), checkmate::test_character(data)))
    checkmate::assert_character(discard_keys, null.ok = TRUE)
    
    allMetadata <- split_SE_components(data, nested_keys = discard_keys)
    
    seColData <- data.table::as.data.table(allMetadata$condition_md)
    seRowData <- data.table::as.data.table(allMetadata$treatment_md)
    
    # Create row_id and col_id directly in data.table for speed
    seColData[, col_id := .I]
    seRowData[, row_id := .I]
    
    cl_entries <- setdiff(colnames(seColData), c("col_id", "name_"))
    cond_entries <- setdiff(colnames(seRowData), c("row_id", "name_"))
    
    # Use data.table's cross-join for speed
    complete <- data.table::CJ(col_id = seColData$col_id, row_id = seRowData$row_id)
    complete[, factor_id := .I]
    
    # Use data.table merge for speed
    completeMerged <- merge(complete, seColData, by = "col_id", all.x = TRUE)
    completeMerged <- merge(completeMerged, seRowData, by = "row_id", all.x = TRUE)
    
    # Merge with original data using data.table merge
    data_assigned <- merge(data, completeMerged, by = c(cond_entries, cl_entries), all.x = TRUE)
    
    # Order by factor_id using data.table's setorder
    data.table::setorder(data_assigned, factor_id)
    
    bm <- BumpyMatrix::splitAsBumpyMatrix(data_assigned[, allMetadata$data_fields, with = FALSE],
                                          row = data_assigned$row_id,
                                          column = data_assigned$col_id)
    
    # Check if the order of rows/cols are correct
    stopifnot(!is.unsorted(as.numeric(rownames(bm))))
    stopifnot(!is.unsorted(as.numeric(colnames(bm))))
    
    bm
  }
