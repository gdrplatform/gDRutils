
#' df_to_bm_assay
#'
#' Convert data.table with dose-reponse data into a BumpyMatrix assay.
#'
#' The 'assay' is simply a BumpyMatrix object with rownames being the treatment ids,
#' colnames being the ids of the cell lines and values with
#' dose-response data for given cell lines under given conditions.
#'
#' @param data data.table with drug-response data
#' @param discard_keys a vector of keys that should be discarded
#' @keywords convert_to_bm
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
    
    data <- methods::as(data, "DataFrame")
    allMetadata <- split_SE_components(data, nested_keys = discard_keys)
    
    seColData <- allMetadata$condition_md
    seColData$col_id <- seq_along(rownames(seColData))
    cl_entries <- setdiff(colnames(seColData), c("col_id", "name_"))
    seRowData <- allMetadata$treatment_md
    seRowData$row_id <- seq_along(rownames(seRowData))
    cond_entries <- setdiff(colnames(seRowData), c("row_id", "name_"))
    
    complete <- S4Vectors::DataFrame(
      expand.grid(col_id = seColData$col_id, row_id = seRowData$row_id, stringsAsFactors = FALSE)
    )
    complete$factor_id <- seq_len(nrow(complete))
    completeMerged <- merge(merge(complete, seColData, by = "col_id"), seRowData, by = "row_id")
    
    data_assigned <- merge(data, completeMerged, by = c(cond_entries, cl_entries))
    data_assigned <- data_assigned[order(data_assigned$factor_id), ]
    bm <- BumpyMatrix::splitAsBumpyMatrix(data_assigned[, allMetadata$data_fields, drop = FALSE], 
                                          row = data_assigned$row_id,
                                          column = data_assigned$col_id)
    
    # Check if the orderd of rows/cols are correct
    stopifnot(!is.unsorted(as.numeric(rownames(bm))))
    stopifnot(!is.unsorted(as.numeric(colnames(bm))))
    
    bm
  }
