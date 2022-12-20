#' df_to_assay
#'
#' Convert data.frame with dose-reponse data to an assay object.
#'
#' The 'assay' object is simply a matrix with rownames being the treatment ids,
#' colnames being the ids of the cell lines and values being the DataFrames with
#' dose-response data for given cell lines under given conditions.
#'
#' @param data tibble or data.frame with drug-response data
#' @param data_type string type of data to be returned: all, for untreated conditions only ('untreated')
#' or for treated conditions only ('treated')
#' @param discard_keys a vector of keys that should be discarded
#'
#' @return matrix
#'
#' @export
df_to_assay <-
  function(data,
           data_type = c("all", "treated", "untreated"),
           discard_keys = NULL) {
    # Assertions:
    stopifnot(any(
      inherits(data, "data.frame"),
      inherits(data, "tbl_df"),
      checkmate::test_character(data),
      inherits(data, "DataFrame")
    ))
    checkmate::assert_character(data_type)
    checkmate::assert_character(discard_keys, null.ok = TRUE)
    ####
    data <- S4Vectors::DataFrame(data)
    
    allMetadata <- getMetaData(data, discard_keys = discard_keys)

    seColData <- allMetadata$colData
    cl_entries <- setdiff(colnames(seColData), c("col_id", "name_"))
    seRowData <- allMetadata$rowData
    cond_entries <-
      setdiff(colnames(seRowData), c("row_id", "name_"))
    
    dataCols <- allMetadata$dataCols

    complete <-
      S4Vectors::DataFrame(
        base::expand.grid(
          row_id = seRowData$row_id,
          col_id = seColData$col_id,
          stringsAsFactors = FALSE
        )
      )
    complete <- S4Vectors::merge(S4Vectors::merge(complete, seRowData, by = "row_id"),
                      seColData, by = "col_id")
    complete <- complete[order(complete$col_id, complete$row_id), ]
    complete$factor_id <- seq_len(nrow(complete))
    
    data_assigned <-
      S4Vectors::merge(data, complete, by = c(cond_entries, cl_entries))
    
    by_factor <- lapply(seq_len(nrow(complete)), function(x)
      data_assigned[data_assigned$factor_id == x, dataCols])
    names(by_factor) <- seq_len(nrow(complete))

    stopifnot(nrow(data) == sum(sapply(by_factor, nrow)))
    stopifnot(length(by_factor) == nrow(complete))

    full.set <- by_factor

    dim(full.set) <- c(nrow(seRowData), nrow(seColData))
    dimnames(full.set) <- list(seRowData$name_, seColData$name_)

    #nolint start
      #add NAs for treatments not present in the given assay
      # ---------------
      # removed as it should be added when combining different assays
      #
      # tNotFound <- setdiff(adrug, rownames(full.set))
      # mNotFound <-
      #   matrix(nrow = length(tNotFound), ncol = ncol(full.set))
      # dimnames(mNotFound) <- list(tNotFound, colnames(full.set))
      #
      # final.set <- rbind(mNotFound, full.set)
      #
      # ---------------
    #nolint end
    final.set <- full.set

    if (data_type == "untreated") {
      untreatedConds <-
        .get_untreated_conditions(seRowData)
      return(final.set[rownames(final.set) %in% untreatedConds, , drop = FALSE])
    } else if (data_type == "treated") {
      treatedConds <- .get_treated_conditions(seRowData)
      return(final.set[rownames(final.set) %in% treatedConds, , drop = FALSE])
    } else if (data_type == "all") {
      return(final.set)
    } else {
      stop(sprintf("bad 'data_type': ('%s')", data_type))
    }
  }

#' .get_untreated_conditions
#'
#' Get untreated conditions
#'
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with untreated conditions
#'
.get_untreated_conditions <-
  function(drug_data) {
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    .untreated_tag_patterns <- vapply(get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
    .untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")
    drugnames <- tolower(as.data.frame(drug_data)[, get_identifier("drugname")])
    drug_data[grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }

#' .get_treated_conditions
#'
#' Get treated conditions
#'
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with treated conditions
#'
.get_treated_conditions <-
  function(drug_data) {
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    .untreated_tag_patterns <- vapply(get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
    .untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")
    drugnames <- tolower(as.data.frame(drug_data)[, get_identifier("drugname")])
    drug_data[!grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }
