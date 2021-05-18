#' Merge multiple Summarized Experiments
#'
#' @param SElist named list of Summarized Experiments
#' @param additional_col_name a string with the name of the column that will be
#' added to assay data for the distinction of possible duplicated metrics
#' that can arise from multiple projects
#' @param discard_keys a character vector of string that will be discarded
#' during creating BumpyMatrix object
#'
#' @return merged Summarized Experiment
#' @export
#'
merge_SE <- function(SElist,
                     additional_col_name = NULL,
                     discard_keys = c("normalization_type",
                                      "fit_source",
                                      "Metrics_rownames")) {
  
  assertthat::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_string(additional_col_name, null.ok = TRUE)
  checkmate::assert_character(discard_keys, null.ok = TRUE)
  
  normalized <- merge_assay(SElist = SElist, assay_name = "Normalized", additional_col_name = additional_col_name)
  averaged <- merge_assay(SElist = SElist, assay_name = "Averaged", additional_col_name = additional_col_name)
  metrics <- merge_assay(SElist = SElist, assay_name = "Metrics", additional_col_name = additional_col_name,
                         discard_keys = c("normalization_type",
                                          "fit_source",
                                          "Metrics_rownames"))
  SEdata <- if (is.null(additional_col_name)) {
    split_SE_components(averaged$DT)
  } else {
    averaged$DT[[additional_col_name]] <- NULL
    split_SE_components(averaged$DT)
  }
  
  SEdata$treatment_md$cId <- NULL
  metadataNames <- unique(unlist(lapply(SElist, function(x) {
    names(S4Vectors::metadata(x))
  })))
  metadataSE <- vector("list", length(metadataNames))
  names(metadataSE) <- metadataNames
  
  metadataSEAll <- lapply(names(metadataSE), function(x) {
    mergedMetadata <- NULL
    do.call(c, lapply(names(SElist), function(SE) {
      metaObj <- list(S4Vectors::metadata(SElist[[SE]])[[x]])
      names(metaObj) <- SE
      metaObj
    }))
  })
  names(metadataSEAll) <- metadataNames
  initialSE <- SummarizedExperiment::SummarizedExperiment(assays = list(Normalized = normalized$BM,
                                                                        Averaged = averaged$BM,
                                                                        Metrics = metrics$BM),
                                                          colData = SEdata$condition_md,
                                                          rowData = SEdata$treatment_md,
                                                          metadata = metadataSE)
  initialSE
}



#' Merge assay data
#'
#' @param SElist named list of Summarized Experiments
#' @param assay_name name of the assay that should be extracted and merged
#' @param additional_col_name a string with the name of the column that will be
#' added to assay data for the distinction of possible duplicated metrics
#' that can arise from multiple projects
#' @param discard_keys a character vector of string that will be discarded
#' during creating BumpyMatrix object
#'
#' @return BumpyMatrix or list with data.table + BumpyMatrix
#' @export
#'
#' @examples
merge_assay <- function(SElist,
                        assay_name,
                        discard_keys = NULL,
                        additional_col_name = NULL) {
  
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_true(all(unlist(lapply(SElist, function(x) {
    assay_name %in% names(SummarizedExperiment::assays(x))
  }))))
  checkmate::assert_string(additional_col_name, null.ok = TRUE)
  checkmate::assert_character(discard_keys, null.ok = TRUE)

  DT <- if (is.null(additional_col_name)) {
    data.table::rbindlist(lapply(names(SElist), function(x)
      convert_se_assay_to_dt(SElist[[x]], assay_name)), fill = TRUE)
  } else {
    data.table::rbindlist(lapply(names(SElist), function(x)
      convert_se_assay_to_dt(SElist[[x]], assay_name)[, eval(additional_col_name) := rep_len(x, .N)]), fill = TRUE)
  }
  DT$rId <- NULL
  DT$cId <- NULL
  BM <- df_to_bm_assay(DT, discard_keys = c(discard_keys, additional_col_name))
  list(DT = DT, BM = BM)
}
