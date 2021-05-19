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
                     additional_col_name = "data_source",
                     discard_keys = c("normalization_type",
                                      "fit_source",
                                      "Metrics_rownames")) {
  
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_string(additional_col_name, null.ok = TRUE)
  checkmate::assert_character(discard_keys, null.ok = TRUE)
  
  normalized <- merge_assay(SElist = SElist, assay_name = "Normalized", additional_col_name = additional_col_name)
  averaged <- merge_assay(SElist = SElist, assay_name = "Averaged", additional_col_name = additional_col_name)
  metrics <- merge_assay(SElist = SElist, assay_name = "Metrics", additional_col_name = additional_col_name,
                         discard_keys = discard_keys)
  SEdata <- if (is.null(additional_col_name)) {
    split_SE_components(averaged$DT)
  } else {
    averaged$DT[[additional_col_name]] <- NULL
    split_SE_components(averaged$DT)
  }
  
  SEdata$treatment_md$cId <- NULL
  metadataNames <- identify_unique_se_metadata_fields(SElist)
  metadataSE <- merge_metadata(SElist, metadataNames)
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(Normalized = normalized$BM,
                                                                        Averaged = averaged$BM,
                                                                        Metrics = metrics$BM),
                                                          colData = SEdata$condition_md,
                                                          rowData = SEdata$treatment_md,
                                                          metadata = metadataSE)
  se
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
merge_assay <- function(SElist,
                        assay_name,
                        additional_col_name = "data_source",
                        discard_keys = NULL) {
  
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  lapply(SElist, function(x) validate_se_assay_name(x, assay_name))
  checkmate::assert_string(additional_col_name, null.ok = TRUE)
  checkmate::assert_character(discard_keys, null.ok = TRUE)

  DT <- if (is.null(additional_col_name)) {
    data.table::rbindlist(lapply(names(SElist), function(x)
      convert_se_assay_to_dt(SElist[[x]], assay_name)), fill = TRUE)
  } else {
    data.table::rbindlist(lapply(names(SElist), function(x)
      convert_se_assay_to_dt(SElist[[x]], assay_name)[, eval(additional_col_name) := rep_len(x, .N)]), fill = TRUE)
  }
  DT$rId <- DT$cId <- NULL
  BM <- df_to_bm_assay(DT, discard_keys = c(discard_keys, additional_col_name))
  list(DT = DT, BM = BM)
}


#' Indetify unique se metadata fields from the list of SEs
#'
#' @param SElist named list of Summarized Experiments
#'
#' @return vector of unique names of metadata
#' @export
#'
identify_unique_se_metadata_fields <- function(SElist) {
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  
  unique(unlist(lapply(SElist, function(x) {
    names(S4Vectors::metadata(x))
  })))
}

#' Merge metadata
#'
#' @param SElist named list of Summarized Experiments
#' @param metadata_fields vector of metadata names that will be merged
#'
#' @return list of merged metadata
#' @export
#'
merge_metadata <- function(SElist,
                           metadata_fields) {
  
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_character(metadata_fields)
  
  metadataSEAll <- lapply(metadata_fields, function(x) {
    do.call(c, lapply(names(SElist), function(SE) {
      metaObj <- list(S4Vectors::metadata(SElist[[SE]])[[x]])
      names(metaObj) <- SE
      metaObj
    }))
  })
  names(metadataSEAll) <- metadata_fields
  metadataSEAll
}