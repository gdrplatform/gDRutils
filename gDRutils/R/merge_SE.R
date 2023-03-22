#' Merge multiple Summarized Experiments
#'
#' @param SElist named list of Summarized Experiments
#' @param additional_col_name string with the name of the column that will be
#' added to assay data for the distinction of possible duplicated metrics
#' that can arise from multiple projects
#' @param discard_keys character vector of string that will be discarded
#' during creating BumpyMatrix object
#'
#' @return merged SummarizedExperiment object
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
  
  normalized <- merge_assay(SElist = SElist, assay_name = "Normalized", additional_col_name = additional_col_name,
                            discard_keys = discard_keys)
  averaged <- merge_assay(SElist = SElist, assay_name = "Averaged", additional_col_name = additional_col_name,
                          discard_keys = discard_keys)
  metrics <- merge_assay(SElist = SElist, assay_name = "Metrics", additional_col_name = additional_col_name,
                         discard_keys = discard_keys)

  if (!is.null(additional_col_name)) {
    data.table::set(averaged$DT, , intersect(names(averaged$DT),
                                             c(additional_col_name, discard_keys)), NULL)
  }
  data <- split_SE_components(averaged$DT)
  
  data$treatment_md$cId <- NULL
  metadataNames <- identify_unique_se_metadata_fields(SElist)
  
  identifiers <- NULL
  identifiersNames <- "identifiers"
  if (identifiersNames %in% metadataNames) {
    metadataNames <- setdiff(metadataNames, identifiersNames)
    identifiers <- S4Vectors::metadata(SElist[[1]])[identifiersNames]
  }
  
  metadata <- merge_metadata(SElist, metadataNames)
  metadata <- c(metadata, identifiers)
  
  # 2021.11.01 - param checkDimnames was added to Summarized Experiment in Bioc 3.14, 
  # to make it compatible with previous solution we re filtering it 'checkDimnames' if not present
  # TODO: remove once all our envs are with Bioc 3.14
  p_list <-
    list(
      assays = list(
        Normalized = normalized$BM,
        Averaged = averaged$BM,
        Metrics = metrics$BM
      ),
      colData = data$condition_md,
      rowData = data$treatment_md,
      metadata = metadata,
      checkDimnames = FALSE
    )
  av_pnames <-
    names(formals(SummarizedExperiment::SummarizedExperiment))
  f_list <- p_list[intersect(names(p_list), av_pnames)]
  do.call(SummarizedExperiment, f_list)
}


#' Merge assay data
#'
#' @param SElist named list of Summarized Experiments
#' @param assay_name name of the assay that should be extracted and merged
#' @param additional_col_name string of column name that will be
#' added to assay data for the distinction of possible duplicated metrics
#' that can arise from multiple projects
#' @param discard_keys character vector of string that will be discarded
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
    data.table::rbindlist(lapply(names(SElist), function(y)
      convert_se_assay_to_dt(SElist[[y]], assay_name)), fill = TRUE)
  } else {
    data.table::rbindlist(lapply(names(SElist), function(y)
      convert_se_assay_to_dt(SElist[[y]], assay_name)[, eval(additional_col_name) := rep_len(y, .N)]), fill = TRUE)
  }
  

  DT$rId <- DT$cId <- NULL
  BM <- df_to_bm_assay(DT, discard_keys = c(discard_keys, additional_col_name))

  list(DT = DT, BM = BM)
}


#' Identify unique metadata fields from a list of \code{SummarizedExperiment}s
#'
#' @param SElist named list of \code{SummarizedExperiment}s
#'
#' @return character vector of unique names of metadata
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
#' @param SElist named list of \code{SummarizedExperiment}s
#' @param metadata_fields vector of metadata names that will be merged
#'
#' @return list of merged metadata
#' @export
#'
merge_metadata <- function(SElist,
                           metadata_fields) {
  
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_character(metadata_fields)
  
  all_metadata <- lapply(metadata_fields, function(x) {
    do.call(c, lapply(names(SElist), function(SE) {
      meta <- list(S4Vectors::metadata(SElist[[SE]])[[x]])
      names(meta) <- SE
      meta
    }))
  })
  names(all_metadata) <- metadata_fields
  all_metadata
}
