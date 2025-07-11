#' Merge multiple MultiAssayExperiment objects
#'
#' @param MAElist Named list of MultiAssayExperiment objects.
#' @param additional_col_name String with the name of the column that will be
#' added to assay data for the distinction of possible duplicated metrics
#' that can arise from multiple projects.
#' @param discard_keys Character vector of strings that will be discarded
#' during creating BumpyMatrix object.
#' @keywords SE_operators
#'
#' @examples
#' mae1 <- get_synthetic_data("finalMAE_combo_2dose_nonoise")
#' mae2 <- get_synthetic_data("finalMAE_combo_2dose_nonoise")
#' merge_MAE(list(mae1 = mae1, mae2 = mae2))
#'
#' @return Merged MultiAssayExperiment object.
#' @export
merge_MAE <- function(MAElist,
                      additional_col_name = "data_source",
                      discard_keys = c("normalization_type",
                                       "fit_source",
                                       "record_id",
                                       "isDay0",
                                       "swap_sa",
                                       "control_type",
                                       "iso_level",
                                       "conc_1",
                                       "conc_2")) {
  
  checkmate::assert_list(MAElist, types = "MultiAssayExperiment")
  
  experiments <- unique(unlist(lapply(MAElist, names)))
  
  merged_SE_assays <- lapply(experiments, function(exp_name) {
    exp_list <- lapply(MAElist, function(mae) {
      if (exp_name %in% names(mae)) {
        mae[[exp_name]]
      } else {
        NULL
      }
    })
    exp_list <- exp_list[!vapply(exp_list, is.null, FUN.VALUE = logical(1))]
    merge_SE(exp_list)
  })
  names(merged_SE_assays) <- experiments
  
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(merged_SE_assays),
    metadata = Reduce(c, lapply(MAElist, S4Vectors::metadata))
  )
}

#' Merge multiple Summarized Experiments
#'
#' @param SElist named list of Summarized Experiments
#' @param additional_col_name string with the name of the column that will be
#' added to assay data for the distinction of possible duplicated metrics
#' that can arise from multiple projects
#' @param discard_keys character vector of string that will be discarded
#' during creating BumpyMatrix object
#' @keywords SE_operators
#'
#' @examples
#' se1 <- get_synthetic_data("finalMAE_small")[[1]]
#' merge_SE(list(se1 = se1, se2 = se1))
#'
#' @return merged SummarizedExperiment object
#' @export
#'
merge_SE <- function(SElist,
                     additional_col_name = "data_source",
                     discard_keys = c("normalization_type",
                                      "fit_source",
                                      "record_id",
                                      "isDay0",
                                      "swap_sa",
                                      "control_type",
                                      "iso_level",
                                      "conc_1",
                                      "conc_2")) {
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_string(additional_col_name, null.ok = TRUE)
  checkmate::assert_character(discard_keys, null.ok = TRUE)
  
  SE_identifiers <- unique(lapply(SElist, get_SE_identifiers))[[1]]
  lapply(names(SE_identifiers), function(x) {
    set_env_identifier(x, SE_identifiers[[x]])
  })
  
  discard_keys <- c(discard_keys, unique(unlist(
    lapply(SElist, get_SE_identifiers,
           c("barcode",
             "concentration",
             "concentration2"),
           simplify = FALSE))))
  se_assays <- unique(unlist(lapply(SElist,
                                    SummarizedExperiment::assayNames)))
  merged_assays <- lapply(se_assays, function(x) {
    merge_assay(SElist = SElist,
                assay_name = x,
                additional_col_name = additional_col_name,
                discard_keys = discard_keys)
  })
  
  names(merged_assays) <- se_assays
  
  if (!is.null(additional_col_name)) {
    data.table::set(merged_assays$Averaged$DT, ,
                    intersect(names(merged_assays$Averaged$DT),
                              c(additional_col_name, discard_keys)), NULL)
  }
  data <- split_SE_components(merged_assays$Averaged$DT,
                              nested_keys = intersect(additional_col_name, names(merged_assays$Averaged$DT)))
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
  
  assays <- lapply(
    merged_assays,
    FUN = function(x) {
      bm_assay <- x[["BM"]]
      colnames(bm_assay) <- rownames(data$condition_md)
      rownames(bm_assay) <- rownames(data$treatment_md)
      bm_assay
    }
  )
  
  p_list <-
    list(
      assays = assays,
      colData = data$condition_md,
      rowData = data$treatment_md,
      metadata = metadata
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
#' @keywords SE_operators
#'
#' @return BumpyMatrix or list with data.table + BumpyMatrix
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_combo_2dose_nonoise")
#' 
#' listSE <- list(
#'   combo1 = mae[[1]], 
#'   sa = mae[[2]]
#' )
#' merge_assay(listSE, "Normalized")
#' 
#' @export
#'
merge_assay <- function(SElist,
                        assay_name,
                        additional_col_name = "data_source",
                        discard_keys = NULL) {
  
  checkmate::assert_list(SElist, types = "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_string(additional_col_name, null.ok = TRUE)
  checkmate::assert_character(discard_keys, null.ok = TRUE)
  
  SElist <- lapply(SElist, function(x) {
    if (assay_name %in% SummarizedExperiment::assayNames(x)) {
      x
    } else {
      SummarizedExperiment::assay(x, assay_name) <-
        BumpyMatrix::splitAsBumpyMatrix(S4Vectors::DataFrame(x = rep(NA, prod(dim(x)))),
                                        row = rownames(x), column = colnames(x))
      x
    }
  })
  
  DT <- data.table::rbindlist(lapply(stats::setNames(names(SElist),
                                                     names(SElist)),
                                     function(y) {
                                       convert_se_assay_to_dt(SElist[[y]], assay_name)
                                     }),  fill = TRUE, idcol = additional_col_name)
  
  
  DT$rId <- DT$cId <- NULL
  discard_keys <- intersect(names(DT), c(discard_keys, additional_col_name))
  BM <- df_to_bm_assay(DT, discard_keys = discard_keys)
  list(DT = DT, BM = BM)
}


#' Identify unique metadata fields from a list of \code{SummarizedExperiment}s
#'
#' @param SElist named list of \code{SummarizedExperiment}s
#' @keywords SE_operators
#'
#' @return character vector of unique names of metadata
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' SElist <- list(
#'   se, 
#'   se
#' )
#' identify_unique_se_metadata_fields(SElist)
#' 
#' @export
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
#' @keywords SE_operators
#'
#' @return list of merged metadata
#' 
#' @examples 
#' mae <- get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' listSE <- list(
#'   se, 
#'   se
#' )
#' metadata_fields <- identify_unique_se_metadata_fields(listSE)
#' merge_metadata(listSE, metadata_fields)
#' 
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
