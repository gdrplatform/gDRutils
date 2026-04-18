#' Merge multiple MultiAssayExperiment objects
#'
#' @param MAElist Named list of MultiAssayExperiment objects.
#' @param additional_col_name String with the name of the column that will be
#' added to assay data for the distinction of possible duplicated metrics
#' that can arise from multiple projects.
#' @param discard_keys Character vector of strings that will be discarded
#' during creating BumpyMatrix object.
#' @param title String specifying the final DataSetDB title. If NULL, auto-generates.
#' @param description String specifying the final DataSetDB description. If NULL, auto-generates.
#' @param source_name String specifying the standard DSDB source name. If NULL, auto-detects or uses "merged_analysis".
#' @param source_id String specifying the unique DSDB source ID. If NULL, uses "merged_dataset".
#' @keywords SE_operators
#'
#' @examples
#' mae1 <- get_synthetic_data("finalMAE_combo_2dose_nonoise.qs2")
#' mae2 <- get_synthetic_data("finalMAE_combo_2dose_nonoise.qs2")
#' merge_MAE(list(mae1 = mae1, mae2 = mae2), title = "Test", description = "Test MAE")
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
                                       "conc_2"),
                      title = NULL,
                      description = NULL,
                      source_name = NULL,
                      source_id = NULL) {
  
  checkmate::assert_list(MAElist, types = "MultiAssayExperiment")
  checkmate::assert_string(title, null.ok = TRUE)
  checkmate::assert_string(description, null.ok = TRUE)
  checkmate::assert_string(source_name, null.ok = TRUE)
  checkmate::assert_string(source_id, null.ok = TRUE)
  
  experiments <- unique(unlist(lapply(MAElist, names)))
  
  merged_SE_assays <- lapply(experiments, function(exp_name) {
    exp_list <- lapply(MAElist, function(mae) {
      if (exp_name %in% names(mae)) mae[[exp_name]] else NULL
    })
    exp_list <- exp_list[!vapply(exp_list, is.null, FUN.VALUE = logical(1))]
    merge_SE(exp_list)
  })
  names(merged_SE_assays) <- experiments
  
  mae_names <- names(MAElist)
  if (is.null(mae_names) || all(trimws(mae_names) == "")) {
    mae_names <- paste0("Dataset_", seq_along(MAElist))
  }
  
  all_sources <- list()
  original_titles <- c()
  
  for (mae in MAElist) {
    for (exp in names(mae)) {
      meta <- as.list(S4Vectors::metadata(mae[[exp]])$experiment_metadata)
      if (length(meta) > 0) {
        if (is.list(meta$sources)) all_sources <- c(all_sources, meta$sources)
        if (!is.null(meta$title)) original_titles <- c(original_titles, meta$title)
      }
    }
  }
  

  if (is.null(title)) {
    title <- sprintf("Merged MAE: %s", paste(mae_names, collapse = " + "))
  }
  
  if (is.null(description)) {
    description <- sprintf("Synthetically merged dataset originating from: %s.", paste(mae_names, collapse = ", "))
    unique_titles <- unique(original_titles)
    if (length(unique_titles) > 0) {
      description <- paste0(description, " Original Titles: [", paste(unique_titles, collapse = " | "), "]")
    }
  }
  
  if (is.null(source_name)) {
    if (length(all_sources) > 0) {
      unique_names <- unique(vapply(all_sources, function(s) {
        if (!is.null(s$name)) s$name else "unknown"
      }, character(1)))
      
      source_name <- if (length(unique_names) == 1 && unique_names[1] != "unknown") {
        unique_names[1]
      } else {
        "merged_analysis"
      }
    } else {
      source_name <- "merged_analysis"
    }
  }
  
  if (is.null(source_id)) {
    source_id <- "merged_dataset"
  }
  
  synthetic_experiment_metadata <- list(
    title = title,
    description = description,
    experimentalist = Sys.info()[["user"]],
    sources = list(list(name = source_name, id = source_id))
  )
  
  for (i in seq_along(merged_SE_assays)) {
    meta_list <- as.list(S4Vectors::metadata(merged_SE_assays[[i]]))
    meta_list$experiment_metadata <- synthetic_experiment_metadata
    S4Vectors::metadata(merged_SE_assays[[i]]) <- meta_list
  }
  
  base_metadata <- as.list(S4Vectors::metadata(MAElist[[1]]))
  if (length(base_metadata) == 0) base_metadata <- list()
  
  if (!is.null(base_metadata$.internal$DataSetDB$dataset)) {
    ds_meta <- as.list(base_metadata$.internal$DataSetDB$dataset)
    ds_meta$title <- synthetic_experiment_metadata$title
    ds_meta$description <- synthetic_experiment_metadata$description
    ds_meta$sources <- synthetic_experiment_metadata$sources
    
    internal_meta <- as.list(base_metadata$.internal)
    internal_meta$DataSetDB <- as.list(internal_meta$DataSetDB)
    internal_meta$DataSetDB$dataset <- ds_meta
    
    base_metadata$.internal <- internal_meta
  }
  
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(merged_SE_assays),
    metadata = base_metadata
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
#' se1 <- get_synthetic_data("finalMAE_small.qs2")[[1]]
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
#' mae <- get_synthetic_data("finalMAE_combo_2dose_nonoise.qs2")
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
  
  drug_cols <- unlist(get_env_identifiers(c("drug", "drug2", "drug3"), simplify = FALSE))
  existing_drug_cols <- intersect(drug_cols, names(DT))
  DT[, (existing_drug_cols) := lapply(.SD, remove_drug_batch), .SDcols = existing_drug_cols]
  
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
#' mae <- get_synthetic_data("finalMAE_small.qs2")
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
#' mae <- get_synthetic_data("finalMAE_small.qs2")
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
    
    if (x %in% c("experiment_metadata", ".internal")) {
      
      valid_metas <- lapply(SElist, function(se) S4Vectors::metadata(se)[[x]])
      valid_metas <- valid_metas[!vapply(valid_metas, is.null, FUN.VALUE = logical(1))]
      
      if (length(valid_metas) == 0) return(list())
      
      if (x == "experiment_metadata") {
        synth <- as.list(valid_metas[[1]]) 
        
        all_sources <- list()
        for (vm in valid_metas) {
          vm_list <- as.list(vm)
          if (is.list(vm_list$sources)) all_sources <- c(all_sources, vm_list$sources)
        }
        
        if (length(all_sources) > 0) {
          unique_names <- unique(vapply(all_sources, function(s) {
            if (!is.null(s$name)) s$name else "unknown"
          }, character(1)))
          
          std_name <- if (length(unique_names) == 1 && unique_names[1] != "unknown") {
            unique_names[1]
          } else {
            "merged_analysis"
          }
          
          synth$sources <- list(list(name = std_name, id = "merged_dataset"))
        } else {
          synth$sources <- list()
        }
        
        return(synth)
      }
      
      return(as.list(valid_metas[[1]]))
    }
    
    do.call(c, lapply(names(SElist), function(SE) {
      meta <- list(S4Vectors::metadata(SElist[[SE]])[[x]])
      names(meta) <- SE
      meta
    }))
  })
  
  names(all_metadata) <- metadata_fields
  all_metadata
}
