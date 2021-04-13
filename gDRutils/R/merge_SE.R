#' Merge multiple Summarized Experiments
#'
#' @param SE list of Summarized Experiments
#'
#' @return Summarized Experiment
#' @export
#'
merge_SE <- function(SE) {
  checkmate::assert_list(SE, types = "SummarizedExperiment", min.len = 1)
  names(SE) <- 
    vapply(SE, function(x) S4Vectors::metadata(x)$experiment_metadata$qcs_id, "QCS")
  normalizedDT <- data.table::rbindlist(lapply(SE, gDRutils::assay_to_dt, "Normalized"), fill = TRUE)
  normalizedDT$rId <- NULL
  normalizedDT$cId <- NULL
  normalizedBM <- gDRutils::df_to_bm_assay(normalizedDT)
  
  averagedDT <- data.table::rbindlist(lapply(SE, gDRutils::assay_to_dt, "Averaged"), fill = TRUE)
  averagedDT$rId <- NULL
  averagedDT$cId <- NULL
  averagedBM <- gDRutils::df_to_bm_assay(averagedDT)
  
  metricsDT <- data.table::rbindlist(lapply(SE, function(x) gDRutils::assay_to_dt(x, "Metrics", merge_metrics = TRUE)), fill = TRUE)
  metricsDT$rId <- NULL
  metricsDT$cId <- NULL
  
  RVcols <- grep("RV|^EC50|IC50|^E_", colnames(metricsDT), value = TRUE)
  GRcols <- grep("GR|GEC", colnames(metricsDT), value = TRUE)
  commonCols <- setdiff(colnames(metricsDT), c(RVcols, GRcols))
  metricsRV <- metricsDT[, RVcols, with = FALSE]
  newCols <- c("x_mean",
               "x_AOC",
               "x_AOC_range",
               "xc50",
               "x_max",
               "c50",
               "x_inf",
               "x_0",
               "h",
               "r2",
               "x_sd_avg",
               "fit_type")
  colnames(metricsRV) <- newCols
  metricsRV$dr_metric <- "RV"
  
  metricsGR <- metricsDT[, GRcols, with = FALSE]
  colnames(metricsGR) <- newCols
  metricsGR$dr_metric <- "GR"
  
  metricsRV <- cbind(metricsRV, metricsDT[, commonCols, with = FALSE])
  metricsGR <- cbind(metricsGR, metricsDT[, commonCols, with = FALSE])
  
  mergedMetrics <- as.data.frame(mapply(function(x,y) rbind(x,y), metricsRV, metricsGR))
  
  numCols <- unname(gDRutils::get_header("metrics_results"))
  numCols <- numCols[numCols != "fit_type"]
  mergedMetrics[, which(colnames(mergedMetrics) %in% numCols)] <- 
    as.data.frame(lapply(mergedMetrics[, which(names(mergedMetrics) %in% numCols)], as.numeric))
  
  metricsBM <- gDRutils::df_to_bm_assay(mergedMetrics, discard_keys = "dr_metric")
  SEdata <- gDR::split_SE_components(metricsDT)
  
  
  if (any(dim(normalizedBM) != dim(metricsBM))) {
    normalizedSplit <- gDR::split_SE_components(normalizedDT)
    metricsSplit <- gDR::split_SE_components(metricsDT)
    treatmentCols <- setdiff(c(colnames(metricsSplit$treatment_md),
                               colnames(normalizedSplit$treatment_md)),
                             c("rId", "cId"))
    
    conditionCols <- setdiff(c(colnames(metricsSplit$condition_md),
                               colnames(normalizedSplit$condition_md)),
                             c("rId", "cId"))
    
    normalizedBM <- normalizedBM[match(metricsSplit$treatment_md[, treatmentCols], 
                                       normalizedSplit$treatment_md[, treatmentCols]),
                                 match(metricsSplit$condition_md[, conditionCols],
                                       normalizedSplit$condition_md[, conditionCols])]
    
  }
  
  SEdata$treatment_md$cId <- NULL
  
  initialSE <- SummarizedExperiment::SummarizedExperiment(assays = list(Normalized = normalizedBM,
                                                                        Averaged = averagedBM,
                                                                        Metrics = metricsBM),
                                                          colData = SEdata$condition_md,
                                                          rowData = SEdata$treatment_md)
  
  initialSE
}