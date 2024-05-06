#' Prettify metric names in flat 'Metrics' assay
#'
#' Map existing column names of a flattened 'Metrics' assay to prettified names.
#'
#' @param x character vector of names to prettify.
#' @param human_readable boolean indicating whether or not to return column names in human readable format.
#' Defaults to \code{FALSE}.
#' @param normalization_type character vector with a specified normalization type.
#' Defaults to \code{c("GR", "RV")}.
#' @keywords identifiers
#'
#' @return character vector of prettified names.
#'
#' @details
#' A common use case for this function is to prettify column names from a flattened version of
#' the \code{"Metrics"} assay.
#' Mode \code{"human_readable" = TRUE} is often used for prettification in the context
#' of front-end applications, whereas \code{"human_readable" = FALSE} is often used for
#' prettification in the context of the command line.
#'
#' @examples 
#' x <- c("CellLineName", "Tissue", "Primary Tissue", "GR_gDR_x_mean", "GR_gDR_xc50", "RV_GDS_x_mean")
#' prettify_flat_metrics(x, human_readable = FALSE)
#'
#' @export
#'
prettify_flat_metrics <- function(x,
                                  human_readable = FALSE,
                                  normalization_type = c("GR", "RV")) {
  
  new_names <- .convert_norm_specific_metrics(x, normalization_type = normalization_type)
  
  if (human_readable) {
    new_names <- .prettify_GDS_columns(new_names)
    new_names <- .prettify_metadata_columns(new_names)
    new_names <- .prettify_metric_columns(new_names)
    new_names <- .prettify_cotreatment_columns(new_names)
    new_names <- gsub("_", " ", new_names)
  }
  
  # gDR is the default name.
  new_names <- gsub("gDR", "", new_names)
  new_names <- gsub("^_+", "", new_names)
  new_names <- gsub("Moa", "MOA", new_names)
  trimws(new_names)
}


####################
# Prettify helpers
####################

#' This function change raw names of metric from long format table into more descriptive
#' names in the wide format table. 
#' It works for metrics: \code{colnames(get_header("metrics_names"))}
#' 
#' @keywords internal
.convert_norm_specific_metrics <- function(x, normalization_type) {
  
  # to do not crush app for unsupported norm type
  normalization_type <- 
    normalization_type[normalization_type %in% rownames(get_header("metrics_names"))]
  
  for (norm in normalization_type) {
    metrics_names <- get_header("metrics_names")[norm, ]
    
    is_norm <- grepl(norm, x)
    
    if (!sum(is_norm)) next # to skip loop below (nothing to change)
    
    # to do not touch combo metrics
    combo_pattern <- paste(c(names(get_combo_score_field_names()), 
                             names(get_combo_excess_field_names())), collapse = "|")
    is_combo <- grepl(combo_pattern, x)
    
    if (!sum(is_norm & !is_combo)) next # to skip loop below (nothing to change)
    
    for (name in names(metrics_names)) {
      replace <- is_norm & grepl(name, x) & !is_combo
      x[replace] <- gsub(name, metrics_names[[name]], gsub(norm, "", x[replace]))
    }
  }
  x
}


#' @keywords internal
.prettify_GDS_columns <- function(cols) {
  # Move the GDS source info to the end as '(GDS)'.
  GDS <- "(GDS)(.*?)_(.*)"
  gsub(GDS, "\\2\\3 (\\1)", cols)
}


#' @keywords internal
.prettify_cotreatment_columns <- function(cols) {
  
  # Replace underscore by space for the Drug/Concentration for co-treatment.
  pattern <- "[0-9]+"
  conc_cotrt <- paste0("^Concentration_", pattern, "$")
  drug_cotrt <- paste0("^", get_env_identifiers("drug", simplify = TRUE), "_", pattern, "$|^drug_.*", pattern,
                       "$|^DrugName_", pattern, "$")
  
  replace <- grepl(paste0(conc_cotrt, "|", drug_cotrt), cols)
  cols[replace] <- gsub("_", " ", cols[replace])
  cols[replace] <- gsub("Name", "", cols[replace])
  cols[replace] <- stringr::str_to_title(gsub("([a-z])([A-Z])", "\\1 \\2", cols[replace]))
  cols
}


#' @keywords internal
.prettify_metric_columns <- function(cols) {
  metric_patterns <- list("E 0" = "E0",
                          "AOC Range" = "AOC within set range",
                          "GRvalue" = "GR value",
                          "RelativeViability" = "Relative Viability",
                          "mean" = "Mean Viability")
  
  for (i in names(metric_patterns)) {
    cols <- gsub(i, metric_patterns[i], cols)
  }
  
  cols
}


#' @keywords internal
.prettify_metadata_columns <- function(cols) {
  
  # prettifying formatting
  
  prettified_cols <- gsub("gDR", "", cols)
  prettified_cols <- gsub("_", " ", prettified_cols)
  prettified_cols <- tools::toTitleCase(prettified_cols)
  # exception
  prettified_cols <- gsub("Hsa ", "HSA ", prettified_cols)
  # adding space between words like “ReferenceDivisionTime”
  prettified_cols <- gsub("([a-z])([A-Z])", "\\1 \\2", prettified_cols)
  
  trimws(prettified_cols)
}
