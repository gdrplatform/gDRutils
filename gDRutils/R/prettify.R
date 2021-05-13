#' Prettify metric names in flat 'Metrics' assay
#'
#' Map existing column names of a flattened 'Metrics' assay to prettified names.
#'
#' @param x character vector of names to prettify.
#' @param human_readable boolean indicating whether or not to return column names in human readable format.
#' Defaults to \code{FALSE}.
#' @param normalization_type a character with a specified normalization type.
#' Defaults to \code{c("GR", "RV")}.
#'
#' @return character vector of prettified names.
#'
#' @details 
#' A common use case for this function is to prettify column names from a flattened version of 
#' the \code{"Metrics"} assay.
#' Mode \code{"human_readable = TRUE"} is often used for prettification in the context
#' of front-end applications, whereas \code{"human_readable" = FALSE} is often used for 
#' prettification in the context of the command line.
#'
#' @export
#'
prettify_flat_metrics <- function(x, 
                                  human_readable = FALSE, 
                                  normalization_type = c("GR", "RV")) {

  new_names <- .convert_norm_specific_metrics(x, normalization_type = normalization_type)

  # gDR is the default name.
  new_names <- gsub("gDR", "", new_names)

  if (human_readable) {
    new_names <- .prettify_GDS_columns(new_names)
    new_names <- .prettify_metadata_columns(new_names)
    new_names <- .prettify_metric_columns(new_names)

    # replace underscore by space for the remaining metrics
    new_names[is_metric] <- gsub("_", " ", new_names[is_metric])

    new_names <- .prettify_cotreatment_columns(new_names)
  }

  new_names
}


####################
# Prettify helpers
####################

.convert_norm_specific_metrics <- function(x, normalization_type) {
  for (norm in normalization_type) {
    metrics_names <- get_header("metrics_names")[norm, ]
    if (length(metrics_names) == 0L) {
      stop(sprintf("missing normalization type: '%s' from header: 'metrics_names'", norm))
    }

    norm_pattern <- paste0("^", norm, "_")
    is_norm <- grepl(norm_pattern, x)
    for (name in names(metrics_names)) {
      replace <- is_norm & grepl(norm, x)
      x[replace] <- gsub(name, metrics_names[name], x[replace])
    }
  }
  x
}


#' @keywords internal
.prettify_GDS_columns <- function(cols) {
  # Move the GDS source info to the end as '(GDS)'.
  GDS <- "(GDS)(.*?)_(.*)"
  cols <- gsub(GDS, "\\2\\3 (\\1)", cols)
}


#' @keywords internal
.prettify_cotreatment_columns <- function(cols) {
    # Replace underscore by space for the Drug/Concentration for co-treatment.
    pattern <- "[0-9]+"
    conc_cotrt <- paste0("^Concentration_", pattern, "$")
    drug_cotrt <- paste0("^", get_identifier("drug"), "_", pattern, "$|^Drug_", pattern, "$")

    replace <- grepl(paste0(conc_cotrt, "|", drug_cotrt), cols)
    cols[replace] <- gsub("_", " ", cols[replace])
}


#' @keywords internal
.prettify_metric_columns <- function(cols) {
  metric_patterns <- list("E_0" = "E0",
                          "AOC_range" = "AOC within set range",
                          "GRvalue" = "GR value",
                          "RelativeViability" = "Relative Viability",
                          "MeanViability" = "Mean Viability")
  # TODO: Need to account for the '_' in '_MeanViability'

  for (i in names(metric_patterns)) {
    cols <- gsub(i, metric_patterns[i], cols)
  }

  cols
}


#' @keywords internal
.prettify_metadata_columns <- function(cols) {
  metadata <- c("Cell line", 
            "Primary Tissue", 
            "Subtype",
            "Parental cell line",
            "Reference cell division time", 
            "cell division time",
            "Drug", 
            "Drug MOA",
            "Nbre of tested conc.", 
            "Highest log10(conc.)")

  # TODO: Eventually, these identifiers can come from the SE object itself.
  names(metadata) <- c(get_identifier("cellline_name"), # CellLineName
            get_identifier("cellline_tissue"), # Tissue
            get_identifier("cellline_subtype"), # subtype
            get_identifier("cellline_parental_identifier"), # parental_identifier
            "ReferenceDivisionTime", 
            "DivisionTime",
            get_identifier("drugname"), # DrugName
            get_identifier("drug_moa"), # drug_moa
            "N_conc", 
            "maxlog10Concentration")
  cols[cols %in% names(metadata)] <- metadata[cols[cols %in% names(metadata)]]
  cols
}
