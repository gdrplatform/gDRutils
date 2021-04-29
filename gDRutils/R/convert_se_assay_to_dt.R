#' Convert a SummarizedExperiment assay to a long data.table
#'
#' Convert an assay within a \linkS4class{SummarizedExperiment} object to a long data.table.
#'
#' @param se A \linkS4class{SummarizedExperiment} object holding raw and/or processed dose-response data in its assays.
#' @param assay_name String of name of the assay to transform within the \code{se}.
#' @param include_metadata Boolean indicating whether or not to include \code{rowData(se)}
#' and \code{colData(se)} in the returned data.table.
#' Defaults to \code{TRUE}.
#' @param retain_nested_rownames Boolean indicating whether or not to retain the rownames 
#' nested within a \code{BumpyMatrix} assay.
#' Defaults to \code{FALSE}.
#' If the \code{assay_name} is not of the \code{BumpyMatrix} class, this argument's value is ignored.
#' If \code{TRUE}, the resulting column in the data.table will be named as \code{"<assay_name>_rownames"}.
#'
#' @return data.table representation of the data in \code{assay_name}.
#'
#' @details NOTE: to extract information about 'Control' data, simply call the
#' function with the name of the assay holding data on controls.
#' To extract the reference data in to same format as 'Averaged' use \code{convert_se_ref_assay_to_dt}.
#'
#' @export
#'
convert_se_assay_to_dt <- function(se,
                                   assay_name,
                                   include_metadata = TRUE,
                                   retain_nested_rownames = FALSE) {

  # Assertions.
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  checkmate::assert_flag(include_metadata)
  
  validate_se_assay_name(se, assay_name)

  dt <- .convert_se_assay_to_dt(se, assay_name, retain_nested_rownames = retain_nested_rownames)

  if (nrow(dt) == 0L) {
    return(dt) # TODO: Should this return something else?
  }

  if (include_metadata) {
    rData <- SummarizedExperiment::rowData(se)
    rData$rId <- rownames(rData)

    cData <- SummarizedExperiment::colData(se)
    cData$cId <- rownames(cData)

    ids <- expand.grid(rData$rId, cData$cId)
    colnames(ids) <- c("rId", "cId")
    ids[] <- lapply(ids, as.character)

    annotations <- merge(ids, rData, by = "rId", all.x = TRUE)
    annotations <- merge(annotations, cData, by = "cId", all.x = TRUE)

    dt <- merge(dt, annotations, by = c("rId", "cId"), all.x = TRUE)
  }

  data.table::as.data.table(dt)
}


#' Convert assay data into data.table.
#' @return data.table of assay data.
#' @keywords internal
#' @noRd
#'
.convert_se_assay_to_dt <- function(se, assay_name, retain_nested_rownames) {
 
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay_name)
  
  object <- SummarizedExperiment::assays(se)[[assay_name]]
  checkmate::assert_true(inherits(object, "BumpyDataFrameMatrix") || inherits(object, "matrix"))

  if (methods::is(object, "BumpyDataFrameMatrix")) {
    as_df <- BumpyMatrix::unsplitAsDataFrame(object, row.field = "rId", column.field = "cId")
    # Retain nested rownames.
    if (retain_nested_rownames) {
      as_df[[paste0(assay_name, "_rownames")]] <- rownames(as_df)
    }

  } else if (methods::is(object, "matrix")) {
    first <- object[1, 1][[1]]
    if (is.numeric(first)) {
      as_df <- reshape2::melt(object, varnames = c("rId", "cId"), value.name = assay_name)
    } else {
      stop(sprintf("matrix with nested objects of class '%s' is not supported", class(first)))
    }
    as_df
  } else {
    stop(sprintf("assay of class '%s' is not supported", class(object)))
  }
  data.table::as.data.table(as_df)
}


########################
# Convenience functions
########################

#' Convert the reference values from a SummarizedExperiment assay to a long data.table
#'
#' Transform the Ref[RelativeViability/GRvalue] within a \linkS4class{SummarizedExperiment} object to a long data.table.
#' Clean up the column names and add columns to match the format of the data.table from the \code{'Averaged'} assay.
#'
#' @param se A \linkS4class{SummarizedExperiment} object holding reference data in its assays.
#' @param ref_relative_viability_assay String of the name of the assay in the \code{se} 
#' holding the reference relative viability data.
#' @param ref_gr_value_assay String of the name of the assay in the \code{se} holding the reference GR value data.
#'
#' @return data.table representation of the reference data.
#'
#' @details This is a convenience function to massage the reference data into the same format as the \code{"Averaged"}
#' assay data.
#'
#' @export
#'
convert_se_ref_assay_to_dt <- function(se,
                                       ref_relative_viability_assay = "RefRelativeViability",
                                       ref_gr_value_assay = "RefGRvalue") {

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(ref_relative_viability_assay)
  checkmate::assert_string(ref_gr_value_assay)

  rv <- convert_se_assay_to_dt(se, ref_relative_viability_assay, include_metadata = TRUE)
  colnames(rv)[colnames(rv) == ref_relative_viability_assay] <- "RelativeViability"
  rv$std_RelativeViability <- NA

  gr <- convert_se_assay_to_dt(se, ref_gr_value_assay)
  colnames(gr)[colnames(gr) == ref_gr_value_assay] <- "GRvalue"
  gr$std_GRvalue <- NA

  dt <- merge(rv, gr, all = TRUE, by = intersect(names(rv), names(gr)))

  # Fill primary drug with 'untreated_tag'.
  dt$Concentration <- 0
  untreated <- get_SE_identifiers(se, "untreated_tag")[1]
  dt[, get_SE_identifiers(se, "drug")] <- untreated
  dt[, get_SE_identifiers(se, "drugname")] <- untreated
  dt[, get_SE_identifiers(se, "drug_moa")] <- untreated

  data.table::as.data.table(dt)
}


#' Flatten a table
#'
#' Flatten a stacked table into a wide format.
#'
#' @param tbl a table to flatten.
#' @param groups character vector of column names representing uniquifying groups in expansion.
#' @param wide_cols character vector of column names to flatten.
#' @param sep string representing separator between \code{wide_cols} columns, used in column renaming.
#' Defaults to \code{"_"}.
#'
#' @return table of flattened data as defined by \code{wide_cols}.
#'
#' @details flattened columns will be named with original column names prefixed by \code{wide_cols} columns,
#' concatenated together and separated by \code{sep}.
#'
#' A common use case for this function is when a flattened version of the \code{"Metrics"} assay is desired.
#'
#' @examples
#'  n <- 4
#'  m <- 5
#'  grid <- expand.grid(normalization_type = c("GR", "RV"),
#'    source = c("GDS", "GDR"))
#'  repgrid <- do.call("rbind", rep(list(grid), m))
#'  repgrid$wide <- seq(m * n)
#'  repgrid$id <- rep(LETTERS[1:m], each = n)
#'
#'  groups <- colnames(grid)
#'  wide_cols <- c("wide")
#'
#'  flatten(repgrid, groups = groups, wide_cols = wide_cols)
#'
#' @export
#'
flatten <- function(tbl, groups, wide_cols, sep = "_") {

  checkmate::assert_character(groups)
  checkmate::assert_character(wide_cols)
  checkmate::assert_string(sep)
  checkmate::assert_true(
    is(tbl, "data.table") ||
      is(tbl, "data.frame") ||
        is(tbl, "DFrame")
  )

  if (!all(groups %in% colnames(tbl))) {
    stop(sprintf("missing expected uniquifying groups: '%s'",
      paste0(setdiff(groups, colnames(tbl)), collapse = ", ")))
  }	
  
  idx <- which(colnames(tbl) %in% groups)
  uniquifying <- subset(tbl, select = idx)
  uniquifying <- unique(uniquifying)

  out <- split(subset(tbl, select = -idx), subset(tbl, select = idx), sep = sep)
  missing <- setdiff(wide_cols, colnames(tbl))
  if (length(missing) != 0L) {
    warning(sprintf("missing listed wide_cols columns: '%s'", paste0(missing, collapse = ", ")))
  }

  rename <- colnames(out[[1]]) %in% wide_cols 
  for (grp in names(out)) {
    group <- out[[grp]]
    colnames(group)[rename] <- paste0(grp, sep, colnames(group)[rename])
    out[[grp]] <- group
  }

  ## Drop empty elements for successful merge.
  filtered <- out[lapply(out, nrow) > 0L]
  Reduce(function(x, y) merge(x, y, by = intersect(names(x), names(y))), filtered)
}


#' Prettify metric names in flat 'Metrics' assay
#'
#' Map existing column names of a flattened 'Metrics' assay to prettified names.
#'
#' @param x a list of names.
#' @param human_readable boolean indicating whether or not to return column names in human readable format.
#' @param normalization_type a character with a specified normalization type.
#' @return character vector of prettified names.
#'
#' @details Rename names that are metrics in a table.
#'
#' A common use case for this function is to convert column names from a flattened version of the \code{"Metrics"} assay.
#' #TODO: Intended to be used both at CL and UI.
#'
#' @export
#'
prettify_flat_metrics <- function(x, human_readable = FALSE, normalization_type = c("GR", "RV")) {

  new_names <- x
  metrics_idx <- c(rep(FALSE, length(x)))

  # convert the metric names into common name for variable
  for (norm in normalization_type) {
    metrics_names <- get_header("metrics_names")[norm, ]
    if (length(metrics_names) == 0L) {
      stop(sprintf("missing normalization type: '%s' from header: 'metrics_names'", norm))
    }

    norm_pattern <- paste0("^", norm, "_")
    
    for (name in metrics_names) {
      name_pattern <- paste0("_", names(metrics_names[metrics_name == name]), "$")

      idx <- grepl(norm_pattern, new_names) & grepl(name_pattern, new_names)
      new_names[idx] <- gsub(name_pattern, paste0("_", name), new_names[idx])
      new_names[idx] <- gsub(norm_pattern, "", new_names[idx])
      metrics_idx[idx] <- TRUE
    }
  }

  # keep track of the non-gDR metrics
  gDR_pattern <- "^gDR_"
  non_gDR_metrics_idx <- metrics_idx & !grepl(gDR_pattern, new_names)

  # gDR is the default name.
  new_names <- gsub(gDR_pattern, "", new_names)

  if (human_readable) {

    # Handle GDS data headers.
    # Move the GDS source info to the end and add '(GDS)'.
    GDS <- "\\(.*?\\)\\(_?GDS_\\)\\(.*?\\)"
    new_names <- gsub(GDS, "\\1 \\3 (\\2)", new_names)

    # rename hard coded metrics and variables
    ## TODO: This looks like it also belongs in the headers.
    display_names <- c("Cell line", 
              "Primary Tissue", 
              "Subtype",
              "Parental cell line",
              "Drug", 
              "Drug MOA", 
              "Nbre of tested conc.", 
              "Highest log10(conc.)",
              "E0", 
              "AOC within set range", 
              "Reference cell division time", 
              "cell division time",
              "GR value", 
              "Relative Viability", 
              "_Mean Viability")
    names(display_names) <- c(get_identifier("cellline_name"), # CellLineName
              get_identifier("cellline_tissue"), # Tissue
              get_identifier("cellline_subtype"), # subtype
              get_identifier("cellline_parental_identifier"), # parental_identifier
              get_identifier("drugname"), # DrugName
              get_identifier("drug_moa"), # drug_moa
              "N_conc", 
              "maxlog10Concentration",
              get_header("metrics_names")["RV","x_0"], # E_0
              "AOC_range", 
              "ReferenceDivisionTime", 
              "DivisionTime",
              "GRvalue",
              "RelativeViability",
              "_mean")

    # TODO: I think this can just be a full vector matching.
    for (i in names(display_names)) {
      new_names <- gsub(paste0("^", i), display_names[i], new_names)
    }

    # replace underscore by space for the remaining metrics
    new_names[metrics_idx] <- gsub("_", " ", new_names[metrics_idx])

    # Replace underscore by space for the Drug/Concentration for co-treatment.
    pattern <- "[0-9]+"
    conc_cotrt <- paste0("^Concentration_", pattern, "$")
    drug_cotrt <- paste0("^", get_identifier("drug"), "_", pattern, "$|^Drug_", pattern, "$")

    replace <- grepl(paste0(conc_cotrt, "|", drug_cotrt), new_names)
    new_names[replace] <- gsub("_", " ", new_names[replace])
  }

  return(new_names)
}
