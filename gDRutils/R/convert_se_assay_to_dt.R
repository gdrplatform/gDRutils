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



#' Pretify metric names in flat 'Metrics' assay
#'
#' Return pretified names for flat 'Metrics' assay
#'
#' @param name_list a list of names.
#' @param human_readable either \code{"variable"} (default) or \code{"gDRviz"}.
#'
#' @return table of data with updated column names.
#'
#' @details Rename names that are  metrics in a table.
#'
#' A common use case for this function is to convert column names from a flattened version of the \code{"Metrics"} assay.
#'
#'
#' @export
#'
prettify_flat_metrics <- function(name_list, human_readable = FALSE) {

  new_names <- name_list
  metrics_idx <- array(FALSE, length(name_list))
  # convert the metric names into common name for variable
  for (normalization_type in c("GR","RV")) {
    metrics_names = get_header("metrics_names")[normalization_type,]
    for (n in metrics_names) {
      idx <- grepl(paste0("^", normalization_type, "_"), new_names) & 
        grepl(paste0("_",names(metrics_names[metrics_names==n]), "$"), new_names)
      new_names[idx] <- gsub(paste0("_",names(metrics_names[metrics_names==n]), "$"), 
                              paste0("_",n), new_names[idx])
      new_names[idx] <- gsub(paste0("^", normalization_type, "_"), 
                              "", new_names[idx])
      metrics_idx[idx] <- TRUE
    }
  }
  # keep track of the non-gDR metrics
  non_gDR_metrics_idx <- metrics_idx & !grepl("^gDR_", new_names)
  # scratch gDR as this is the default name
  new_names <- gsub("^gDR_", "", new_names)

  

  # convert into user friendly string for visualization
  if (gDRviz_format) {

    source_type_prefix <- sapply(strsplit(new_names, "_"), function(x) x[[1]][1])
    for (i in which(non_gDR_metrics_idx)) {
      # move the fit source at the end and add ( )
      new_names[i] <- paste0(gsub(paste0("^", source_type_prefix[i], "_"), "", new_names[i]),
                      " (", source_type_prefix[i], ")")
    }

    # rename hard coded metrics and variables
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
    for (i in names(display_names)) {
      new_names <- gsub(paste0("^", i), display_names[i], new_names)
    }

    # replace underscore by space for the remaining metrics
    new_names[metrics_idx] = gsub("_", " ", new_names[metrics_idx])

    # replace underscore by space for the Drug/Concentration for co-treatment
    for (i in 2:20) {
      idx <- grepl(paste0("^Concentration_", i, "$"), new_names)
      new_names[idx] <- gsub("_", " ", new_names[idx])
      idx <- grepl(paste0("^", get_identifier("drug"), "_", i, "$"), new_names)
      new_names[idx] <- gsub("_", " ", new_names[idx])
      idx <- grepl(paste0("^Drug_", i, "$"), new_names)
      new_names[idx] <- gsub("_", " ", new_names[idx])
    }
  }
  return(new_names)
}