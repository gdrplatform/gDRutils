#' Create JSON document.
#'
#' Convert a MultiAssayExperiment object to a JSON document.
#'
#' @param mae SummarizedExperiment object.
#' @param with_experiments logical convert experiment metadata as well?
#'
#' @return String representation of a JSON document.
#'
#' @export
convert_mae_to_json <- function(mae, with_experiments = TRUE) {

  ljson <- list()

  ljson[["mae"]] <- .convert_mae_summary_to_json(mae)

  if (with_experiments) {

    ljson[["se"]] <- MAEpply(mae, convert_se_to_json)
  }

  ljson
}

#' Create JSON document with MAE summary
#'
#' Create JSON document with MAE summary,
#'  currently only experiment names
#'
#' @param mae MultiAssayExperiment object.
#'
#' @return String representation of a JSON document.
#'
.convert_mae_summary_to_json <- function(mae) {
  ml <- list()
  ml[["experiment_names"]] <-
    names(mae)

  # currently there is no colData and metadata on MAE level
  # if data model changes one has to add validation of these properties below

  jsonlite::toJSON(ml)
}


#' Create JSON document.
#'
#' Convert a SummarizedExperiment object to a JSON document.
#'
#' @param se SummarizedExperiment object.
#'
#' @return String representation of a JSON document.
#'
#' @examples
#' md <- list(title = "my awesome experiment",
#'   description = "description of experiment",
#'   source = list(name = "GeneData_Screener", id = "QCS-12345"))
#' rdata <- data.table::setDT(mydrug = letters, mydrugname = letters, mydrugmoa = letters, Duration = 1)
#' cdata <- data.table::setDT(mycellline = letters, mycelllinename = letters,
#'  mycelllinetissue = letters, cellline_ref_div_time = letters)
#' identifiers <- list(cellline = "mycellline",
#'                     cellline_name = "mycelllinename",
#'                     cellline_tissue = "mycelllinetissue",
#'                     cellline_ref_div_time = "cellline_ref_div_time",
#'                     drug = "mydrug",
#'                     drug_name = "mydrugname",
#'                     drug_moa = "mydrugmoa",
#'                     duration = "Duration")
#' se <- SummarizedExperiment::SummarizedExperiment(rowData = rdata,
#'                                                  colData = cdata)
#' se <- gDRutils::set_SE_experiment_metadata(se, md)
#' se <- gDRutils::set_SE_identifiers(se, identifiers)
#' convert_se_to_json(se)
#'
#' @export
convert_se_to_json <- function(se) {
  ml <- list(
    mjson = .convert_metadata_to_json(se),
    rjson = .convert_rowData_to_json(rowData(se), get_SE_identifiers(se)),
    cjson = .convert_colData_to_json(colData(se), get_SE_identifiers(se))
  )

  # filter out empty strings
  # otherwise JSON data will be invalid
  # and we want be able to use it with jsonvalidate
  ml <- ml[vapply(ml, nchar, integer(1)) > 0]

  json <- sprintf("{%s}", toString(ml))
  stopifnot(jsonlite::validate(json))
  # minify is used here to get rid of all the backslashing that
  # happens when using paste0/paste/sprintf on json objects and or strings
  json <- jsonlite::minify(json)
  json
}


#' Convert experiment metadata to JSON
#'
#' Convert experiment metadata to JSON format for elasticsearch indexing.
#'
#' @param se SummarizedExperiment object.
#'
#' @return JSON string capturing experiment metadata.
#'
#' @examples
#' md <- list(title = "my awesome experiment",
#'   description = "description of experiment",
#'   sources = list(list(name = "GeneData_Screener", id = "QCS-12345")))
#' se <- SummarizedExperiment::SummarizedExperiment(metadata = md)
#' gDRutils:::.convert_metadata_to_json(se)
#'
#' @keywords internal
.convert_metadata_to_json <- function(se) {

  md <- get_SE_experiment_metadata(se)

  if (inherits(md, "data.table") || inherits(md, "DataFrame")) {
    # data.table after converting to JSON have strange format, which causing error in `jsonlite::validate`
    md <- as.list(md)
  }

  mjson <- jsonlite::toJSON(md, auto_unbox = TRUE)
  strip_first_and_last_char(mjson)
}


#' Convert rowData to JSON
#'
#' Convert rowData to JSON format for elasticsearch indexing.
#'
#' @param rdata data.table of \code{rowData}.
#' @param identifiers charvec with identifiers
#' @param req_cols charvec required columns
#'
#' @return JSON string capturing the \code{rdata}.
#'
#' @examples
#' rdata <- data.table::setDT(mydrug = letters, mydrugname = letters, mydrugmoa = letters, Duration = 1)
#' identifiers <- list(drug = "mydrug", drug_name = "mydrugname", drug_moa = "mydrugmoa",
#' duration = "Duration")
#' gDRutils:::.convert_rowData_to_json(rdata, identifiers)
#'
#' @details Standardizes the \code{rdata} to common schema fields
#' and tidies formatting to be condusive to joining
#' with other JSON responses.
#' @keywords internal
.convert_rowData_to_json <-
  function(rdata,
           identifiers,
           req_cols = c("drug", "drug_name", "drug_moa", "duration")) {

  json <- .standardize_and_convert_element_metadata_to_json(rdata, identifiers, req_cols)
  sprintf("%s, \"misc_rowdata\": %s", json$main, json$opt)
}


#' Convert colData to JSON
#'
#' Convert colData to JSON format for elasticsearch indexing.
#'
#' @param cdata data.table of \code{colData}.
#' @param identifiers charvec with identifiers
#' @param req_cols charvec required columns
#'
#' @return JSON string capturing the \code{cdata}.
#'
#' @examples
#' cdata <- data.table::setDT(mycellline = letters, mycelllinename = letters, mycelllinetissue = letters,
#'                    cellline_ref_div_time = "cellline_ref_div_time")
#' identifiers <- list(cellline = "mycellline",
#'                     cellline_name = "mycelllinename",
#'                     cellline_ref_div_time = "cellline_ref_div_time",
#'                     cellline_tissue = "mycelllinetissue")
#' gDRutils:::.convert_colData_to_json(cdata, identifiers)
#'
#' @details Standardizes the \code{cdata} to common schema fields
#' and tidies formatting to be condusive to joining
#' with other JSON responses.
#' @keywords internal
.convert_colData_to_json <-
  function(cdata,
           identifiers,
           req_cols = c("cellline", "cellline_name", "cellline_tissue", "cellline_ref_div_time")) {

  json <- .standardize_and_convert_element_metadata_to_json(cdata, identifiers, req_cols)
  sprintf("%s, \"misc_coldata\": %s", json$main, json$opt)
}


#' @keywords internal
.standardize_and_convert_element_metadata_to_json <- function(data, identifiers, req_cols) {
  data <- .standardize_column_names(data, identifiers[req_cols])
  .convert_element_metadata_to_json(data, req_cols)
}


#' @keywords internal
.standardize_column_names <- function(mdata, identifiers) {
  stopifnot(all(identifiers %in% names(mdata)))
  colnames(mdata)[match(identifiers, colnames(mdata))] <- names(identifiers)
  mdata
}


#' @keywords internal
.convert_element_metadata_to_json <- function(mdata, req_cols) {
  stopifnot(all(req_cols %in% names(mdata)))

  mdata <- data.table::setDT(data.frame(data.matrix(mdata)))
  rownames(mdata) <- NULL

  main_mdata <- mdata[, ..req_cols]
  mjson <- jsonlite::toJSON(main_mdata, "columns")
  mjson <- strip_first_and_last_char(mjson)

  opt_cols <- setdiff(colnames(mdata), req_cols)
  opt_mdata <- mdata[, ..opt_cols]
  ojson <- jsonlite::toJSON(opt_mdata, "columns")

  list(main = mjson, opt = ojson)
}


#' String first and last characters of a string.
#'
#' String first and last characters of a string.
#'
#' @param jstring String of any number of characters greater than 1.
#'
#' @return String with first and last characters stripped.
#' @details This is most often used to remove the JSON brackets \code{'{'} and \code{'}'}.
strip_first_and_last_char <- function(jstring) {
  stopifnot(nchar(jstring) > 1)
  substr(jstring, 2, nchar(jstring) - 1)
}
