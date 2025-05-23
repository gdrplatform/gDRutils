#' check if data.table contains duplicated data
#' 
#' An auxiliary function that checks for duplicates in the data.table (or its subset)
#' 
#' @param dt data.table
#' @param col_names charvec with columns to be used for subsetting
#' @examples
#' dt <- data.table::data.table(a = c(1, 2, 3), b = c(3, 2, 2))
#' has_dt_duplicated_rows(dt, "b")
#' @return logical flag indicating if a dt contains duplicated rows or not
#' @keywords duplicates
#'
#' @export
#' 
has_dt_duplicated_rows <- function(dt, col_names = NULL) {
  checkmate::assert_data_table(dt)
  checkmate::assert_character(col_names, null.ok = TRUE)
  
  if (is.null(col_names)) {
    anyDuplicated(dt) != 0
  } else {
    checkmate::assert_subset(col_names, colnames(dt))
    anyDuplicated(dt, by = col_names) != 0
  }
  
}

#' get columns in the assay data required to have unique data
#' 
#' get columns in the assay data required to have unique (non-duplicated) data
#' 
#' @param dt data.table with assay data
#' @examples
#' sdata <- get_synthetic_data("finalMAE_small")
#' smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
#' get_assay_req_uniq_cols(smetrics_data)
#' @return charvec with columns required to have unique data
#' @keywords duplicates
#' @export
#' 
get_assay_req_uniq_cols <- function(dt) {
  
  checkmate::assert_data_table(dt)
  col_ids <- get_settings_from_json(
    "assay_dt_req_uniq_col_ids",
    system.file(package = "gDRutils", "settings.json")
  )
 
  # check with both pretiffied and unprettified version of ids
  col_names_p <- unlist(get_prettified_identifiers(col_ids, simplify = FALSE))
  col_names_up <- as.character(get_env_identifiers(col_ids, simplify = FALSE))
  
  # add additional variables
  ## there are columns that should not be considered additional variables 
  skip_col_ids <- get_settings_from_json(
    "assay_dt_skip_uniq_col_ids",
    system.file(package = "gDRutils", "settings.json")
  )
  skip_col_ids <- intersect(colnames(dt), skip_col_ids)
  if (NROW(skip_col_ids)) {
    sdt <- dt[, -skip_col_ids, with = FALSE]
    add_v_p <- gDRutils::get_additional_variables(sdt)
    add_v_up <- gDRutils::get_additional_variables(sdt, prettified = FALSE)
  } else {
    add_v_p <- gDRutils::get_additional_variables(dt)
    add_v_up <- gDRutils::get_additional_variables(dt, prettified = FALSE)
  }
  
  col_names <- unique(c(col_names_p, col_names_up, add_v_p, add_v_up))
  
  intersect(col_names, names(dt))
}

#' check if assay data contains duplicated data
#' 
#' An auxiliary function that checks for duplicates in the assay data
#' 
#' @param dt data.table with assay data
#' 
#' @return logical flag indicating if a dt contains duplicated rows or not
#' @keywords duplicates
#' @examples
#' sdata <- get_synthetic_data("finalMAE_small")
#' smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
#' has_assay_dt_duplicated_rows(smetrics_data)
#' @export
#' 
has_assay_dt_duplicated_rows <- function(dt) { 

  checkmate::assert_data_table(dt)
  
  col_names <- get_assay_req_uniq_cols(dt)
  has_dt_duplicated_rows(dt, col_names)

}


#' Helper function to find duplicated rows
#'
#' @param x DataFrame or data.table
#' @param col_names character vector, columns in which duplication are searched for
#' @param output string with the output format to be returned - 
#' one of "index" (index of duplicates) or "data" (subset of input data with duplicates)
#' @examples
#' dt <- data.table::data.table(a = c(1, 2, 3), b = c(3, 2, 2))
#' get_duplicated_rows(dt, "b")
#' get_duplicated_rows(dt, "b", output = "data")
#' @return integer vector or data.table with duplicated rows
#' @keywords duplicates
#' @export
get_duplicated_rows <- function(x, 
                                col_names = NULL, 
                                output = "index") {
  
  checkmate::assertMultiClass(x, c("data.table", "DataFrame"))
  checkmate::assert_true(all(col_names %in% colnames(x)))
  checkmate::assert_choice(output, c("index", "data"))
  
  
  if (!is.null(col_names)) {
    sub_x <- subset(x, select = col_names)
  } else {
  sub_x <- x
 }
  idx <- which(duplicated(sub_x) | duplicated(sub_x, fromLast = TRUE))
  
  out <- if (output == "index") {
    idx
  } else {
    if (length(idx)) {
      x[idx, ]
    } else {
      x[0, ]
    }
  }
  out
}

#' Helper function to find duplicated rows in assay data
#'
#' @param dt data.table
#' @param output string with the output format to be returned
#' @return integer vector or data.table with duplicated rows
#' @examples
#' sdata <- get_synthetic_data("finalMAE_small")
#' smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
#' get_assay_dt_duplicated_rows(smetrics_data, output = "data")
#' get_assay_dt_duplicated_rows(smetrics_data)
#' @keywords duplicates
#' @export
get_assay_dt_duplicated_rows <- function(dt, output = "index") {
  
  checkmate::assert_data_table(dt)
  
  col_names <- get_assay_req_uniq_cols(dt)
  
  get_duplicated_rows(dt, col_names, output = output)
}


#' throw message if assay data.table contains duplicated rows
#' 
#' An auxiliary function that checks for duplicated rows in assay data.table, 
#' In case of duplicates it throws a message. The messsage function is by default `stop()` 
#' The message function can be customized with `msg_f` parameter
#' 
#' @param dt data.table with assay data
#' @param assay_name string with the name of the assay
#' @param msg_f function to be used to throw the message
#' @param preview_max_numb number of rows to preview if duplicates found
#' @examples
#' sdata <- get_synthetic_data("finalMAE_small")
#' smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
#' throw_msg_if_duplicates(smetrics_data, assay_name = "Metrics", msg_f = futile.logger::flog.info)
#' @return NULL
#' @keywords duplicates
#'
#' @export
#' 
throw_msg_if_duplicates <- function(dt, assay_name = "unknown", msg_f = stop, preview_max_numb = 4) {

  checkmate::assert_data_table(dt)
  checkmate::assert_string(assay_name)
  checkmate::assert_function(msg_f)
  checkmate::assert_number(preview_max_numb)

  if (has_assay_dt_duplicated_rows(dt)) { 

    dup_dt <- get_assay_dt_duplicated_rows(dt, output = "data")
    preview_numb <- min(c(preview_max_numb, NROW(dup_dt)))
   
     msg <- sprintf(
          "The %i ouf of %i rows are duplicated in the assay '%s'",
          NROW(dup_dt),
          NROW(dt),
          assay_name)
     msg2 <- sprintf(" when checking uniquness with the following set of columns: '%s'. ",
          toString(get_assay_req_uniq_cols(dt)))
     msg3 <- sprintf("Here is the preview of the first %i duplicated rows in JSON format: '%s'",
          preview_numb,
          jsonlite::toJSON(dup_dt[seq(preview_numb), ]))
     msg_f(paste0(msg, msg2, msg3))
  }
}

