
#' get settings from JSON file

#' In most common scenario the settings are stored in JSON file 
#' to avoid hardcoding
#'
#' @param s charvec with setting entry/entries
#' @param json_path string with the path to the JSON file
#'
#' @return value/values for entry/entries from JSON file
#' @keywords json_const
#' @export
#'
get_settings_from_json <-
  function(s = NULL,
           json_path = system.file(package = "gDRutils", "cache.json")) {
    
    checkmate::assert_character(s, null.ok = TRUE)
    checkmate::assert_file_exists(json_path)
    
    cache_l <- jsonlite::fromJSON(json_path)
    
    if (!is.null(s)) {
      checkmate::assert_subset(s, names(cache_l))
      cache_l[[s]]
    } else {
      cache_l
    }
  }


#' Get isobologram column names
#'
#' @param k key
#' @param prettify change to upper case and add underscore, iso_level --> Iso_Level
#' 
#' @examples
#' get_isobologram_columns()
#' get_isobologram_columns("iso_level", prettify = TRUE)
#'
#' @keywords json_const
#' @return character vector of isobologram column names for combination data
#' @export
get_isobologram_columns <- function(k = NULL, prettify = TRUE) {
  checkmate::assert_character(k, null.ok = TRUE)
  checkmate::assert_flag(prettify)
  
  ic <- get_settings_from_json("ISOBOLOGRAM_COLUMNS",
                               system.file(package = "gDRutils", "settings.json"))
  if (!is.null(k)) {
    out <- ic[[k]]
  } else {
    out <- ic
  }
  
  gsub(" ", "_", prettify_flat_metrics(out, human_readable = prettify))
}