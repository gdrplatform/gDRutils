#' Validate JSON against a schema.
#'
#' Validate JSON describing an object against a schema.
#'
#' @param json String of JSON in memory.
#' @param schema_path String of the schema to validate against.
#'
#' @return Boolean of whether or not JSON successfully validated.
#'
#' @examples
#' json <- '{}'
#'
#' @details This is most often used to validate JSON
#' before passing it in as a document to an ElasticSearch index.
#' @export
validate_json <- function(json, schema_path) {
  # asserts that should be valid regardless user data
  stopifnot(file.exists(schema_path))
  
  # status list returned by function
  
  # errors/issues dependent on the user data
  # let's start with the basic JSON validation
  vjson <- jsonlite::validate(json)
  stl <- if (!isTRUE(vjson)) {
    derror <-
      data.frame(field = "-", message = attributes(vjson)[["err"]])
    list(error = "global JSON validation failed",
         exit_code = 1,
         derror = derror)
  } else {
    v <- jsonvalidate::json_validator(schema_path, engine = "ajv")
    dvjson <- v(json,
                verbose = TRUE,
                greedy = TRUE,
                error = FALSE)
    
    if (!isTRUE(dvjson)) {
      list(
        error = "JSON validation of data model failed",
        exit_code = 2,
        derror = attributes(dvjson)[["errors"]]
      )
    }
    else {
      list(error = NULL,
           exit_code = 0,
           derror = NULL)
    }
  }
 
  st <- stl$exit_code == 0

  attributes(st) <- stl
  st
}

#' Validate MAE against a schema.
#'
#' Validate MAE object against a schema.
#' Currently only SEs are validated
#' TODO: add mae.json schema and validate full MAE object
#'
#' @param mae MultiAssayExperiment object
#' @param schema_package string name of the package with JSON schema files
#' @param schema_dir_path path to the dir with JSON schema files
#' @param schema named charvec with filenames of schemas to validate against.
#'
#' @return Boolean of whether or not mae is valid
#'
#' @export
validate_mae_with_schema <-
  function(mae,
           schema_package = Sys.getenv("SCHEMA_PACKAGE", "gDRutils"),
           schema_dir_path = Sys.getenv("SCHEMA_DIR_PATH", "schemas"),
           schema = c(se = "se.json", mae = "mae.json")) {
    
    checkmate::assert_class(mae, "MultiAssayExperiment")
    checkmate::assert_string(schema_package)
    checkmate::assert_character(schema_dir_path)
    checkmate::assert_character(schema, names = "unique")
    
    schema_dir_apath <- system.file(package = schema_package, schema_dir_path)
    checkmate::assert_directory_exists(schema_dir_apath)
    
    experiments <- names(mae)
     
    # firstly convert mae to json
    ljson <- convert_mae_to_json(mae)
    
    # validate on the mae level
   
    mae_schema_path <- file.path(schema_dir_apath, schema[["mae"]])
    mtl <-
      list(mae = validate_json(ljson[["mae"]], mae_schema_path))
    
    # validate se experiments
    se_schema_path <- file.path(schema_dir_apath, schema[["se"]])
    stl <- lapply(experiments, function(x) {
      validate_json(ljson[["se"]][[x]], se_schema_path)
    })
    names(stl) <- paste0("experiment:", experiments)
    
    st <- if (isTRUE(mtl) && all(vapply(stl, isTRUE, logical(1)))) {
      TRUE
    } else {
      c(mtl, stl)
    }
    st
  }
