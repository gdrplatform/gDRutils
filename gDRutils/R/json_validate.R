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
#' @importFrom jsonvalidate json_validate
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
    dvjson <-
      jsonvalidate::json_validate(
        json,
        schema_path,
        verbose = TRUE,
        greedy = TRUE,
        error = FALSE
      )
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
#' @param schema_path path to the dir with JSON schema files
#' @param schema named charvec with filenames of schemas to validate against.
#'
#' @return Boolean of whether or not mae is valid
#'
#' @importFrom jsonvalidate json_validator
#' @export
validate_mae_with_schema <-
  function(mae,
           schema_path = system.file(package = "gDRutils", "schemas"),
           schema = c(se = "se.json", mae = "mae.json")) {
    checkmate::assert_class(mae, "MultiAssayExperiment")
    experiments <- names(mae)
    checkmate::assert_subset(experiments,
                             c(
                               "single-agent",
                               "co-dilution",
                               "matrix",
                               "cotreatment",
                               "other"
                             ))
    se_schema_path <- file.path(schema_path, schema[["se"]])
    v_se <- json_validator(se_schema_path)
    stl <- lapply(experiments, function(x) {
      se_json <- convert_se_to_json(mae[[x]])
      validate_json(se_json, se_schema_path)
    })
    names(stl) <- paste0("experiment:", experiments)
    
    st <- if (all(vapply(stl, isTRUE, logical(1)))) {
      TRUE
    } else {
      stl
    }
    st
  }
