ASSAY_INFO_TBL <- data.table::data.table(
  name = c("rawTreated", "Controls", "Normalized", "Averaged", "Metrics",
           "excess", "scores", "isobolograms"),
  pname = c("Raw Treated", "Controls", "Normalized", "Averaged", "Metrics",
            "Excess", "Scores", "isobolograms"),
  type = c("raw", "control", "normalized", "averaged", "metrics", "excess", "scores", "isobolograms"),
  group = c("raw", "core", "core", "core", "core", "combo_excess", "combo_score", "combo_iso"),
  data_type = c("single-agent", "single-agent", "single-agent", "single-agent",
                "single-agent", "combo", "combo", "combo")
)


#' get default assay names
#' for the specified filters,
#' i.e. set of assay types, assay groups and assay data types
#'
#' @param type charvec of assay types
#' @param group charvec of assay groups
#' @param data_type charvec assay of data types
#' @param prettify logical flag, prettify the assay name?
#' @param simplify logical flag, simplify the output?
#'    will return single string instead of named vector with single element
#'    useful when function is expected to return single element/assay only
#'    
#' @return charvec
#' 
#' @examples 
#' get_env_assay_names()
#'
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @export
get_env_assay_names <-
  function(type = NULL,
           group = NULL,
           data_type = NULL,
           prettify = FALSE,
           simplify = TRUE) {
    assert_choices(type, choices = ASSAY_INFO_TBL$type, null.ok = TRUE)
    assert_choices(group, choices = ASSAY_INFO_TBL$group, null.ok = TRUE)
    assert_choices(data_type,
                   choices = ASSAY_INFO_TBL$data_type,
                   null.ok = TRUE)
    checkmate::assert_flag(prettify)

    fname <- if (prettify) {
      "pname"
    } else {
      "name"
    }

    filters <-
      list(type = type,
           group = group,
           data_type = data_type)
    are_null <-
      vapply(filters, is.null, logical(1))
    v_filters <- filters[!are_null]

    if (all(are_null)) {
      return(structure(ASSAY_INFO_TBL[[fname]], names = ASSAY_INFO_TBL$type))
    }
    df <- ASSAY_INFO_TBL
    for (filter in names(v_filters)) {
      df <- df[df[[eval(filter)]] %in% v_filters[[filter]], ]
    }

    if ((nrow(df)) == 0)  {
      v_filters_str <-
        paste(names(v_filters),
              v_filters,
              sep = "=",
              collapse = "&")
      stop(sprintf("Assay name not found for '%s'", v_filters_str))
    }
    res <- structure(df[[fname]], names = df$type)
    if (length(res) == 1L && simplify) {
      res <- as.character(res)
    }
    res
  }

#' get assay names of the given se/dataset
#' fetch the data from the se if provided as metadata
#' use predefined values from `get_env_assay_names` otherwise
#'
#' @param se SummarizedExperiment or NULL
#' @param ... Additional arguments to pass to \code{get_env_assay_names}.
#'
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @return  charvec
#'
#' @examples 
#' get_assay_names()
#' 
#' @export
get_assay_names <- function(se = NULL, ...) {
  if (!is.null(se) &&
      !is.null(S4Vectors::metadata(se = NULL)[["assay_info_tbl"]])) {
    # TODO: extend the logic to support metadata from se
    # i.e. sth like `get_se_assay_names(group = "combo_base")`
    # more details: https://jira.gene.com/jira/browse/GDR-1116
    stop("'get_assay_name' currently does not support non-null 'se'")
  } else {
    get_env_assay_names(...)
  }
}

#' get names of combo assays
#'
#' @param se SummarizedExperiment or NULL
#' @param ... Additional arguments to pass to \code{get_assay_names}.
#' @return charvec of combo assay names.
#' @export
#' @examples 
#' get_combo_assay_names()
#' 
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
get_combo_assay_names <- function(se = NULL, ...) {
  get_assay_names(se, data_type = "combo", simplify = FALSE, ...)
}

#' get names of combo base assays
#'
#' @param se SummarizedExperiment or NULL
#' @param ... Additional arguments to pass to \code{get_combo_assay_names}.
#'
#' @return  charvec
#' @export
#'
#' @examples 
#' get_combo_base_assay_names()
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
get_combo_base_assay_names <- function(se = NULL, ...) {
  get_combo_assay_names(group = "combo_excess", ...)
}

#' get names of combo score assays
#'
#' @param se SummarizedExperiment or NULL
#' @param ... Additional arguments to pass to \code{get_combo_assay_names}.
#'
#' @return  charvec
#'
#' @export
#' 
#' @examples 
#' get_combo_score_assay_names()
#' 
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
get_combo_score_assay_names <- function(se = NULL, ...) {
  get_combo_assay_names(group = "combo_score", ...)
}
