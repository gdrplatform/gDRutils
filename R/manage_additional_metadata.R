#' add arbitrary S3 class to an object
#'
#' Modify and object's \code{class} attribute.
#'
#' This is a simple convenience function that an item to the \code{class} attribute of an object
#' so that it can be dispatched to a proper S3 method. This is purely for code clarity,
#' so that individual methods do not clutter the definitions of higher order functions.
#'
#' @param x an object
#' @param newClass character string; class to be added
#' @keywords metadata_management
#'
#' @return The same object with an added S3 class.
#' 
#' @examples 
#' addClass(data.table::data.table(), "someClass")
#'
#' @export
#'
addClass <- function(x, newClass) {
  checkmate::assert_string(newClass)
  if (!is(x, newClass)) {
    class(x) <- c(newClass, class(x))
  }
  return(x)
}

#' modify assay with additional data
#'
#' Reduce dimensionality of an assay by dropping extra data or combining variables.
#'
#' If an essay extracted from a \code{SummarizedExperiment} contains additional information,
#' i.e. factors beyond \code{DrugName} and \code{CellLineName}, that information will be treated
#' in one of three ways, depending on the value of \code{option}:
#'
#' \itemize{
#'   \item{\code{drop}: Some information will be discarded and only one value
#'                      of the additional variable (chosen by the user) will be kept.
#'   }
#'   \item{\code{toDrug}: The information is pasted together with the primary drug name.
#'                        All observations are kept.
#'   }
#'   \item{\code{toCellLine}: As above, but pasting is done with cell line name.
#'   }

#' }
#'
#' Depending on the type of additional information, the exact details will differ.
#' This is handled in the app by adding special classes to the data tables and dispatching to S3 methods.
#'
#' Following modification, the additional columns are discarded.
#'
#' @param x a \code{data.table} containing an assay
#' @param ... additional arguments passed to methods
#' @param option character string specifying the action to be taken, see \code{Details}
#' @param keep character string specifying the value of the active variable that will be kept
#' @keywords metadata_management
#'
#' @return modified object
#'
#' @examples
#' dt <- data.table::data.table(a = as.character(1:10), b = "data")
#' dt <- addClass(dt, "a")
#' modifyData(dt, "average", keep = "b")
#'
#' @export
#'
modifyData <- function(x, ...) {
  UseMethod("modifyData")
}

#' @keywords metadata_management
#' @export
#' @describeIn modifyData includes the name and concentration of the second drug
modifyData.drug_name2 <- function(x, option, keep, ...) {
  checkmate::assert_data_table(x)
  checkmate::assert_string(option)
  checkmate::assert_choice(option, c("average", "toDrug", "toCellLine"))
  checkmate::assert_string(keep, null.ok = TRUE)
  
  pidfs <- get_prettified_identifiers(simplify = TRUE)
  drug_name <- pidfs[["drug_name"]]
  drug_name2 <- pidfs[["drug_name2"]]
  conc2 <- pidfs[["concentration2"]]
  drug2 <- pidfs[["drug2"]]
  cell_name <- pidfs[["cellline_name"]]
  
  if (option == "average") {
    # drop data and keep only the requested value
    x <- average_biological_replicates_dt(x, drug_name2, prettified = TRUE)
  } else {
    # ensure concentration of co-drug is a numeric value
    if (is.factor(x[[conc2]])) {
      x[[conc2]] <- as.character(x[[conc2]])
    }
    if (is.character(x[[conc2]])) {
      x[[conc2]] <- as.numeric(x[[conc2]])
    }
    if (option == "toDrug") {
      x[[drug_name]] <-
        sprintf("%s (%s = %s at %.4f &mu;M)", x[[drug_name]], drug_name2, x[[drug_name2]], x[[conc2]])
      x[[drug_name]] <-
        sub(" \\(.*? at 0\\.?0* &mu;M\\)", "", x[[drug_name]])
    } else if (option == "toCellLine") {
      x[[cell_name]] <-
        sprintf("%s (%s = %s at %.4f &mu;M)", x[[cell_name]], drug_name2, x[[drug_name2]], x[[conc2]])
      x[[cell_name]] <-
        sub(" \\(.*? at 0\\.?0* &mu;M\\)", "", x[[cell_name]])
    }
  }
  
  # drop the additional columns
  x[c(drug_name2, conc2, drug2)] <- NULL
  # remove special class
  class(x) <- setdiff(class(x), "drug_name2")
  return(x)
}

#' @keywords metadata_management
#' @export
#' @describeIn modifyData includes the data source
modifyData.data_source <- function(x, option, keep, ...) {
  checkmate::assert_data_table(x)
  checkmate::assert_string(option)
  checkmate::assert_choice(option, c("average", "toDrug", "toCellLine"))
  checkmate::assert_string(keep, null.ok = TRUE)
  
  pidfs <- get_prettified_identifiers(simplify = TRUE)
  dt_src <- pidfs[["data_source"]]
  drug <- pidfs[["drug_name"]]
  clid <- pidfs[["cellline"]]
  cl_name <- pidfs[["cellline_name"]]
  
  if (option == "average") {
    # drop data and keep only the requested value
    x <- average_biological_replicates_dt(x, dt_src, prettified = TRUE)
  } else {
    duplicated_rows <- get_duplicated_rows(x, c(drug, clid))
    if (option == "toDrug") {
      drugs_to_combine <- unique(x[duplicated_rows, drug])
      drug_idx <- which(x[[drug]] %in% drugs_to_combine)
      drug_to_replace <- x[drug_idx, drug]
      x[drug_idx, drug] <-
        vapply(
          seq_len(length(drug_to_replace)),
          function(y) sprintf("%s (%s)", drug_to_replace[y], x[, dt_src][y]), "string")
    } else if (option == "toCellLine") {
      cell_lines_to_combine <- unique(x[duplicated_rows, cl_name])
      cell_line_idx <- which(x[[cl_name]] %in% cell_lines_to_combine)
      x[cell_line_idx, cl_name] <-
        sprintf("%s (%s)", x[cell_line_idx, cl_name], x[cell_line_idx, dt_src])
    }
  }
  # drop the additional columns
  x[, c(dt_src) := NULL]
  # remove special class
  class(x) <- setdiff(class(x), dt_src)
  return(x)
}


#' @keywords metadata_management
#' @export
#' @describeIn modifyData includes the name of other additional variables
modifyData.default <- function(x, option, keep, ...) {
  checkmate::assert_data_table(x)
  checkmate::assert_string(option)
  checkmate::assert_choice(option, c("average", "toDrug", "toCellLine"))
  checkmate::assert_string(keep, null.ok = TRUE)
  pidfs <- get_prettified_identifiers(simplify = TRUE)
  additional_var_names <- class(x)[[1]]
  additional_var <- ifelse(additional_var_names %in% names(pidfs),
                           pidfs[[additional_var_names]],
                           additional_var_names)
  cell_name <- pidfs[["cellline_name"]]
  drug_name <- pidfs[["drug_name"]]
  
  if (option == "average") {
    # drop data and keep only the requested value
    x <- average_biological_replicates_dt(x, additional_var, prettified = TRUE)
  } else {
    if (option == "toDrug") {
      x <- modify_label(x, drug_name, additional_var)
    } else if (option == "toCellLine") {
      x <- modify_label(x, cell_name, additional_var)
    }
  }
  # drop the additional columns
  if ("data.table" %in% class(x)) {
    x[[additional_var]] <- NULL
  } else {
    x[additional_var] <- NULL
  }
  # remove special class
  class(x) <- setdiff(class(x), additional_var_names)
  return(x)
}

#' @keywords internal
modify_label <- function(x, option, var_name) {
  x[[option]] <-
    sprintf("%s (%s = %s)", x[[option]], var_name, as.character(x[[var_name]]))
  x[[option]] <-
    sub(" \\(.*? at 0\\.?0* &mu;M\\)", "", x[[option]])
  x
}
