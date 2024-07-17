#' convert combo assays from SummarizedExperiments to the list of data.tables
#'
#' @param se \code{SummarizedExperiment} object with dose-response data
#' @param c_assays charvec of combo assays to be used
#' @param normalization_type charvec of normalization_types expected in the data
#' @param prettify boolean flag indicating whether or not to prettify the colnames of the returned data 
#' @keywords combination_data
#' 
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @return list of data.table(s) with combo data
#'
#' @examples 
#' mae <- get_synthetic_data("finalMAE_combo_matrix_small.qs")
#' convert_combo_data_to_dt(mae[[1]])
#'
#' @export
convert_combo_data_to_dt <-
  function(se,
           c_assays = get_combo_assay_names(),
           normalization_type = c("RV", "GR"),
           prettify = TRUE) {

    checkmate::assert_class(se, "SummarizedExperiment")
    checkmate::assert_character(normalization_type)
    assert_choices(c_assays, get_combo_assay_names())
    checkmate::assert_flag(prettify)

    ntype_name <- prettify_flat_metrics("normalization_type", human_readable = FALSE)

    my_l <-
      lapply(names(c_assays), function(x) {
        dt <- convert_se_assay_to_dt(se, c_assays[[x]])
        # TODO: let's discuss how to handle names in isobologram assay (pos_x, iso_level)
        if (prettify) {
          colnames(dt) <-
            prettify_flat_metrics(colnames(dt), human_readable = TRUE)
        }
        dt
      })
    
    # TODO: discuss what should be returned: assay_name or maybe assay_type?
    names(my_l) <- as.character(c_assays)
    my_l
  }


DATA_COMBO_INFO_TBL <- data.table::data.table(
  name = c("hsa_score", "bliss_score", "CIScore_50", "CIScore_80",
           "smooth", "hsa_excess", "bliss_excess"),
  pname = c("HSA Score", "Bliss Score", "log2(CI) @ GR/IC50", "log2(CI) @ GR/IC80",
            "MX full", "HSA excess", "Bliss excess"),
  type = c("scores", "scores", "scores", "scores",
           "excess", "excess", "excess")
)

#' get names of combo score fields
#'
#' @keywords combination_data
#' @return  charvec
#'
#' @export
#' 
#' @examples 
#' get_combo_score_assay_names()
#' 
get_combo_score_field_names <- function() {
  dt <- DATA_COMBO_INFO_TBL[type == "scores", c("name", "pname"), with = FALSE]
  stats::setNames(dt$pname, dt$name)
}

#' get names of combo excess fields
#'
#' @keywords combination_data
#' @return charvec
#'
#' @export
#' 
#' @examples 
#' get_combo_excess_field_names()
#' 
get_combo_excess_field_names <- function() {
  dt <- DATA_COMBO_INFO_TBL[type == "excess", c("name", "pname"), with = FALSE]
  stats::setNames(dt$pname, dt$name)
}


#' get combo assay names based on the field name
#'
#'
#' @param field String containing name of the field for which the assay name should be returned
#' @keywords combination_data
#' @return charvec
#'
#' @export
#' 
#' @examples 
#' convert_combo_field_to_assay("hsa_score")
#' 
convert_combo_field_to_assay <- function(field) {
  checkmate::assert_string(field)
  DATA_COMBO_INFO_TBL[name == field, ][["type"]]
}


#' Define matrix grid positions
#'
#' @param conc1 drug_1 concentration
#' @param conc2 drug_2 concentration
#'
#' @details
#' \code{drug_1} is diluted along the rows as the y-axis and
#' \code{drug_2} is diluted along the columns and will be the x-axis.
#' 
#' @keywords combination_data
#' @return list with axis grid positions
#' 
#' @examples
#' cl_name <- "cellline_BC"
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' se <- get_synthetic_data("combo_matrix_small")[["combination"]]
#' dt_average <- convert_se_assay_to_dt(se, "Averaged")[normalization_type == "GR"]
#' 
#' ls_axes <- define_matrix_grid_positions(
#'    dt_average[["Concentration"]], dt_average[["Concentration_2"]])
#' 
#' @export
define_matrix_grid_positions <- function(conc1, conc2) {
  checkmate::assert_numeric(conc1)
  checkmate::assert_numeric(conc2)
  
  .generate_gap_for_single_agent <- function(x) {
    if (NROW(x) == 1) {
      x
    } else if (NROW(x) > 2) {
    2 * x[2] - x[3] - log10(1.5)
    } else {
      x[2] - 0.5
    }
  } 
  
  conc_1 <- sort(unique(round_concentration(conc1)))
  pos_y <- log10conc_1 <- log10(conc_1)
  pos_y[1] <- .generate_gap_for_single_agent(log10conc_1)
  axis_1 <- data.table::data.table(conc_1 = conc_1,
                                   log10conc_1 = log10conc_1,
                                   pos_y = pos_y,
                                   marks_y = sprintf("%.2g", conc_1)
  )
  
  conc_2 <- sort(unique(round_concentration(conc2)))
  pos_x <- log10conc_2 <- log10(conc_2)
  pos_x[1] <- .generate_gap_for_single_agent(log10conc_2)
  axis_2 <- data.table::data.table(conc_2 = conc_2,
                                   log10conc_2 = log10conc_2,
                                   pos_x = pos_x,
                                   marks_x = sprintf("%.2g", conc_2)
  )
  
  list(axis_1 = axis_1, axis_2 = axis_2)
}

#' Round concentration to ndigit significant digits
#'
#' @param x value to be rounded.
#' @param ndigit number of significant digits (default = 4).
#' 
#' @examples 
#' round_concentration(x = c(0.00175,0.00324,0.0091), ndigit = 1)
#'
#' @return rounded x
#' @keywords combination_data
#' @export
round_concentration <- function(x, ndigit = 3) {
  checkmate::assert_numeric(x)
  checkmate::assert_integerish(ndigit)
  
  round(10 ^ (round(log10(x), ndigit)), ndigit - 1 - floor(log10(x)))
}
