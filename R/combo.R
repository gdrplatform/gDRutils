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


#' get_iso_colors
#' 
#' 
#' @param  normalization_type charvec normalization_types expected in the data
#' @keywords combination_data
#'
#' @return named charvec with iso colors
#' 
#' @examples 
#' get_iso_colors()
#' 
#' @export
get_iso_colors <-
  function(normalization_type = c("RV", "GR")) {
    normalization_type <- match.arg(normalization_type)
    iso_cutoff <- seq(0, 1, 0.05)
    breaks <- iso_cutoff
    if (normalization_type == "GR") {
      colors <- vapply(iso_cutoff, function(x) {
        color_vector <- c(70, round((0.85 - x * 0.7) * 170), round((1.1 - x * 0.7) * 200))
        assert_RGB_format(color_vector)
        sprintf("#%s", paste(as.hexmode(color_vector), collapse = ""))
      },
      character(1)
      )
    } else {
      colors <- vapply(iso_cutoff, function(x) {
        color_vector <- c(70, round((1 - x * .85) * 170), round((1.1 - x * .85) * 232))
        assert_RGB_format(color_vector)
        sprintf("#%s", paste(as.hexmode(color_vector), collapse = ""))
      },
      character(1)
      )
    }
    names(colors) <- iso_cutoff
    colors
  }

assert_RGB_format <- function(x) {
  if (any(x > 255)) {
    stop("Some value is greater than 255. Not valid RGB format.")
  }
}

#' Get colorscale data for given combo assay and growth metric
#'
#' @param g_metric growth metric
#' @param assay_type assay type
#' @keywords combination_data
#'
#' @return list with colors, breaks and limits
#' @examples 
#' get_combo_col_settings("GR", "excess")
#' 
#' @export
get_combo_col_settings <-
  function(g_metric,
           assay_type) {
    assay_names <- get_combo_assay_names()
    assay_types <- names(assay_names)
    checkmate::assert_choice(g_metric, c("RV", "GR"))
    checkmate::assert_choice(assay_type, assay_types)

    colors <- breaks <- limits <- NULL
    if (assay_type %in% names(get_combo_assay_names(group = "combo_iso"))) {
      myv <- get_iso_colors(g_metric)
      colors <- as.character(myv)
      breaks <- names(myv)
    } else if (assay_type %in% c(names(get_combo_assay_names(group = "combo_excess")),
                                 names(get_combo_assay_names(group = "combo_score")))) {
      colors <- c("#003355", "#4488dd", "#eeeedd", "#CC8844", "#662200")
      if (g_metric == "GR") {
        breaks <- c(-0.5, -0.25, 0, 0.25, 0.5)
      } else if (g_metric == "RV") {
        breaks <- c(-0.3, -0.15, 0, 0.15, 0.3)
      } else {
        stop(sprintf("unexpected 'g_metric' type: '%s'", g_metric))
      }
    } else if (assay_type %in% names(get_combo_assay_names(group = "combo_score_mx"))) {
      colors <- c("#003355", "#4488dd", "#eeeedd", "#CC8844", "#662200")
      breaks <- c(-4, -2, 0, 2, 4)
    } else if (assay_type %in% names(get_combo_assay_names(group = "combo_base_mx"))) {
      if (g_metric == "GR") {
        colors <- c("#001155", "#1122AA", "#AA4400", "#FF7711", "#ffffee")
        breaks <- c(-0.6, -0.2, 0.2, 0.6, 1)
      } else if (g_metric == "RV") {
        colors <- c("#330033", "#770033", "#BB6633", "#ffaa11", "#ffffee")
        breaks <- c(0, 0.25, 0.5, 0.75, 1)
      } else {
        stop(sprintf("unexpected 'g_metric' type: '%s'", g_metric))
      }
    } else {
      stop(sprintf("no logic found for the assay_type: '%s'", assay_type))
    }
    stopifnot(
      "unexpected error when determining combo color settings - either 'colors' or 'breaks' is NULL" = 
        !(is.null(colors) || is.null(breaks))
    )
    list(
      colors = colors,
      breaks = breaks,
      limits = c(min(breaks), max(breaks))
    )
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
#' @export
define_matrix_grid_positions <- function(conc1, conc2) {
  .generate_gap_for_single_agent <- function(x) {
    2 * x[2] - x[3] - log10(1.5)
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
#' @keywords utils
#' @export
round_concentration <- function(x, ndigit = 3) {
  round(10 ^ (round(log10(x), ndigit)), ndigit - 1 - floor(log10(x)))
}

#' @keywords internal
#' @noRd
rbindParallelList <- function(x, name) {
  S4Vectors::DataFrame(
    do.call(
      rbind, 
      c(lapply(x, function(x) {
        dt <- data.table::as.data.table("[[" (x, name))
        data.table::setorder(dt)
        dt
      }), fill = TRUE)
    )
  )
}
