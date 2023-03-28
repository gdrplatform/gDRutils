#' convert combo assays from SummarizedExperiments to the list of data.tables
#'
#' @param se \code{SummarizedExperiment} object with dose-response data
#' @param c_assays charvec of combo assays to be used
#' @param normalization_type charvec of normalization_types expected in the data
#' @param prettify boolean flag indicating whether or not to prettify the colnames of the returned data 
#' 
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @return list of data.table(s) with combo data
#'
#' @examples 
#' mae <- get_synthetic_data("finalMAE_combo_matrix_small.RDS")
#' convert_combo_data_to_dt(mae[[1]])
#'
#' @export
convert_combo_data_to_dt <-
  function(se,
           c_assays = get_combo_assay_names(),
           normalization_type = c("RelativeViability", "GRvalue"),
           prettify = TRUE) {

    checkmate::assert_class(se, "SummarizedExperiment")
    checkmate::assert_character(normalization_type)
    assert_choices(c_assays, get_combo_assay_names())
    checkmate::assert_flag(prettify)

    ntype_dict <- vapply(normalization_type, shorten_normalization_type_name, character(1))
    ntype_name <- prettify_flat_metrics("normalization_type", human_readable = FALSE)


    my_l <-
      lapply(names(c_assays), function(x) {
        if (x %in% names(get_combo_base_assay_names())) {
          if (x == names(get_combo_assay_names(group = "combo_base_mx"))) {
            dt <- convert_se_assay_to_dt(se, c_assays[[x]])
            for (n in names(ntype_dict)) {
              data.table::setnames(dt, n, ntype_dict[n], skip_absent = TRUE)
            }
            dt <- data.table::melt(
              dt,
              measure.vars = as.character(ntype_dict),
              variable.name = ntype_name,
              value.name = "base-value"
            )
          } else {
            dt <- convert_se_assay_to_dt(se, c_assays[[x]])
            dt[["base-value"]] <- dt[["excess"]]
          }
        } else {
          dt <- convert_se_assay_to_dt(se, c_assays[[x]])
        }
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
#'
#' @return named charvec with iso colors
#' 
#' @examples 
#' get_iso_colors()
#' 
#' @export
get_iso_colors <-
  function(normalization_type = c("RelativeViability", "GRvalue")) {
    normalization_type <- match.arg(normalization_type)
    iso_cutoff <- seq(0, 1, 0.05)
    breaks <- iso_cutoff
    if (normalization_type == "GRvalue") {
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
#'
#' @return list with colors, breaks and limits
#' @examples 
#' get_combo_col_settings("GRvalue", "smooth")
#' 
#' @export
get_combo_col_settings <-
  function(g_metric,
           assay_type) {
    assay_names <- get_combo_assay_names()
    assay_types <- names(assay_names)
    checkmate::assert_choice(g_metric, c("RelativeViability", "GRvalue"))
    checkmate::assert_choice(assay_type, assay_types)

    colors <- breaks <- limits <- NULL
    if (assay_type %in% names(get_combo_assay_names(group = "combo_iso"))) {
      myv <- get_iso_colors(g_metric)
      colors <- as.character(myv)
      breaks <- names(myv)
    } else if (assay_type %in% c(names(get_combo_assay_names(group = "combo_score_excess")),
                                 names(get_combo_assay_names(group = "combo_base_excess")))) {
      colors <- c("#003355", "#4488dd", "#eeeedd", "#CC8844", "#662200")
      if (g_metric == "GRvalue") {
        breaks <- c(-0.5, -0.25, 0, 0.25, 0.5)
      } else if (g_metric == "RelativeViability") {
        breaks <- c(-0.3, -0.15, 0, 0.15, 0.3)
      } else {
        stop(sprintf("unexpected 'g_metric' type: '%s'", g_metric))
      }
    } else if (assay_type %in% names(get_combo_assay_names(group = "combo_score_mx"))) {
      colors <- c("#003355", "#4488dd", "#eeeedd", "#CC8844", "#662200")
      breaks <- c(-4, -2, 0, 2, 4)
    } else if (assay_type %in% names(get_combo_assay_names(group = "combo_base_mx"))) {
      if (g_metric == "GRvalue") {
        colors <- c("#001155", "#1122AA", "#AA4400", "#FF7711", "#ffffee")
        breaks <- c(-0.6, -0.2, 0.2, 0.6, 1)
      } else if (g_metric == "RelativeViability") {
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
