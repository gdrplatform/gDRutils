#' Transform a SummarizedExperiment assay to a long data.table
#'
#' Transform a SummarizedExperiment assay to a long data.table with a single entry for each row and column combination.
#'
#' @param se \linkS4class{SummarizedExperiment} object with dose-response data.
#' @param assay_name String of name of the assay or index of the assay in the \code{se}.
#' @param merge_metrics Logical indicating whether the metrics should be merged.
#' Defaults to \code{FALSE}.
#' @param include_metadata Boolean indicating whether to include the metadata on the SummarizedExperiment.
#' Defaults to \code{TRUE}.
#'
#' @return data.table with dose-response data
#'
#' @export
#'
assay_to_dt <- function(se,
                        assay_name,
                        merge_metrics = FALSE,
                        include_metadata = TRUE) {

  # Assertions.
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assertTRUE(checkmate::test_count(assay_name) ||
                          checkmate::test_string(assay_name))
  checkmate::assert_flag(merge_metrics)

  .Deprecated(new = convert_se_assay_to_dt,
    msg = "support for 'assay_to_dt' will be dropped next release cycle. See 'convert_se_assay_to_dt' instead")

  if (is.integer(assay_name)) {
    assay_name <- SummarizedExperiment::assayNames(se)[assay_name]
  }

  as_dt <- convert_se_assay_to_dt(se, assay_name, include_metadata = include_metadata)
  if (assay_name == "Metrics") {
    ## TODO: Put in issue to BumpyMatrix::unsplitAsBumpyMatrix to also return nested rownames.
    ## Then can remove all hard-coded logic below regarding metrics.
    as_dt$dr_metric <- rep_len(c("RV", "GR"), nrow(as_dt))

    ## NOTE: assay_to_dt function is deprecated for convert_se_assay_to_dt function,
    ## which noteably does NOT support the merge_metrics argument.
    if (merge_metrics) {
      metric_headers <- get_header("response_metrics")
      if (!all(metric_headers %in% colnames(as_dt))) {
        stop(sprintf("missing expected metric headers: '%s'",
          paste0(setdiff(metric_headers, colnames(as_dt)), collapse = ", ")))
      }

      id_headers <- c("rId", "cId")
      headers <- c(id_headers, metric_headers)

      metric_headers <- colnames(as_dt)[colnames(as_dt) %in% metric_headers]

      Df_RV <- as_dt[dr_metric == "RV", ..headers]
      rv_map <- get_header("RV_metrics")

      Df_GR <- as_dt[dr_metric == "GR", - "dr_metric"]
      gr_map <- get_header("GR_metrics")

      data.table::setnames(Df_RV,
               old = metric_headers,
               new = unname(rv_map[metric_headers]))

      data.table::setnames(Df_GR,
               old = metric_headers,
               new = unname(gr_map[metric_headers]))

      as_dt <- merge(Df_RV,
                     Df_GR,
                     by = id_headers,
                     all = TRUE)
    }
  }
  return(as_dt)
}


#' @export
#'
RVGRfits <- function(df_,
                     e_0 = 1,
                     GR_0 = 1,
                     n_point_cutoff = 4,
                     range_conc = c(5e-3, 5),
                     force = FALSE,
                     pcutoff = 0.05) {
  .Deprecated("fit_curves", package = "gDRutils")
  fit_curves(
    df_ = df_,
    e_0 = e_0,
    GR_0 = GR_0,
    n_point_cutoff = n_point_cutoff,
    range_conc = range_conc,
    force_fit = force,
    pcutoff = pcutoff
  )
}


