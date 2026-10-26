####
# Public API
####

#' apply_fit
#'
#' Generic layer for applying a user-supplied fit function across all
#' (row x column x slicing_col) triplets in a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} and persisting
#' the results into a named output assay.
#'
#' @details
#' The experiment type is declared via \code{data_type}, which resolves a
#' built-in slicing profile (default \code{slicing_cols}, \code{slicing_values},
#' and \code{input_assay}).  Any profile field can be overridden with the
#' corresponding argument.
#'
#' An optional \code{summary_fn} may be provided.  It is called once per
#' (row x column) cell and receives all rows that \code{fit_fn} produced for
#' that cell, returning a single summary row written to \code{summary_assay}.
#' This is the right place for aggregated metrics (e.g. mean synergy across
#' normalization types) that span multiple slices.
#'
#' Use the pipe to apply several custom fits to the same SE:
#' \preformatted{
#'   se |>
#'     apply_fit(bliss_fn, "combination",
#'                      output_assay = "custom_bliss", ...) |>
#'     apply_fit(musyc_fn,  "combination",
#'                      output_assay = "musyc_params",
#'                      summary_fn   = musyc_summary_fn,
#'                      summary_assay = "musyc_summary", ...)
#' }
#'
#' Each call records the recipe of the fit in \code{metadata(se)$provenance}:
#' which function produced it (package and version, or its source code when it
#' does not come from a package), the arguments in \code{fit_args} and the call
#' settings.  Pass settings through \code{fit_args} rather than wrapping the
#' function in a closure, so they are visible in the recipe.  See
#' \code{\link{get_fit_provenance}} and \code{\link{check_reproducibility}}.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @param fit_fn function(\code{data.table}) returning a named list or
#'   single-row data.frame.  Receives data already filtered to one
#'   \code{slicing_cols} combination.
#' @param data_type one of \code{"single-agent"}, \code{"co-dilution"},
#'   \code{"combination"}, \code{"time-course"}.  Resolves profile defaults.
#' @param slicing_cols character vector of column name(s) used to split each
#'   BumpyMatrix cell into sub-experiments.  \code{NULL} uses the profile
#'   default.
#' @param slicing_values character vector of values to iterate over in
#'   \code{slicing_cols}.  \code{NULL} iterates all unique values found in
#'   each cell.
#' @param input_assay name of the source BumpyMatrix assay.  \code{NULL} uses
#'   the profile default.
#' @param output_assay name of the assay to write fit results into.  Required.
#'   Any assay name is accepted, including the native \code{"Metrics"} assay.
#'   With \code{merge = "merge"} (default), existing rows keyed by
#'   \code{fit_source} are replaced while all other rows are preserved —
#'   safe for co-existing alongside native gDR metrics.
#'   With \code{merge = "replace"}, the entire assay is overwritten.
#' @param summary_fn optional function(\code{data.table}) → named list called
#'   once per (row x column) cell on all rows returned by \code{fit_fn} for
#'   that cell.  Requires \code{summary_assay}.
#' @param summary_assay name of the assay to write summary results into.
#'   Required when \code{summary_fn} is provided.
#' @param merge \code{"merge"} (idempotent upsert keyed by \code{fit_source}
#'   + \code{slicing_cols}) or \code{"replace"} (overwrite the whole assay).
#' @param on_error \code{"warn"} (skip failed cells with a warning) or
#'   \code{"stop"} (propagate the error).
#' @param fit_source character string recorded in the \code{fit_source} column
#'   of every output row.  Forms part of the upsert key.
#' @param fit_args named list of extra arguments passed to \code{fit_fn} after
#'   the data, e.g. \code{list(pcutoff = 0.1)}.  Recorded in the provenance.
#'
#' @return updated \code{SummarizedExperiment} with \code{output_assay} (and
#'   optionally \code{summary_assay}) added or updated, and the recipe of the
#'   fit in \code{metadata(se)$provenance}.
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' mean_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
#' se_out <- apply_fit(
#'   se, mean_fn, "single-agent",
#'   output_assay = "custom_mean", fit_source = "demo"
#' )
#'
#' @keywords metrics
#' @export
#'
apply_fit <- function(se,
                      fit_fn,
                      data_type,
                      slicing_cols = NULL,
                      slicing_values = NULL,
                      input_assay = NULL,
                      output_assay,
                      summary_fn = NULL,
                      summary_assay = NULL,
                      merge = "merge",
                      on_error = "warn",
                      fit_source,
                      fit_args = list()) {

  # names as the caller wrote them, for the recipe; read before anything is evaluated
  fit_fn_name <- .call_name(substitute(fit_fn))
  summary_fn_name <- .call_name(substitute(summary_fn))

  checkmate::assert_string(data_type)
  profile <- get_fit_profile(data_type)  # errors with list of valid profiles if unknown

  if (is.null(slicing_cols)) slicing_cols <- profile$slicing_cols
  if (is.null(slicing_values)) slicing_values <- profile$slicing_values
  if (is.null(input_assay)) input_assay <- profile$input_assay

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_function(fit_fn)
  checkmate::assert_character(slicing_cols, min.len = 1L)
  checkmate::assert_character(slicing_values, null.ok = TRUE, min.len = 1L)
  checkmate::assert_string(input_assay)
  checkmate::assert_string(output_assay)
  checkmate::assert_choice(merge, c("merge", "replace"))
  checkmate::assert_choice(on_error, c("warn", "stop"))
  checkmate::assert_string(fit_source)
  checkmate::assert_list(fit_args, names = "unique")
  if (!is.null(summary_fn)) {
    checkmate::assert_function(summary_fn)
    checkmate::assert_string(summary_assay)
  }

  if (!input_assay %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("%s assay is required", input_assay))
  }
  input_hash <- .assay_hash(se, input_assay)

  tmp_assay <- "__custom_fit_tmp__"
  wrapper_fn <- .make_fit_wrapper(
    fit_fn, slicing_cols, slicing_values, on_error, fit_source, fit_args
  )

  se_out <- apply_bumpy_function(
    se = se,
    FUN = wrapper_fn,
    req_assay_name = input_assay,
    out_assay_name = tmp_assay
  )

  if (tmp_assay %in% SummarizedExperiment::assayNames(se_out)) {
    new_df <- BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se_out, tmp_assay),
      row.field = "row", column.field = "column"
    )
    new_metrics <- data.table::as.data.table(new_df)

    if (NROW(new_metrics) > 0L) {
      upsert_key <- c("fit_source", slicing_cols)
      se <- persist_fit_assay(se, new_metrics, merge, output_assay, "row", "column", upsert_key)

      if (!is.null(summary_fn)) {
        se <- .apply_summary_fn(
          se, new_metrics, summary_fn, summary_assay, merge, fit_source, on_error
        )
      }

      recipe <- .fit_recipe(
        called_by = "apply_fit",
        fn_record = .describe_fn(fit_fn, fit_fn_name),
        fit_args = fit_args,
        summary_record = if (!is.null(summary_fn)) .describe_fn(summary_fn, summary_fn_name),
        data_type = data_type,
        input_assay = input_assay,
        output_assay = output_assay,
        summary_assay = summary_assay,
        slicing_cols = slicing_cols,
        slicing_values = slicing_values,
        merge = merge,
        input_hash = input_hash,
        n_rows = NROW(new_metrics)
      )
      se <- .record_fit(se, output_assay, fit_source, recipe,
                        replace = merge == "replace", summary_assay = summary_assay)
    }
  }

  se
}


#' apply_fit_to_se
#'
#' Apply a user-supplied fit function per (drug x cell line x normalization_type)
#' triplet and persist results into the Metrics assay.
#'
#' This is a convenience wrapper around \code{\link{apply_fit}} for the
#' single-agent use case with the standard Metrics assay.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   an Averaged assay
#' @param fit_fn function taking a \code{data.table} and returning a named list
#' @param normalization_types character vector of types to iterate
#' @param averaged_assay name of the input assay
#' @param metrics_assay name of the output assay
#' @param merge \code{"merge"} (idempotent upsert) or \code{"replace"}
#' @param on_error \code{"warn"} (skip + warning) or \code{"stop"}
#' @param fit_source string recorded in the \code{fit_source} column
#' @param fit_args named list of extra arguments passed to \code{fit_fn}
#'
#' @return updated \code{SummarizedExperiment}
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' simple_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
#' se_out <- apply_fit_to_se(se, simple_fn, fit_source = "demo")
#'
#' @keywords metrics
#' @export
#'
apply_fit_to_se <- function(se,
                            fit_fn,
                            normalization_types = c("GR", "RV"),
                            averaged_assay = "Averaged",
                            metrics_assay = "Metrics",
                            merge = "merge",
                            on_error = "warn",
                            fit_source = "custom",
                            fit_args = list()) {

  # Warn for missing standard metric columns before delegating
  # (apply_fit does not impose column expectations on custom assays)
  result <- apply_fit(
    se = se,
    fit_fn = fit_fn,
    data_type = "single-agent",
    slicing_values = normalization_types,
    input_assay = averaged_assay,
    output_assay = metrics_assay,
    merge = merge,
    on_error = on_error,
    fit_source = fit_source,
    fit_args = fit_args
  )

  # Warn only for the standard Metrics assay — custom assays (e.g. "musyc_params")
  # intentionally use a different column schema and should not be checked here.
  if (metrics_assay == "Metrics" &&
      metrics_assay %in% SummarizedExperiment::assayNames(result)) {
    new_df <- BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(result, metrics_assay),
      row.field = "row", column.field = "column"
    )
    missing <- setdiff(
      get_header("response_metrics"),
      names(new_df)
    )
    if (length(missing) > 0L) {
      warning(sprintf(
        "fit_fn output missing columns: %s",
        toString(missing)
      ))
    }
  }

  result
}


####
# Multi-fit (single-pass)
####

#' apply_fits
#'
#' Apply multiple fit functions in a \strong{single pass} over each BumpyMatrix
#' cell, writing each function's results into its own named output assay.
#'
#' @details
#' \code{apply_fits()} is the performance-efficient alternative to
#' chaining multiple \code{\link{apply_fit}} calls when two or more fit
#' functions operate on the \strong{same input assay}.  Instead of unsplitting
#' the BumpyMatrix K times (once per function), it traverses each cell once and
#' applies all functions in that single pass.
#'
#' Use this when:
#' \itemize{
#'   \item You have two or more independent fit functions on the same data
#'     (e.g. Bliss + HSS on combination data).
#'   \item A single fit function produces results for multiple output assays
#'     (shared pre-computation pattern — see below).
#' }
#'
#' \subsection{Independent fit functions (most common)}{
#' Names of \code{fit_fns} become the output assay names:
#' \preformatted{
#'   apply_fits(
#'     combo_se,
#'     fit_fns = list(
#'       custom_bliss = bliss_fit_fn,
#'       custom_hss   = hss_fit_fn
#'     ),
#'     data_type  = "combination",
#'     fit_source = "synergy"
#'   )
#' }
#' }
#'
#' \subsection{Shared pre-computation (advanced)}{
#' A single function can return a \strong{named list of named lists} to write
#' multiple assays while computing expensive intermediates only once:
#' \preformatted{
#'   bliss_and_hss <- function(dt) {
#'     sa_curves <- fit_hill_curves(dt)  # expensive — done once
#'     list(
#'       custom_bliss = list(bliss_score = compute_bliss(sa_curves, dt)),
#'       custom_hss   = list(hss_score   = compute_hss(sa_curves, dt))
#'     )
#'   }
#'   apply_fits(
#'     combo_se,
#'     fit_fns    = list(bliss_and_hss = bliss_and_hss),
#'     output_assay_map = c(bliss_and_hss = NA),  # ignored; keys from return value
#'     ...
#'   )
#' }
#' The multi-output pattern is detected automatically when a fit function
#' returns a named list whose values are themselves named lists.  Each top-level
#' name maps to an assay; the inner named list provides the row columns.
#' }
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @param fit_fns named list of functions.  Each name becomes an output assay
#'   name; each value is a function(\code{data.table}) → named list (or a named
#'   list of named lists for the shared pre-computation pattern).
#' @param data_type one of \code{"single-agent"}, \code{"co-dilution"},
#'   \code{"combination"}, \code{"time-course"}.
#' @param slicing_cols character vector; \code{NULL} uses profile default.
#' @param slicing_values character vector; \code{NULL} uses profile default.
#' @param input_assay string; \code{NULL} uses profile default.
#' @param merge \code{"merge"} or \code{"replace"}.
#' @param on_error \code{"warn"} or \code{"stop"}.
#' @param fit_source character string stamped as \code{fit_source} in every
#'   output row.
#' @param fit_args named list of argument lists, keyed by names of
#'   \code{fit_fns}, e.g. \code{list(custom_bliss = list(q = 0.9))}.  Each list
#'   is passed to its function after the data and recorded in the provenance.
#'
#' @return updated \code{SummarizedExperiment} with one new (or updated) assay
#'   per entry in \code{fit_fns}, and the recipe of each fit in
#'   \code{metadata(se)$provenance} (see \code{\link{get_fit_provenance}}).
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' fn_a <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
#' fn_b <- function(dt) list(x_sd = sd(dt$x, na.rm = TRUE))
#' se_out <- apply_fits(
#'   se,
#'   fit_fns = list(mean_metrics = fn_a, sd_metrics = fn_b),
#'   data_type = "single-agent",
#'   fit_source = "demo"
#' )
#'
#' @keywords metrics
#' @export
#'
apply_fits <- function(se,
                       fit_fns,
                       data_type,
                       slicing_cols = NULL,
                       slicing_values = NULL,
                       input_assay = NULL,
                       merge = "merge",
                       on_error = "warn",
                       fit_source,
                       fit_args = list()) {

  # names as the caller wrote them, for the recipes
  fn_names <- .list_call_names(substitute(fit_fns))

  checkmate::assert_string(data_type)
  profile <- get_fit_profile(data_type)

  if (is.null(slicing_cols)) slicing_cols <- profile$slicing_cols
  if (is.null(slicing_values)) slicing_values <- profile$slicing_values
  if (is.null(input_assay)) input_assay <- profile$input_assay

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_list(fit_fns, min.len = 1L, types = "function", names = "unique")
  checkmate::assert_character(slicing_cols, min.len = 1L)
  checkmate::assert_character(slicing_values, null.ok = TRUE, min.len = 1L)
  checkmate::assert_string(input_assay)
  checkmate::assert_choice(merge, c("merge", "replace"))
  checkmate::assert_choice(on_error, c("warn", "stop"))
  checkmate::assert_string(fit_source)
  checkmate::assert_list(fit_args, types = "list", names = "unique")
  checkmate::assert_subset(names(fit_args), names(fit_fns))

  if (!input_assay %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("%s assay is required", input_assay))
  }

  input_hash <- .assay_hash(se, input_assay)
  out_assay_names <- names(fit_fns)
  # function that wrote each assay: with the shared pre-computation pattern one
  # function can fill several assays
  producers <- character(0)
  slice_col <- slicing_cols[1L]

  # Unsplit the BumpyMatrix ONCE into a flat data.table, then split by
  # (row × column) — avoids repeated S4 dispatch per cell that dominates
  # wall-clock time when fit_fn is fast (e.g. Bliss, HSS, summary stats).
  asy <- SummarizedExperiment::assay(se, input_assay)
  dt_all <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(asy, row.field = "row", column.field = "column")
  )
  # BumpyMatrix may store slicing_col (e.g. normalization_type) as a factor;
  # coerce to character so downstream comparisons and labels are always strings.
  if (slice_col %in% names(dt_all) && is.factor(dt_all[[slice_col]])) {
    dt_all[, (slice_col) := as.character(get(slice_col))]
  }
  cell_list <- split(dt_all, by = c("row", "column"), sorted = FALSE)

  # Accumulate per-assay results: list(assay_name -> list of row data.tables)
  assay_results <- stats::setNames(
    lapply(out_assay_names, function(nm) list()),
    out_assay_names
  )

  for (cell_dt in cell_list) {
    r <- cell_dt[["row"]][1L]
    cc <- cell_dt[["column"]][1L]

    vals <- if (!is.null(slicing_values)) {
      slicing_values
    } else if (slice_col %in% names(cell_dt)) {
      as.character(unique(cell_dt[[slice_col]]))
    } else {
      NA_character_
    }

    for (val in vals) {
      sub_dt <- if (slice_col %in% names(cell_dt) && !is.na(val)) {
        cell_dt[cell_dt[[slice_col]] == val, ]
      } else {
        cell_dt
      }
      if (NROW(sub_dt) == 0L) next
      label <- if (!is.na(val)) sprintf("%s=%s", slice_col, val) else "all"
      # a multi-output function is listed once per assay it fills, but one call
      # already returns rows for all of them; calling it again duplicated the rows
      multi_done <- list()

      for (fn_nm in out_assay_names) {
        fn <- fit_fns[[fn_nm]]
        if (any(vapply(multi_done, identical, logical(1), fn))) next
        result_dt <- tryCatch({
          res <- .call_fit_fn(fn, sub_dt, fit_args[[fn_nm]])
          # Detect multi-output pattern: named list of named lists
          if (is.list(res) && length(res) > 0L &&
              is.list(res[[1L]]) && !is.null(names(res))) {
            for (assay_nm in names(res)) {
              row_dt <- data.table::as.data.table(as.list(res[[assay_nm]]))
              row_dt[["fit_source"]] <- fit_source
              row_dt[["row"]] <- r
              row_dt[["column"]] <- cc
              if (!is.na(val) && !slice_col %in% names(row_dt)) {
                row_dt[[slice_col]] <- val
              }
              if (assay_nm %in% out_assay_names) {
                assay_results[[assay_nm]] <- c(
                  assay_results[[assay_nm]], list(row_dt)
                )
                if (is.na(producers[assay_nm])) producers[[assay_nm]] <- fn_nm
              }
            }
            multi_done <- c(multi_done, list(fn))
            next
          }
          data.table::as.data.table(as.list(res))
        }, error = function(e) {
          if (on_error == "stop") stop(e)
          warning(sprintf("fit_fn '%s' failed for %s row=%s col=%s: %s",
                          fn_nm, label, r, cc, conditionMessage(e)))
          NULL
        })

        if (is.null(result_dt)) next

        result_dt[["fit_source"]] <- fit_source
        result_dt[["row"]] <- r
        result_dt[["column"]] <- cc
        if (!is.na(val) && !slice_col %in% names(result_dt)) {
          result_dt[[slice_col]] <- val
        }
        assay_results[[fn_nm]] <- c(assay_results[[fn_nm]], list(result_dt))
        if (is.na(producers[fn_nm])) producers[[fn_nm]] <- fn_nm
      }
    }
  }

  # Persist each assay's accumulated rows
  upsert_key <- c("fit_source", slicing_cols)
  for (assay_nm in out_assay_names) {
    rows <- assay_results[[assay_nm]]
    if (length(rows) == 0L) next
    merged_dt <- data.table::rbindlist(rows, fill = TRUE)
    se <- persist_fit_assay(se, merged_dt, merge, assay_nm, "row", "column", upsert_key)

    producer <- producers[[assay_nm]]
    recipe <- .fit_recipe(
      called_by = "apply_fits",
      fn_record = .describe_fn(fit_fns[[producer]], fn_names[[producer]]),
      fit_args = fit_args[[producer]],
      summary_record = NULL,
      data_type = data_type,
      input_assay = input_assay,
      output_assay = assay_nm,
      summary_assay = NULL,
      slicing_cols = slicing_cols,
      slicing_values = slicing_values,
      merge = merge,
      input_hash = input_hash,
      n_rows = NROW(merged_dt)
    )
    se <- .record_fit(se, assay_nm, fit_source, recipe, replace = merge == "replace")
  }

  se
}


####
# Reference fit functions
####

#' fit_drug_response_metrics
#'
#' Reference fit function replicating standard gDR Hill curve fitting
#' (3-parameter model, \code{x_0} fixed at 1).  Produces output
#' column-compatible with \code{fit_SE()} in gDRcore and
#' \code{\link{logisticFit}}, including \code{p_value},
#' \code{rss}, \code{x_AOC_range}, \code{x_max}, and \code{x_sd_avg}.
#' For use with \code{\link{apply_fit_to_se}} or \code{\link{apply_fit}}
#' on single-agent or time-course data.
#'
#' @param avg_dt \code{data.table} of averaged data for one
#'   (drug x cell line x normalization_type) triplet.  Should contain
#'   columns \code{Concentration}, the response column named by \code{x_col},
#'   \code{normalization_type}, and optionally \code{x_std}
#'   (used for \code{x_sd_avg}).
#' @param x_col string; name of the response column in \code{avg_dt}.
#'   Defaults to \code{"x"} (the standard \code{Averaged} assay column name).
#'   Set to e.g. \code{"NormalizedGrowthRate"} when fitting time-course
#'   growth-rate data.
#' @param capping_fold numeric capping fold passed to
#'   \code{\link{cap_xc50}}
#' @param range_conc numeric vector of length 2 specifying the concentration
#'   range used to compute \code{x_AOC_range}. Default \code{c(5e-3, 5)},
#'   matching the \code{fit_SE()} default.
#' @param pcutoff numeric p-value cutoff for the fit. Default \code{0.05}.
#' @param n_point_cutoff integer minimum number of concentration points required
#'   to attempt fitting. Default \code{4L}.
#' @param force_fit logical; if \code{TRUE}, force fitting even when the number
#'   of points is below \code{n_point_cutoff}. Default \code{FALSE}.
#' @param cap numeric; upper cap on \code{xc50} relative to the concentration
#'   range. Default \code{0.1}.
#'
#' @return Named list of fit metrics fully compatible with the standard gDR
#'   \code{Metrics} assay schema: \code{ec50}, \code{xc50}, \code{h},
#'   \code{x_inf}, \code{x_0}, \code{r2}, \code{rss}, \code{p_value},
#'   \code{x_mean}, \code{x_AOC}, \code{x_AOC_range}, \code{x_max},
#'   \code{x_sd_avg}, \code{N_conc}, \code{maxlog10Concentration},
#'   \code{fit_type}, \code{normalization_type}. \code{fit_type} is one of:
#'   \itemize{
#'     \item \code{"DRC3pHillFitModelFixS0"} (\code{"DRC4pHillFitModel"} for the
#'       4-parameter variant): the curve was fitted;
#'     \item \code{"DRCConstantFitResult"}: all responses are equal, or the fit is not
#'       significant at \code{pcutoff} (its \code{p_value} and \code{rss} are kept);
#'       the curve is flat at the mean response;
#'     \item \code{"DRCTooFewPointsToFit"}: fewer than \code{n_point_cutoff}
#'       concentrations with a response; no curve parameters;
#'     \item \code{"DRCInvalidFitResult"}: the solver failed; no curve parameters.
#'   }
#'   The same labels as \code{\link{logisticFit}} gives for the same data.
#'
#' @seealso \code{\link{fit_drug_response_metrics_4p}} for a 4-parameter
#'   variant that fits \code{x_0} freely.
#'
#' @examples
#' dt <- data.table::data.table(
#'   Concentration = c(0.001, 0.01, 0.1, 1, 10),
#'   x = c(0.95, 0.8, 0.5, 0.2, 0.1),
#'   normalization_type = "RV"
#' )
#' fit_drug_response_metrics(dt)
#'
#' # Custom response column (e.g. time-course growth rates):
#' dt2 <- data.table::data.table(
#'   Concentration = c(0.001, 0.01, 0.1, 1, 10),
#'   NormalizedGrowthRate = c(0.98, 0.82, 0.55, 0.21, 0.09),
#'   normalization_type = "GR"
#' )
#' fit_drug_response_metrics(dt2, x_col = "NormalizedGrowthRate")
#'
#' @export
fit_drug_response_metrics <- function(avg_dt, x_col = "x",
                                      capping_fold = 5,
                                      range_conc = c(5e-3, 5),
                                      pcutoff = 0.05,
                                      n_point_cutoff = 4L,
                                      force_fit = FALSE,
                                      cap = 0.1) {
  .fit_drug_response_metrics_impl(avg_dt, x_0 = 1,
                                  x_col = x_col,
                                  capping_fold = capping_fold,
                                  range_conc = range_conc,
                                  pcutoff = pcutoff,
                                  n_point_cutoff = n_point_cutoff,
                                  force_fit = force_fit,
                                  cap = cap)
}


#' fit_drug_response_metrics_4p
#'
#' Variant of \code{\link{fit_drug_response_metrics}} using a 4-parameter
#' Hill model — \code{x_0} is fitted freely rather than fixed at 1.
#' Useful when the upper asymptote is expected to deviate from 1 (e.g.
#' partial agonists, non-standard assay readouts).
#'
#' @inheritParams fit_drug_response_metrics
#'
#' @return Named list of fit metrics.  \code{fit_type} will be
#'   \code{"DRC4pHillFitModel"} on success.
#'
#' @seealso \code{\link{fit_drug_response_metrics}} for the default
#'   3-parameter variant that matches \code{fit_SE()} output.
#'
#' @examples
#' dt <- data.table::data.table(
#'   Concentration = c(0.001, 0.01, 0.1, 1, 10),
#'   x = c(0.95, 0.8, 0.5, 0.2, 0.1),
#'   normalization_type = "RV"
#' )
#' fit_drug_response_metrics_4p(dt)
#'
#' @export
fit_drug_response_metrics_4p <- function(avg_dt, x_col = "x",
                                         capping_fold = 5,
                                         range_conc = c(5e-3, 5),
                                         pcutoff = 0.05,
                                         n_point_cutoff = 4L,
                                         force_fit = FALSE,
                                         cap = 0.1) {
  .fit_drug_response_metrics_impl(avg_dt, x_0 = NA_real_,
                                  x_col = x_col,
                                  capping_fold = capping_fold,
                                  range_conc = range_conc,
                                  pcutoff = pcutoff,
                                  n_point_cutoff = n_point_cutoff,
                                  force_fit = force_fit,
                                  cap = cap)
}

# What the reference fits do, copied into the recipe apply_fit() records.
fit_drug_response_metrics <- set_fit_fn_info(
  fit_drug_response_metrics,
  model = "3-parameter log-logistic (drc::LL.3u), upper asymptote fixed at 1",
  description = paste(
    "Hill curve of the response against concentration, fitted per normalization type.",
    "Falls back to a constant fit when fewer than n_point_cutoff concentrations are",
    "available or, unless force_fit = TRUE, when the F-test p-value is at least pcutoff."
  ),
  reference = "Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021"
)

fit_drug_response_metrics_4p <- set_fit_fn_info(
  fit_drug_response_metrics_4p,
  model = "4-parameter log-logistic (drc::LL.4), upper asymptote fitted",
  description = paste(
    "Hill curve of the response against concentration with a free upper asymptote,",
    "fitted per normalization type. Falls back to a constant fit like",
    "fit_drug_response_metrics()."
  ),
  reference = "Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021"
)


# Fit configuration for one normalization type. The values, the fallback and the
# rationale for each live in the fit_config block of inst/extdata/fit_profiles.json —
# adding or retuning a type is a data edit, not a code change. Keyed by normalization
# type rather than nested per profile because fit_fn receives only the data.table,
# never the resolved profile.
#' @keywords internal
.fit_config_for <- function(norm_type) {
  if (length(norm_type) != 1L || is.na(norm_type)) {
    stop("normalization_type must be a single non-missing value, got: ",
         toString(norm_type))
  }
  # the config is read from JSON on first use; without this, a caller that reaches the
  # fitting layer directly rather than through apply_fit() sees an empty environment
  .ensure_fit_profiles()
  cfg <- .fit_config_env[[norm_type]] %||% .fit_config_env[[.fit_config_env[["_default"]]]]
  if (is.null(cfg)) {
    stop("no fit configuration for normalization_type '", norm_type,
         "' and no usable default in fit_profiles.json")
  }
  lapply(cfg, as.numeric)
}


# Starting value and lower bounds for the asymptote parameters, per normalization type.
#' @keywords internal
.fit_bounds_for <- function(norm_type) {
  .fit_config_for(norm_type)
}


# Response level at which a curve counts as having crossed 50%. Hardcoded as 0.5 in two
# places before; carried in configuration so the value is visible per normalization type
# — see the note in fit_profiles.json for why it is not self-evident on the NGR scale.
#' @keywords internal
.xc50_threshold_for <- function(norm_type) {
  .fit_config_for(norm_type)[["xc50_threshold"]] %||% 0.5
}


#' @keywords internal
.fit_drug_response_metrics_impl <- function(avg_dt, x_0 = 1, x_col = "x",
                                            capping_fold = 5,
                                            range_conc = c(5e-3, 5),
                                            pcutoff = 0.05,
                                            n_point_cutoff = 4L,
                                            force_fit = FALSE,
                                            cap = 0.1) {
  checkmate::assert_string(x_col)
  if (!x_col %in% names(avg_dt)) {
    stop(sprintf("Column '%s' not found in avg_dt. Available: %s",
                 x_col, toString(names(avg_dt))))
  }
  norm_type <- as.character(avg_dt$normalization_type[1])
  conc_col <- get_env_identifiers("concentration")
  conc <- avg_dt[[conc_col]]
  x <- avg_dt[[x_col]]
  x_std <- if ("x_std" %in% names(avg_dt)) avg_dt$x_std else rep(NA_real_, length(x))

  keep <- !is.na(x) & !is.na(conc) & conc != 0
  x_std <- x_std[keep]
  x <- x[keep]
  conc <- conc[keep]

  # Cap x values as logisticFit does (pmin to x_0 + cap), before anything is decided
  # or summarized from them: a flat response above the cap is then one value, and
  # x_mean and x_max of every outcome use the capped values
  checkmate::assert_number(cap, lower = 0)
  if (is.na(x_0)) {
    x_cap_limit <- 1 + cap
  } else {
    x_cap_limit <- x_0 + cap
  }
  x <- pmin(x, x_cap_limit)

  # The three ways a curve ends up without a dose-response fit, as in logisticFit():
  # - fewer than n_point_cutoff concentrations (none at all included): nothing to fit
  #   a curve to, so no curve parameters -> DRCTooFewPointsToFit
  # - one distinct response: nothing to fit (drc would stop on it) -> DRCConstantFitResult
  # - the solver fails -> DRCInvalidFitResult, further down
  # (is.na(NaN) is TRUE, so NaN responses were dropped with the NAs above)
  if (length(unique(conc)) < n_point_cutoff) {
    .too_few_fit_result(norm_type, x, conc, x_std)
  } else if (length(unique(x)) == 1L) {
    .constant_fit_result(norm_type, x, conc, x_std, x_0, range_conc)
  } else {

  x_sd_avg <- mean(x_std, na.rm = TRUE)         # matches logisticFit line 239
  N_conc <- length(unique(conc))
  maxlog10Conc <- log10(max(conc, na.rm = TRUE))

  # x_max: min of the two highest-concentration x values (matches .calculate_x_max)
  conc_ord <- order(conc)
  x_ordered <- x[conc_ord]
  n_x <- length(x_ordered)
  if (n_x >= 2L) {
    x_max <- min(x_ordered[(n_x - 1L):n_x])
  } else {
    x_max <- x_ordered[n_x]
  }

  bounds <- .fit_bounds_for(norm_type)
  x_inf_prior <- bounds$x_inf_prior
  lower_x_inf <- bounds$lower_x_inf
  lower_x_0 <- bounds$lower_x_0
  med_conc <- stats::median(conc)
  min_conc <- min(conc)
  max_conc <- max(conc)
  controls <- drc::drmc(relTol = 1e-04, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)

  # Compute mean on capped x (must be after capping so rss1 uses consistent data)
  x_mean_obs <- mean(x, na.rm = TRUE)

  four_param <- is.na(x_0)
  fit_param <- c("h", "x_inf", "x_0", "ec50")

  if (four_param) {
    fct <- drc::LL.4(names = fit_param)
    start <- c(2, x_inf_prior, 1, med_conc)
    lowerl <- c(0.1, lower_x_inf, lower_x_0, min_conc / 10)
    upperl <- c(5, 1, 1 + cap, max_conc * 10)
  } else {
    fit_param_3p <- fit_param[-3]   # drop x_0
    fct <- drc::LL.3u(upper = x_0, names = fit_param_3p)
    start <- c(2, x_inf_prior, med_conc)
    lowerl <- c(0.1, lower_x_inf, min_conc / 10)
    upperl <- c(5, min(x_0 + cap, 1), max_conc * 10)
  }

  # See the note in .prepareFitModel() - drc leaks na.action into the global
  # option, and a qualified name there breaks every later stats::lm().
  old_na_action <- getOption("na.action")
  on.exit(options(na.action = old_na_action), add = TRUE)

  fit <- tryCatch(
    drc::drm(x ~ conc,
             data = data.table::data.table(x = x, conc = conc),
             logDose = NULL,
             fct = fct,
             start = start,
             lowerl = lowerl,
             upperl = upperl,
             control = controls,
             # Unqualified on purpose - see the note in .prepareFitModel().
             na.action = na.omit),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    coefs <- stats::coef(fit)
    rss <- sum(stats::residuals(fit)^2, na.rm = TRUE)
    rss1 <- sum((x - x_mean_obs)^2)
    r2 <- 1 - rss / rss1

    # F-test p-value (matches .calculate_f_pval in gDRutils/fit_curves.R)
    if (four_param) {
      n_param <- 4L
    } else {
      n_param <- 3L
    }
    df1 <- n_param - 1L
    df2 <- length(x) - n_param + 1L
    f_value <- ((rss1 - rss) / df1) / (rss / df2)
    p_value <- stats::pf(f_value, df1, df2, lower.tail = FALSE)

    # pcutoff: a fit that is not statistically significant falls back to a constant
    # fit, which keeps the p-value and rss of the rejected fit, as logisticFit() does
    if (!force_fit && !is.na(p_value) && p_value >= pcutoff) {
      .constant_fit_result(norm_type, x, conc, x_std, x_0, range_conc,
                           p_value = p_value, rss = rss)
    } else {
      # Helper: mean of model predictions over a concentration range
      .predict_mean <- function(model, lo, hi) {
        mean(stats::predict(model, data.table::data.table(
          conc = 10^seq(log10(lo), log10(hi), length.out = 101L)
        )), na.rm = TRUE)
      }

      # x_mean: model-predicted mean over observed range (matches logisticFit)
      x_mean_model <- .predict_mean(fit, min(conc), max(conc))

      # x_AOC_range: model-predicted mean over full range_conc (matches logisticFit line 290)
      # logisticFit integrates over [range_conc[1], range_conc[2]] without
      # clamping to the observed concentration range — we must do the same for parity.
      if (range_conc[1] < range_conc[2]) {
        x_AOC_range <- 1 - .predict_mean(fit, range_conc[1], range_conc[2])
      } else {
        x_AOC_range <- 1 - x_mean_model
      }

      if (four_param) {
        ec50_val <- coefs[["ec50:(Intercept)"]]
        x_0_val <- coefs[["x_0:(Intercept)"]]
        x_inf_val <- coefs[["x_inf:(Intercept)"]]
        fit_type <- "DRC4pHillFitModel"
      } else {
        ec50_val <- coefs[["ec50:(Intercept)"]]
        x_0_val <- x_0
        x_inf_val <- coefs[["x_inf:(Intercept)"]]
        fit_type <- "DRC3pHillFitModelFixS0"
      }
      list(
        normalization_type = norm_type,
        x_mean = x_mean_model,
        x_AOC = 1 - x_mean_model,
        x_AOC_range = x_AOC_range,
        x_max = x_max,
        x_sd_avg = x_sd_avg,
        N_conc = N_conc,
        maxlog10Concentration = maxlog10Conc,
        ec50 = ec50_val,
        xc50 = {
          # xc50 = concentration at which response = 0.5
          # (not ec50 — matches .calculate_xc50 + .extendWithXc50 in gDRutils/fit_curves.R)
          xc50_raw <- tryCatch(
            predict_conc_from_efficacy(0.5, x_inf_val, x_0_val, ec50_val, coefs[1]),
            error = function(e) NA_real_
          )
          if (is.na(xc50_raw)) {
            cap_xc50(ec50_val, max_conc,
                               min_conc = min(conc[conc > 0]),
                               capping_fold = capping_fold)
          } else {
            cap_xc50(xc50_raw, max_conc,
                               min_conc = min(conc[conc > 0]),
                               capping_fold = capping_fold)
          }
        },
        h = coefs[1],
        r2 = r2,
        rss = rss,
        p_value = p_value,
        x_0 = x_0_val,
        x_inf = x_inf_val,
        fit_type = fit_type
      )
    }
  } else {
    list(
      normalization_type = norm_type,
      x_mean = x_mean_obs,
      x_AOC = 1 - x_mean_obs,
      x_AOC_range = 1 - x_mean_obs,
      x_max = x_max,
      x_sd_avg = x_sd_avg,
      N_conc = N_conc,
      maxlog10Concentration = maxlog10Conc,
      ec50 = NA_real_,
      xc50 = .estimate_xc50_fallback(x, .xc50_threshold_for(norm_type)),
      h = NA_real_,
      r2 = NA_real_,
      rss = NA_real_,
      p_value = NA_real_,
      x_0 = if (four_param) NA_real_ else x_0,
      x_inf = NA_real_,
      fit_type = "DRCInvalidFitResult"
    )
  }
  } # end else (enough conc points)
}


####
# Internal helpers
####

#' Build a per-cell wrapper for apply_bumpy_function
#'
#' Iterates over each unique combination of slicing_col values, filters the
#' cell data.table, calls fit_fn, stamps fit_source and slicing column values,
#' and rbinds all results.  Respects on_error semantics.
#'
#' @param fit_fn user-supplied fit function
#' @param slicing_cols character vector of column names to slice by
#' @param slicing_values character vector of values to iterate, or NULL for all
#' @param on_error \code{"warn"} or \code{"stop"}
#' @param fit_source fit_source label
#' @param fit_args named list of extra arguments passed to \code{fit_fn}
#'
#' @return function(DFrame) -> data.table
#'
#' @keywords internal
.make_fit_wrapper <- function(fit_fn, slicing_cols, slicing_values,
                              on_error, fit_source, fit_args = list()) {
  function(avg_df) {
    avg_dt <- data.table::as.data.table(avg_df)
    results <- list()

    # Determine the actual values to iterate for the primary slicing column.
    # Multi-column slicing is deferred to a future iteration; for the POC the
    # profiles always supply a single slicing column.
    slice_col <- slicing_cols[1L]

    # BumpyMatrix may store slicing columns (e.g. normalization_type) as factors.
    # Coerce to character so subset comparisons and output labels are always strings.
    if (slice_col %in% names(avg_dt) && is.factor(avg_dt[[slice_col]])) {
      avg_dt[, (slice_col) := as.character(get(slice_col))]
    }

    vals <- if (!is.null(slicing_values)) {
      slicing_values
    } else if (slice_col %in% names(avg_dt)) {
      unique(avg_dt[[slice_col]])
    } else {
      NA_character_   # no slicing column present → one call on whole cell
    }

    for (val in vals) {
      sub_dt <- if (slice_col %in% names(avg_dt) && !is.na(val)) {
        avg_dt[avg_dt[[slice_col]] == val, ]
      } else {
        avg_dt
      }

      if (NROW(sub_dt) == 0L) next

      label <- if (!is.na(val)) sprintf("%s=%s", slice_col, val) else "all"
      result_dt <- tryCatch({
        res <- .call_fit_fn(fit_fn, sub_dt, fit_args)
        data.table::as.data.table(as.list(res))
      }, error = function(e) {
        if (on_error == "stop") stop(e)
        warning(sprintf("fit_fn failed for %s: %s", label, conditionMessage(e)))
        NULL
      })

      if (is.null(result_dt)) next

      # Stamp mandatory identifier columns (caller's value always wins)
      result_dt[["fit_source"]] <- fit_source
      if (!is.na(val) && !slice_col %in% names(result_dt)) {
        result_dt[[slice_col]] <- val
      }

      results[[length(results) + 1L]] <- result_dt
    }

    if (length(results) == 0L) {
      data.table::data.table()
    } else {
      data.table::rbindlist(results, fill = TRUE)
    }
  }
}


#' Apply summary_fn across cells and persist into summary_assay
#'
#' @param se SummarizedExperiment
#' @param fit_dt flat data.table with "row", "column", and all fit columns
#' @param summary_fn function(data.table) -> named list
#' @param summary_assay assay name for summary output
#' @param merge merge strategy
#' @param fit_source fit_source label
#' @param on_error error handling strategy
#'
#' @return updated SummarizedExperiment
#'
#' @keywords internal
.apply_summary_fn <- function(se, fit_dt, summary_fn, summary_assay,
                              merge, fit_source, on_error) {
  cell_list <- split(fit_dt, by = c("row", "column"), sorted = FALSE)

  summary_rows <- lapply(cell_list, function(cell_fit) {
    r <- cell_fit[["row"]][1L]
    cc <- cell_fit[["column"]][1L]

    result <- tryCatch({
      res <- summary_fn(cell_fit)
      data.table::as.data.table(as.list(res))
    }, error = function(e) {
      if (on_error == "stop") stop(e)
      warning(sprintf(
        "summary_fn failed for row=%s col=%s: %s", r, cc, conditionMessage(e)
      ))
      NULL
    })

    if (!is.null(result)) {
      result[["row"]] <- r
      result[["column"]] <- cc
      result[["fit_source"]] <- fit_source
    }
    result
  })

  summary_rows <- Filter(Negate(is.null), summary_rows)
  if (length(summary_rows) > 0L) {
    summary_dt <- data.table::rbindlist(summary_rows, fill = TRUE)
    se <- persist_fit_assay(se, summary_dt, merge, summary_assay, "row", "column", "fit_source")
  }
  se
}


#' Persist a data.table into a BumpyMatrix assay
#'
#' Performs an idempotent upsert (merge mode) keyed by \code{upsert_key_cols},
#' or a full overwrite (replace mode).
#'
#' @param se SummarizedExperiment
#' @param new_dt data.table including "row" and "column" index fields
#' @param merge \code{"merge"} or \code{"replace"}
#' @param assay_name assay name to write into
#' @param row row index column name in new_dt
#' @param col column index column name in new_dt
#' @param upsert_key_cols character vector of column(s) forming the upsert key
#'
#' @return updated SummarizedExperiment
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' dt <- data.table::data.table(
#'   row = rownames(se)[1], column = colnames(se)[1],
#'   fit_source = "demo", x_mean = 0.5
#' )
#' se <- persist_fit_assay(se, dt, "merge", "custom_mean", "row", "column", "fit_source")
#' SummarizedExperiment::assayNames(se)
#'
#' @keywords metrics
#' @export
persist_fit_assay <- function(se, new_dt, merge, assay_name, row, col,
                              upsert_key_cols) {
  has_assay <- assay_name %in% SummarizedExperiment::assayNames(se)

  if (merge == "merge" && has_assay) {
    existing_df <- BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se, assay_name),
      row.field = row, column.field = col
    )
    existing_dt <- data.table::as.data.table(existing_df)

    # BumpyMatrix may store character columns as factors; coerce key columns
    # to character so rbindlist does not silently coerce new character values.
    for (kc in upsert_key_cols) {
      if (kc %in% names(existing_dt) && is.factor(existing_dt[[kc]])) {
        existing_dt[, (kc) := as.character(get(kc))]
      }
    }

    # Ensure all key columns exist in existing data (fill with NA if absent)
    for (kc in upsert_key_cols) {
      if (!kc %in% names(existing_dt)) existing_dt[[kc]] <- NA_character_
    }

    key_cols_present <- intersect(upsert_key_cols, names(new_dt))
    # data.table anti-join: keep existing rows whose key is not in new_dt
    existing_pruned <- existing_dt[!new_dt, on = key_cols_present]

    merged_dt <- data.table::rbindlist(list(existing_pruned, new_dt), fill = TRUE)
  } else {
    merged_dt <- new_dt
  }

  data_cols <- names(merged_dt)[!names(merged_dt) %in% c(row, col)]

  # A name absent from the SE would be dropped by the factor levels below, so
  # reject it here rather than silently losing the rows.
  unknown <- c(setdiff(unique(merged_dt[[row]]), rownames(se)),
               setdiff(unique(merged_dt[[col]]), colnames(se)))
  if (length(unknown) > 0L) {
    stop(sprintf("new_dt references row/column names absent from the SE: %s",
                 toString(unknown)))
  }

  # Build over the full SE dimensions: splitAsBumpyMatrix takes them from the
  # factor levels, so a cell with no results comes out empty rather than absent
  # and assay<- does not trip on "rownames not identical". Reshaping afterwards
  # cannot do this — `[` selects existing entries, it does not create them.
  mx <- BumpyMatrix::splitAsBumpyMatrix(
    merged_dt[, data_cols, with = FALSE],
    row = factor(merged_dt[[row]], levels = rownames(se)),
    column = factor(merged_dt[[col]], levels = colnames(se))
  )
  SummarizedExperiment::assay(se, assay_name) <- mx
  se
}


# Constant fit: p_value >= pcutoff or too few points — matches set_constant_fit_params()
# in gDRutils/fit_curves.R (ec50=0, h=0.0001, r2=0, x_mean=x_inf=x_0=mean(x))
#' @keywords internal
.constant_fit_result <- function(norm_type, x, conc, x_std, x_0, range_conc,
                                 p_value = NA_real_, rss = NA_real_) {
  mn <- mean(x, na.rm = TRUE)
  N_conc <- length(unique(conc))
  maxlog10Conc <- if (length(conc) > 0L) log10(max(conc)) else NA_real_
  conc_ord <- order(conc)
  x_ordered <- x[conc_ord]
  n_x <- length(x_ordered)
  x_max_val <- if (n_x >= 2L) {
    min(x_ordered[(n_x - 1L):n_x])
  } else if (n_x == 1L) {
    x_ordered
  } else {
    NA_real_
  }
  # xc50 sign follows the mean normalized value, as .set_mean_params() does in
  # gDRutils/fit_curves.R: above the threshold the curve never reaches 50% (Inf), at or
  # below it is already past 50% at the lowest dose (-Inf). The threshold is per
  # normalization type — see fit_config in fit_profiles.json
  xc50_val <- .estimate_xc50_fallback(mn, .xc50_threshold_for(norm_type))
  list(
    normalization_type = norm_type,
    x_mean = mn,
    x_AOC = 1 - mn,
    x_AOC_range = 1 - mn,
    x_max = x_max_val,
    x_sd_avg = mean(x_std, na.rm = TRUE),
    N_conc = N_conc,
    maxlog10Concentration = maxlog10Conc,
    ec50 = 0,
    xc50 = xc50_val,
    h = 0.0001,
    r2 = 0,
    rss = rss,
    p_value = p_value,
    x_0 = if (is.na(x_0)) mn else x_0,
    x_inf = mn,
    fit_type = "DRCConstantFitResult"
  )
}


# Too few concentrations to fit a curve: the summary of the measured values, without
# curve parameters. xc50 is +/-Inf only when every response is on the same side of the
# threshold, as .set_too_few_fit_params() does in fit_curves.R.
#' @keywords internal
.too_few_fit_result <- function(norm_type, x, conc, x_std) {
  out <- .empty_fit_result(norm_type)
  out$fit_type <- "DRCTooFewPointsToFit"
  if (length(x) == 0L) return(out)
  mn <- mean(x)
  x_ordered <- x[order(conc)]
  n_x <- length(x_ordered)
  out$x_mean <- mn
  out$x_AOC <- 1 - mn
  out$x_max <- if (n_x >= 2L) min(x_ordered[(n_x - 1L):n_x]) else x_ordered
  out$x_sd_avg <- mean(x_std, na.rm = TRUE)
  out$N_conc <- length(unique(conc))
  out$maxlog10Concentration <- log10(max(conc))
  out$xc50 <- .estimate_xc50_fallback(x, .xc50_threshold_for(norm_type))
  out
}


#' @keywords internal
.empty_fit_result <- function(norm_type) {
  list(
    normalization_type = norm_type,
    x_mean = NA_real_,
    x_AOC = NA_real_,
    x_AOC_range = NA_real_,
    x_max = NA_real_,
    x_sd_avg = NA_real_,
    N_conc = 0L,
    maxlog10Concentration = NA_real_,
    ec50 = NA_real_,
    xc50 = NA_real_,
    h = NA_real_,
    r2 = NA_real_,
    rss = NA_real_,
    p_value = NA_real_,
    x_0 = NA_real_,
    x_inf = NA_real_,
    fit_type = "DRCInvalidFitResult"
  )
}


#' @keywords internal
.estimate_xc50_fallback <- function(x, threshold = 0.5) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0L) {
    NA_real_
  } else if (all(x_clean > threshold)) {
    Inf
  } else if (all(x_clean <= threshold)) {
    -Inf
  } else {
    NA_real_
  }
}
