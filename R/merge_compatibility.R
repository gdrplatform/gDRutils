#' Compare the SummarizedExperiments that \code{merge_SE()} would combine
#'
#' \code{merge_SE()} stacks the assays of its inputs. This reports what the inputs
#' have in common and where they differ, before anything is stacked, so a merge
#' that would mix things up can be seen first.
#'
#' @details
#' The report has one table per topic, each with a \code{data_source} column naming
#' the input (the name in \code{SElist}):
#' \describe{
#'   \item{\code{identifiers}}{the column each input uses for every identifier, and whether
#'     it differs from the first input's. \code{merge_SE()} renames the columns of the
#'     other inputs to the first input's.}
#'   \item{\code{assays}}{which assays each input has.}
#'   \item{\code{normalization_types}}{the normalization types in each \code{Metrics}
#'     assay.}
#'   \item{\code{fit_recipes}}{the fits recorded in each input's provenance: assay,
#'     \code{fit_source}, \code{fit_id}, function and solver.}
#'   \item{\code{fit_parameters}, \code{versions}}{the fit parameters and the package
#'     versions each input was processed with.}
#'   \item{\code{issues}}{what the others add up to, with a \code{severity}.}
#' }
#'
#' A \code{"contradiction"} is the same \code{fit_source} in the same assay produced
#' by different methods in different inputs: another fit function or another solver.
#' The merged rows would carry one label for two methods, so \code{merge_SE()} stops
#' unless \code{allow_mismatch = TRUE}. The same method with other settings
#' (\code{fit_args}, e.g. the concentration range, which depends on the data) or in
#' another package version is reported as a \code{"difference"}, since data processed
#' with different settings or in different releases is often merged on purpose.
#' Everything else (other identifiers, missing assays, other normalization types,
#' other fit parameters or versions, an input without provenance) is a
#' \code{"difference"} too and does not stop the merge.
#'
#' @param SElist named list of \code{SummarizedExperiment}s, as for \code{merge_SE()}.
#'
#' @return object of class \code{merge_compatibility}: a list of \code{data.table}s
#'   (see Details).
#'
#' @examples
#' se <- get_synthetic_data("finalMAE_small.qs2")[[1]]
#' merge_compatibility(list(a = se, b = se))
#'
#' @seealso \code{\link{merge_SE}}
#' @keywords SE_operators
#' @export
merge_compatibility <- function(SElist) {
  checkmate::assert_list(SElist, types = "SummarizedExperiment", min.len = 1L)
  nms <- names(SElist)
  if (is.null(nms)) nms <- rep("", length(SElist))
  nms[nms == ""] <- as.character(which(nms == ""))
  names(SElist) <- nms
  first <- nms[[1L]]

  paste_value <- function(x) paste(format(x), collapse = " | ")

  ids <- data.table::rbindlist(lapply(nms, function(n) {
    x <- get_SE_identifiers(SElist[[n]])
    data.table::data.table(data_source = n, identifier = names(x),
                           value = vapply(x, paste_value, character(1)))
  }))
  target <- ids[ids[["data_source"]] == first, ]
  ids$differs <- is.na(match(paste(ids$identifier, ids$value), paste(target$identifier, target$value)))

  assay_names <- unique(unlist(lapply(SElist, SummarizedExperiment::assayNames)))
  assays <- data.table::CJ(data_source = nms, assay = assay_names, sorted = FALSE)
  assays$present <- mapply(function(n, a) a %in% SummarizedExperiment::assayNames(SElist[[n]]),
                           assays$data_source, assays$assay)

  norm <- data.table::rbindlist(lapply(nms, function(n) {
    se <- SElist[[n]]
    if (!"Metrics" %in% SummarizedExperiment::assayNames(se)) return(NULL)
    metrics <- .assay_dt(se, "Metrics")
    if (!"normalization_type" %in% names(metrics)) return(NULL)
    data.table::data.table(
      data_source = n,
      normalization_type = sort(unique(as.character(metrics[["normalization_type"]])))
    )
  }))
  if (NROW(norm) == 0L) norm <- data.table::data.table(data_source = character(0),
                                                       normalization_type = character(0))

  recipes <- data.table::rbindlist(lapply(nms, function(n) {
    prov <- .normalize_provenance(get_SE_provenance(SElist[[n]]))
    rows <- list()
    for (a in names(prov$fits)) {
      for (id in names(prov$fits[[a]])) {
        r <- prov$fits[[a]][[id]]
        if (!is.null(r$summary_of)) next
        rows[[length(rows) + 1L]] <- data.table::data.table(
          data_source = n, assay = a, fit_source = r$fit_source, fit_id = id,
          called_by = if (is.null(r$called_by)) NA_character_ else r$called_by,
          data_type = if (is.null(r$data_type)) NA_character_ else r$data_type,
          fn = .format_fn(r$fn),
          solver = if (is.null(r$solver)) NA_character_ else r$solver,
          signature = .recipe_signature(r),
          settings = digest::digest(r[["fit_args"]], algo = "xxhash64")
        )
      }
    }
    data.table::rbindlist(rows)
  }), fill = TRUE)
  if (NROW(recipes) == 0L) {
    recipes <- data.table::data.table(
      data_source = character(0), assay = character(0), fit_source = character(0),
      fit_id = character(0), called_by = character(0), data_type = character(0),
      fn = character(0), solver = character(0), signature = character(0), settings = character(0)
    )
  }

  as_table <- function(get, value_name) {
    out <- data.table::rbindlist(lapply(nms, function(n) {
      x <- tryCatch(get(SElist[[n]]), error = function(e) NULL)
      if (length(x) == 0L) return(NULL)
      data.table::data.table(data_source = n, name = names(x),
                             value = vapply(x, paste_value, character(1)))
    }))
    if (NROW(out) == 0L) out <- data.table::data.table(data_source = character(0), name = character(0),
                                                       value = character(0))
    data.table::setnames(out, "name", value_name)
    out
  }
  fit_parameters <- as_table(get_SE_fit_parameters, "parameter")
  versions <- as_table(function(se) get_SE_processing_metadata(se)$versions, "package")

  issues <- .compatibility_issues(
    nms, ids, assays, norm, recipes, fit_parameters, versions,
    recorded = vapply(SElist, function(se) !is.null(get_SE_provenance(se)), logical(1))
  )

  structure(
    list(identifiers = ids, assays = assays, normalization_types = norm,
         fit_recipes = recipes, fit_parameters = fit_parameters, versions = versions,
         issues = issues),
    class = "merge_compatibility"
  )
}


#' @param x a \code{merge_compatibility} object.
#' @param ... unused.
#' @rdname merge_compatibility
#' @export
print.merge_compatibility <- function(x, ...) {
  n_sources <- data.table::uniqueN(x$identifiers$data_source)
  cat(sprintf("Merge compatibility of %d input%s\n", n_sources, if (n_sources == 1L) "" else "s"))
  if (NROW(x$issues) == 0L) {
    cat("  no differences\n")
  } else {
    for (i in seq_len(NROW(x$issues))) {
      cat(sprintf("  [%s] %s: %s\n", x$issues$severity[[i]], x$issues$topic[[i]], x$issues$detail[[i]]))
    }
  }
  invisible(x)
}


# The method of a fit as a signature: the function and any function passed to it (a
# solver), not its settings, the package version it ran in, the label of the function
# or the data it ran on.
#' @keywords internal
.recipe_signature <- function(recipe) {
  keep <- function(rec) {
    if (is.null(rec)) NULL else rec[setdiff(names(rec), c(
      "name", "info", "env", "not_stored", "version", "git_sha", "hash", "dependencies"
    ))]
  }
  digest::digest(list(
    called_by = recipe$called_by,
    fn = keep(recipe$fn),
    summary_fn = keep(recipe$summary_fn),
    fn_args = lapply(recipe[["fn_args"]], keep)
  ), algo = "xxhash64")
}


#' @keywords internal
.compatibility_issues <- function(nms, ids, assays, norm, recipes, fit_parameters, versions,
                                  recorded) {
  rows <- list()
  add <- function(severity, topic, detail) {
    rows[[length(rows) + 1L]] <<- data.table::data.table(
      severity = severity, topic = topic, detail = detail
    )
  }
  first <- nms[[1L]]

  # fits: same fit_source and assay, different fit in different inputs. Only fits of the
  # same kind are compared: the Metrics of a combination experiment are made by another
  # step than those of a single-agent experiment, under the same label, on purpose.
  if (NROW(recipes) > 0L) {
    by_label <- split(recipes, paste(recipes$assay, recipes$fit_source, sep = "\r"))
    for (g in by_label) {
      kinds <- unique(paste(g$called_by, g$data_type))
      if (length(kinds) > 1L && data.table::uniqueN(g$data_source) > 1L) {
        add("difference", "fit_recipes", sprintf(
          "fit_source '%s' in assay '%s' comes from different steps: %s",
          g$fit_source[[1L]], g$assay[[1L]],
          paste(unique(sprintf("%s (%s)", g$data_source, g$fn)), collapse = ", ")
        ))
      }
    }
    for (g in split(recipes, paste(recipes$assay, recipes$fit_source, recipes$called_by,
                                   recipes$data_type, sep = "\r"))) {
      per_input <- vapply(split(g$signature, g$data_source),
                          function(s) paste(sort(unique(s)), collapse = ","), character(1))
      if (length(per_input) > 1L && length(unique(per_input)) > 1L) {
        add("contradiction", "fit_recipes", sprintf(
          "fit_source '%s' in assay '%s' was fitted differently: %s",
          g$fit_source[[1L]], g$assay[[1L]],
          paste(sprintf("%s (%s)", unique(g[, c("data_source", "fn")])$data_source,
                        unique(g[, c("data_source", "fn")])$fn), collapse = ", ")
        ))
      } else if (length(per_input) > 1L) {
        same_settings <- length(unique(vapply(
          split(g$settings, g$data_source), function(s) paste(sort(unique(s)), collapse = ","), character(1)
        ))) == 1L
        if (!same_settings) {
          add("difference", "fit_recipes", sprintf(
            "fit_source '%s' in assay '%s' was fitted by the same method with other settings (fit_id %s)",
            g$fit_source[[1L]], g$assay[[1L]], toString(unique(g$fit_id))
          ))
        } else if (length(unique(g$fit_id)) > 1L) {
          add("difference", "fit_recipes", sprintf(
            "fit_source '%s' in assay '%s' was fitted by the same method in other package versions (fit_id %s)",
            g$fit_source[[1L]], g$assay[[1L]], toString(unique(g$fit_id))
          ))
        }
      }
    }
  }
  for (n in nms[!recorded]) {
    add("difference", "provenance", sprintf("'%s' has no provenance: its fits cannot be compared", n))
  }

  differing_keys <- unique(ids$identifier[ids$differs])
  for (k in differing_keys) {
    d <- ids[ids$identifier == k & ids$differs, ]
    add("difference", "identifiers", sprintf(
      "identifier '%s' differs from '%s' in %s: %s", k, first, toString(sQuote(d$data_source, FALSE)),
      paste(sprintf("%s = %s", d$data_source, d$value), collapse = "; ")
    ))
  }
  for (a in unique(assays$assay[!assays$present])) {
    add("difference", "assays", sprintf(
      "assay '%s' is missing in %s", a, toString(sQuote(assays$data_source[assays$assay == a & !assays$present], FALSE))
    ))
  }
  if (data.table::uniqueN(norm$normalization_type) > 0L) {
    per_input <- vapply(split(norm$normalization_type, norm$data_source),
                        function(s) paste(sort(s), collapse = ","), character(1))
    if (length(unique(per_input)) > 1L) {
      add("difference", "normalization_types", sprintf(
        "normalization types differ: %s", paste(sprintf("%s = %s", names(per_input), per_input), collapse = "; ")
      ))
    }
  }
  for (tbl in list(list(fit_parameters, "parameter", "fit_parameters"), list(versions, "package", "versions"))) {
    dt <- tbl[[1L]]
    if (NROW(dt) == 0L) next
    for (p in unique(dt[[tbl[[2L]]]])) {
      v <- dt[dt[[tbl[[2L]]]] == p, ]
      if (length(unique(v$value)) > 1L) {
        add("difference", tbl[[3L]], sprintf(
          "%s '%s' differs: %s", sub("s$", "", tbl[[3L]]), p,
          paste(sprintf("%s = %s", v$data_source, v$value), collapse = "; ")
        ))
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.table::data.table(severity = character(0), topic = character(0), detail = character(0)))
  }
  data.table::rbindlist(rows)
}


# The report as plain data, for metadata(se)$merge_report.
#' @keywords internal
.merge_report <- function(report, target_identifiers, allow_mismatch) {
  # plain data.frames, which serialize and compare as data, unlike data.tables
  tables <- lapply(unclass(report), function(x) data.table::setDF(data.table::copy(x)))
  c(tables, list(
    target_identifiers = target_identifiers,
    allow_mismatch = allow_mismatch
  ))
}


# Rename the identifier columns of one input to the target's, so inputs that name the
# same thing differently stack into one column. Only identifiers that name one column
# in both are renamed.
#' @keywords internal
.to_target_identifiers <- function(dt, own, target) {
  for (k in setdiff(intersect(names(own), names(target)), "data_source")) {
    from <- own[[k]]
    to <- target[[k]]
    if (length(from) != 1L || length(to) != 1L || identical(from, to)) next
    if (from %in% names(dt) && !to %in% names(dt)) data.table::setnames(dt, from, to)
  }
  dt
}
