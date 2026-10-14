####
# Fit provenance
####
#
# Every fit written by apply_fit() / apply_fits() leaves its recipe in
# metadata(se)$provenance, so a stored object says which code produced each result
# and the fit can be re-run later. The entry holds plain lists, strings and numbers
# only: it serializes to qs2 and DSDB as is and can be exported to JSON.
#
#   provenance
#     schema_version  integer, bumped when the layout changes
#     fits            <assay> -> <fit_source> -> recipe
#     steps           processing steps before the fit (filled by gDRcore)
#     sources         <data_source> -> provenance of each SE combined by merge_SE()
#
# The fit_source column of an assay is the key into fits[[assay]]. A recipe holds a
# record of the fit function (.describe_fn()), the arguments passed in fit_args, the
# call settings and a short run block.
#
# Functions from installed packages are stored as a reference (package, name,
# version, hash): their code ships with the package, and copying it into every
# object would only make the object bigger. Anything else (a function written in the
# session, an inline function, a closure) is stored as source code together with
# the values and helper functions it uses, since it cannot be recovered otherwise.

.provenance_schema_version <- 1L

# Captured values and fit_args are stored only when they are plain data no larger
# than this; anything bigger is listed by name instead, so a fit never grows the SE
# by more than a few kilobytes.
.provenance_max_value_bytes <- 10240

# Packages that ship with R: their version is the R version, already in the run block.
.provenance_base_packages <- c(
  "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods",
  "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"
)


####
# Public API
####

#' Describe a fit function
#'
#' Attach a short human-readable description to a fit function. It is copied
#' into the recipe that \code{\link{apply_fit}} and \code{\link{apply_fits}}
#' record in \code{metadata(se)$provenance}, so a stored object says not only
#' which function produced a fit but also what that function does.
#'
#' @param fn the fit function to describe.
#' @param model string; name of the model, e.g.
#'   \code{"3-parameter log-logistic"}.
#' @param description string; what the function computes, in a sentence or
#'   two.
#' @param reference string; citation of the method, ideally with a DOI.
#' @param env optional function without arguments returning a named list of
#'   strings that describe the external runtime the fit depends on, e.g. the
#'   version of a Python package called through reticulate. It is called each
#'   time a fit is recorded.
#'
#' @return \code{fn} with the description attached.
#'
#' @examples
#' mean_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
#' mean_fn <- set_fit_fn_info(
#'   mean_fn,
#'   model = "mean response",
#'   description = "Mean of the averaged response over all concentrations."
#' )
#'
#' @seealso \code{\link{get_fit_provenance}}
#' @keywords provenance
#' @export
#'
set_fit_fn_info <- function(fn,
                            model = NULL,
                            description = NULL,
                            reference = NULL,
                            env = NULL) {
  checkmate::assert_function(fn)
  checkmate::assert_string(model, null.ok = TRUE)
  checkmate::assert_string(description, null.ok = TRUE)
  checkmate::assert_string(reference, null.ok = TRUE)
  checkmate::assert_function(env, nargs = 0L, null.ok = TRUE)

  info <- .drop_null(list(model = model, description = description, reference = reference))
  attr(fn, "fit_fn_info") <- if (length(info) > 0L) info else NULL
  attr(fn, "fit_fn_env") <- env
  fn
}


#' Get the provenance of fits stored in a SummarizedExperiment
#'
#' \code{\link{apply_fit}} and \code{\link{apply_fits}} record a recipe for
#' every fit in \code{metadata(se)$provenance}. These functions read it back.
#'
#' @details
#' A recipe contains:
#' \describe{
#'   \item{\code{fn}, \code{summary_fn}}{record of the function: \code{name},
#'     \code{package}, \code{version}, \code{git_sha} and \code{hash} for a
#'     function from an installed package; \code{source}, \code{captured}
#'     (values the function uses from its environment), \code{helpers}
#'     (source of helper functions it calls) and \code{dependencies} (versions
#'     of the packages it calls) for any other function; \code{info} set with
#'     \code{\link{set_fit_fn_info}}; \code{env} for declared runtime
#'     dependencies.}
#'   \item{\code{fit_args}}{arguments passed to the fit function.}
#'   \item{\code{called_by}, \code{data_type}, \code{input_assay},
#'     \code{output_assay}, \code{summary_assay}, \code{slicing_cols},
#'     \code{slicing_values}, \code{merge}}{settings of the call.}
#'   \item{\code{run}}{R and gDRutils versions, hash of the input assay and
#'     the number of rows written.}
#' }
#' Names listed under \code{not_stored} or \code{args_not_stored} were too
#' large or not plain data, so the fit cannot be fully reproduced from the
#' object alone.
#'
#' For a summary assay, \code{get_fit_provenance()} returns the recipe of the
#' fit that wrote it. For an SE built by \code{\link{merge_SE}}, the recipes of
#' the merged objects are kept per \code{data_source}.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @param assay string; name of the assay holding the fit results.
#' @param fit_source string; value of the \code{fit_source} column.
#' @param data_source string; for an SE built by \code{\link{merge_SE}}, the
#'   name of the merged object the rows come from (the \code{data_source}
#'   column). \code{NULL} reads the fits done on \code{se} itself.
#'
#' @return \code{get_SE_provenance()}: the whole provenance entry, or
#'   \code{NULL} for objects without one. \code{get_fit_provenance()}: one
#'   recipe. \code{summarize_fit_provenance()}: a \code{data.table} with one
#'   row per recorded fit.
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' se <- apply_fit(
#'   se, fit_drug_response_metrics, "single-agent",
#'   output_assay = "custom_hill", fit_source = "hill_loose",
#'   fit_args = list(pcutoff = 0.1)
#' )
#' summarize_fit_provenance(se)
#' get_fit_provenance(se, "custom_hill", "hill_loose")$fit_args
#'
#' @name fit_provenance
#' @keywords provenance
NULL


#' @rdname fit_provenance
#' @export
#'
get_SE_provenance <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  S4Vectors::metadata(se)[["provenance"]]
}


#' @rdname fit_provenance
#' @export
#'
get_fit_provenance <- function(se, assay, fit_source, data_source = NULL) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay)
  checkmate::assert_string(fit_source)
  checkmate::assert_string(data_source, null.ok = TRUE)

  prov <- get_SE_provenance(se)
  if (!is.null(data_source)) {
    prov <- prov$sources[[data_source]]
    if (is.null(prov)) {
      stop(sprintf("no provenance recorded for data_source '%s'", data_source))
    }
  }
  .find_recipe(prov, assay, fit_source)
}


#' @rdname fit_provenance
#' @export
#'
summarize_fit_provenance <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  .summarize_provenance(get_SE_provenance(se), data_source = NA_character_)
}


#' Re-run a recorded fit
#'
#' Rebuild the fit function from its recipe in \code{metadata(se)$provenance}
#' and fit the data again. The rows of \code{fit_source} in the output assay
#' (and in the summary assay, if any) are replaced; other rows are kept.
#'
#' @details
#' A function from an installed package is taken from that package. A warning
#' is raised when the installed version or the code differs from the recorded
#' one, since results may then differ.
#'
#' A function stored as source code is evaluated only with
#' \code{allow_source = TRUE}: the object may come from someone else, and
#' evaluating its code runs that code in your session.
#'
#' @inheritParams get_fit_provenance
#' @param allow_source logical; evaluate fit functions stored as source code.
#'
#' @return \code{se} with the fit re-run.
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' se <- apply_fit(
#'   se, fit_drug_response_metrics, "single-agent",
#'   output_assay = "custom_hill", fit_source = "hill"
#' )
#' se <- refit_from_provenance(se, "custom_hill", "hill")
#'
#' @seealso \code{\link{check_reproducibility}}
#' @keywords provenance
#' @export
#'
refit_from_provenance <- function(se, assay, fit_source, allow_source = FALSE) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_flag(allow_source)
  recipe <- get_fit_provenance(se, assay, fit_source)

  fn <- .restore_fn(recipe$fn, allow_source)
  summary_fn <- if (!is.null(recipe$summary_fn)) .restore_fn(recipe$summary_fn, allow_source)
  if (length(recipe$args_not_stored) > 0L) {
    .provenance_warning(sprintf(
      "fit_args %s were not stored, so the fit runs with their default values",
      toString(sQuote(recipe$args_not_stored, FALSE))
    ))
  }
  fit_args <- recipe[["fit_args"]]
  if (is.null(fit_args)) fit_args <- list()
  # functions passed as arguments (e.g. a score function) are stored as records
  fn_args <- lapply(recipe[["fn_args"]], .restore_fn, allow_source = allow_source)
  fit_args <- c(fit_args, fn_args)
  output_assay <- recipe$output_assay

  # Re-run with merge = "merge" whatever the original call used: "replace" would
  # also remove the rows of every other fit_source in the assay.
  se <- if (identical(recipe$called_by, "se_function")) {
    # a function of the whole SE (e.g. a combination step) records its own recipe again
    do.call(fn, c(list(se), fit_args))
  } else if (identical(recipe$called_by, "apply_fits")) {
    apply_fits(
      se,
      fit_fns = stats::setNames(list(fn), output_assay),
      data_type = recipe$data_type,
      slicing_cols = recipe$slicing_cols,
      slicing_values = recipe$slicing_values,
      input_assay = recipe$input_assay,
      merge = "merge",
      fit_args = stats::setNames(list(fit_args), output_assay),
      fit_source = fit_source
    )
  } else {
    apply_fit(
      se,
      fit_fn = fn,
      data_type = recipe$data_type,
      slicing_cols = recipe$slicing_cols,
      slicing_values = recipe$slicing_values,
      input_assay = recipe$input_assay,
      output_assay = output_assay,
      summary_fn = summary_fn,
      summary_assay = recipe$summary_assay,
      merge = "merge",
      fit_args = fit_args,
      fit_source = fit_source
    )
  }

  # The new recipe describes the rebuilt function, which no longer carries the
  # name it was called by or its runtime hook; keep those from the original.
  new_recipe <- S4Vectors::metadata(se)$provenance$fits[[output_assay]][[fit_source]]
  if (!is.null(new_recipe)) {
    for (slot in c("fn", "summary_fn")) {
      for (field in c("name", "env", "info")) {
        if (!is.null(recipe[[slot]][[field]])) new_recipe[[slot]][[field]] <- recipe[[slot]][[field]]
      }
    }
    S4Vectors::metadata(se)$provenance$fits[[output_assay]][[fit_source]] <- new_recipe
  }
  se
}


#' Check that recorded fits can be reproduced
#'
#' Re-run each recorded fit with \code{\link{refit_from_provenance}} and
#' compare the new rows with the stored ones.
#'
#' @inheritParams refit_from_provenance
#' @param assay string or \code{NULL}; check only fits written to this assay.
#' @param fit_source string or \code{NULL}; check only this fit source.
#' @param tolerance numeric; relative tolerance for numeric columns, as in
#'   \code{\link[base]{all.equal}}.
#'
#' @return \code{data.table} with one row per checked fit: \code{assay},
#'   \code{fit_source}, \code{reproducible} (\code{NA} when the fit could not
#'   be re-run), \code{n_rows} and \code{detail}.
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' se <- apply_fit(
#'   se, fit_drug_response_metrics, "single-agent",
#'   output_assay = "custom_hill", fit_source = "hill"
#' )
#' check_reproducibility(se)
#'
#' @keywords provenance
#' @export
#'
check_reproducibility <- function(se,
                                  assay = NULL,
                                  fit_source = NULL,
                                  tolerance = sqrt(.Machine$double.eps),
                                  allow_source = FALSE) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay, null.ok = TRUE)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_number(tolerance, lower = 0)
  checkmate::assert_flag(allow_source)

  fits <- get_SE_provenance(se)$fits
  targets <- list()
  for (a in names(fits)) {
    for (s in names(fits[[a]])) {
      # summary assays point at the fit that wrote them and are checked with it
      if (!is.null(fits[[a]][[s]]$summary_of)) next
      if (!is.null(assay) && a != assay) next
      if (!is.null(fit_source) && s != fit_source) next
      targets[[length(targets) + 1L]] <- c(assay = a, fit_source = s)
    }
  }
  if (length(targets) == 0L) {
    return(data.table::data.table(
      assay = character(0), fit_source = character(0), reproducible = logical(0),
      n_rows = integer(0), detail = character(0)
    ))
  }

  data.table::rbindlist(lapply(targets, function(t) {
    .check_one_fit(se, t[["assay"]], t[["fit_source"]], tolerance, allow_source)
  }))
}


#' Record provenance from code that does not use apply_fit()
#'
#' \code{\link{apply_fit}} and \code{\link{apply_fits}} record their recipes
#' themselves. Code that writes results another way uses these functions to
#' add to the same \code{metadata(se)$provenance} entry.
#'
#' @details
#' \code{record_fit_provenance()} is for a function of the whole SE that
#' writes an assay, e.g. a combination scoring step. The recipe keeps the
#' function and its arguments, so \code{\link{refit_from_provenance}} can call
#' it again as \code{fn(se, ...)}. Functions passed among the arguments (e.g.
#' a score function) are stored as function records.
#'
#' \code{record_processing_step()} is for a step whose result cannot be
#' recomputed from the object it returns, e.g. growth rates computed from an
#' assay that the output no longer carries. It stores the functions and
#' settings of the step under \code{provenance$steps[[step]]}.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @param assay string; assay the function wrote.
#' @param fit_source string; \code{fit_source} of the written rows. Assays
#'   without a \code{fit_source} column use it only as a key.
#' @param fn function of the whole SE that wrote \code{assay}, called as
#'   \code{fn(se, ...)}.
#' @param fit_args named list of the arguments \code{fn} was called with,
#'   other than \code{se}.
#' @param input_assay string or \code{NULL}; assay the results were computed
#'   from, hashed to tell later whether a re-run used the same data.
#' @param step string; name of the processing step.
#' @param settings named list of plain-data settings of the step.
#' @param fns named list of the functions the step used, e.g.
#'   \code{list(rate_fn = rate_fn)}.
#'
#' @return \code{se} with the provenance entry updated.
#'
#' @examples
#' mae <- get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' se <- record_processing_step(
#'   se, "growth_rates",
#'   settings = list(periods = list(early = c(0, 48))),
#'   fns = list(rate_fn = function(dt) 0)
#' )
#' names(get_SE_provenance(se)$steps)
#'
#' @name record_provenance
#' @keywords provenance
NULL


#' @rdname record_provenance
#' @export
#'
record_fit_provenance <- function(se, assay, fit_source, fn, fit_args = list(), input_assay = NULL) {
  fn_name <- .call_name(substitute(fn))
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(assay)
  checkmate::assert_string(fit_source)
  checkmate::assert_function(fn)
  checkmate::assert_list(fit_args, names = "unique")
  checkmate::assert_string(input_assay, null.ok = TRUE)

  recipe <- .fit_recipe(
    called_by = "se_function",
    fn_record = .describe_fn(fn, fn_name),
    fit_args = fit_args,
    summary_record = NULL,
    data_type = NULL,
    input_assay = input_assay,
    output_assay = assay,
    summary_assay = NULL,
    slicing_cols = NULL,
    slicing_values = NULL,
    merge = NULL,
    input_hash = if (!is.null(input_assay)) .assay_hash(se, input_assay),
    n_rows = NROW(.fit_rows(se, assay, fit_source))
  )
  .record_fit(se, assay, fit_source, recipe)
}


#' @rdname record_provenance
#' @export
#'
record_processing_step <- function(se, step, settings = list(), fns = list()) {
  fn_names <- .list_call_names(substitute(fns))
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(step)
  checkmate::assert_list(settings, names = "unique")
  checkmate::assert_list(fns, types = "function", names = "unique")

  stored <- .split_storable(settings)
  functions <- lapply(names(fns), function(role) .describe_fn(fns[[role]], fn_names[[role]]))
  names(functions) <- names(fns)
  entry <- .drop_null(list(
    functions = if (length(functions) > 0L) functions,
    settings = if (length(stored$stored) > 0L) stored$stored,
    settings_not_stored = if (length(stored$not_stored) > 0L) stored$not_stored,
    run = list(
      R = paste(R.version$major, R.version$minor, sep = "."),
      gDRutils = as.character(utils::packageVersion("gDRutils"))
    )
  ))

  prov <- S4Vectors::metadata(se)[["provenance"]]
  if (is.null(prov)) {
    prov <- list(schema_version = .provenance_schema_version, fits = list(), steps = list())
  }
  prov$steps[[step]] <- entry
  S4Vectors::metadata(se)[["provenance"]] <- prov

  not_stored <- c(entry$settings_not_stored,
                  unlist(lapply(functions, `[[`, "not_stored"), use.names = FALSE))
  if (length(not_stored) > 0L) {
    .provenance_warning(sprintf(
      "step '%s' cannot be fully reproduced from the object: %s not stored",
      step, toString(sQuote(unique(not_stored), FALSE))
    ))
  }
  se
}


####
# Recording
####

# Build the recipe of one fit.
#' @keywords internal
.fit_recipe <- function(called_by,
                        fn_record,
                        fit_args,
                        summary_record,
                        data_type,
                        input_assay,
                        output_assay,
                        summary_assay,
                        slicing_cols,
                        slicing_values,
                        merge,
                        input_hash,
                        n_rows) {
  is_fn <- vapply(fit_args, is.function, logical(1))
  fn_args <- lapply(fit_args[is_fn], .describe_fn)
  stored_args <- .split_storable(fit_args[!is_fn])
  .drop_null(list(
    called_by = called_by,
    fn = fn_record,
    fit_args = if (length(stored_args$stored) > 0L) stored_args$stored,
    fn_args = if (length(fn_args) > 0L) fn_args,
    args_not_stored = if (length(stored_args$not_stored) > 0L) stored_args$not_stored,
    summary_fn = summary_record,
    data_type = data_type,
    input_assay = input_assay,
    output_assay = output_assay,
    summary_assay = summary_assay,
    slicing_cols = slicing_cols,
    slicing_values = slicing_values,
    merge = merge,
    # no timestamp: the same data fitted with the same code gives an identical
    # object, which downstream checks rely on
    run = list(
      R = paste(R.version$major, R.version$minor, sep = "."),
      gDRutils = as.character(utils::packageVersion("gDRutils")),
      input_hash = input_hash,
      n_rows = as.integer(n_rows)
    )
  ))
}


# Store a recipe under fits[[assay]][[fit_source]]. merge = "replace" overwrote the
# whole assay, so the recipes of the other fit sources in it no longer apply.
#' @keywords internal
.record_fit <- function(se, assay, fit_source, recipe, replace = FALSE, summary_assay = NULL) {
  prov <- S4Vectors::metadata(se)[["provenance"]]
  if (is.null(prov)) {
    prov <- list(schema_version = .provenance_schema_version, fits = list(), steps = list())
  }
  if (replace || is.null(prov$fits[[assay]])) prov$fits[[assay]] <- list()
  prov$fits[[assay]][[fit_source]] <- recipe
  if (!is.null(summary_assay)) {
    if (is.null(prov$fits[[summary_assay]])) prov$fits[[summary_assay]] <- list()
    prov$fits[[summary_assay]][[fit_source]] <- list(summary_of = assay)
  }
  S4Vectors::metadata(se)[["provenance"]] <- prov

  not_stored <- c(recipe$fn$not_stored, recipe$summary_fn$not_stored, recipe$args_not_stored,
                  unlist(lapply(recipe[["fn_args"]], `[[`, "not_stored"), use.names = FALSE))
  if (length(not_stored) > 0L) {
    .provenance_warning(sprintf(
      paste0("fit '%s' cannot be fully reproduced from the object: %s not stored because ",
             "they are not plain data or exceed %d bytes; pass small settings through fit_args"),
      fit_source, toString(sQuote(unique(not_stored), FALSE)), .provenance_max_value_bytes
    ))
  }
  se
}


# Provenance of an SE built by merge_SE(): nothing was fitted on the merged object
# itself, so each input keeps its own provenance under its data_source name.
#' @keywords internal
.merge_provenance <- function(SElist) {
  sources <- lapply(SElist, function(se) S4Vectors::metadata(se)[["provenance"]])
  sources <- sources[!vapply(sources, is.null, logical(1))]
  list(
    schema_version = .provenance_schema_version,
    fits = list(),
    steps = list(),
    sources = sources
  )
}


####
# Function records
####

# Record a function so that it can be identified and, if possible, rebuilt.
# `name` is the name the caller used for it, when known.
#' @keywords internal
.describe_fn <- function(fn, name = NULL) {
  code <- .deparse_fn(fn)
  pkg <- .fn_package(fn)
  ref_name <- if (!is.na(pkg)) .find_in_namespace(fn, pkg, name) else NA_character_

  if (is.na(pkg)) name <- .global_name(fn, name)

  record <- list(
    name = if (!is.na(ref_name)) ref_name else name,
    package = if (!is.na(pkg)) pkg,
    version = if (!is.na(pkg)) .package_version(pkg),
    git_sha = if (!is.na(pkg)) .package_git_sha(pkg),
    hash = digest::digest(code, algo = "xxhash64", serialize = FALSE)
  )

  if (is.na(ref_name)) {
    # Not reachable by name in a package namespace (a function written in the
    # session, an inline function or a closure built at run time): keep its code
    # and the state it reads from its environment.
    state <- .capture_state(fn)
    record$source <- code
    record$captured <- if (length(state$values) > 0L) state$values
    record$helpers <- if (length(state$helpers) > 0L) state$helpers
    record$aliases <- if (length(state$aliases) > 0L) state$aliases
    record$not_stored <- if (length(state$not_stored) > 0L) state$not_stored
    deps <- state$packages
  } else {
    # the package's own imports pin the code it calls (e.g. drc for a Hill fit)
    deps <- c(pkg, .package_imports(pkg))
  }
  record$dependencies <- .package_versions(deps)
  record$info <- attr(fn, "fit_fn_info")
  record$env <- .fit_fn_env(fn)
  .drop_null(record)
}


# Code of a function without source references or attributes, so the text (and its
# hash) does not depend on how the function was created.
#' @keywords internal
.deparse_fn <- function(fn) {
  if (is.primitive(fn)) return(paste(deparse(fn), collapse = "\n"))
  fn <- utils::removeSource(fn)
  attributes(fn) <- NULL
  paste(deparse(fn, width.cutoff = 80L), collapse = "\n")
}


# Name of the package namespace a function was defined in, or NA.
#' @keywords internal
.fn_package <- function(fn) {
  if (is.primitive(fn)) return("base")
  top <- topenv(environment(fn))
  if (isNamespace(top)) getNamespaceName(top)[[1L]] else NA_character_
}


# Name under which `fn` can be retrieved from the namespace of `pkg`, or NA when it
# cannot (e.g. a closure built inside a package function).
#' @keywords internal
.find_in_namespace <- function(fn, pkg, name = NULL) {
  ns <- asNamespace(pkg)
  if (!is.null(name) && exists(name, envir = ns, inherits = FALSE) &&
      identical(get(name, envir = ns), fn)) {
    return(name)
  }
  # called through a variable (e.g. a wrapper's fit_fn argument): look it up
  for (nm in ls(ns, all.names = TRUE)) {
    obj <- get(nm, envir = ns)
    if (is.function(obj) && identical(obj, fn)) return(nm)
  }
  NA_character_
}


# Values, helper functions and packages a function reads from its environment.
# Walks the free variables of the function and, recursively, of its helpers.
#' @keywords internal
.capture_state <- function(fn) {
  values <- list()
  helpers <- list()
  aliases <- character(0)
  not_stored <- character(0)
  packages <- character(0)
  seen <- character(0)
  queue <- list(fn)

  while (length(queue) > 0L) {
    f <- queue[[1L]]
    queue <- queue[-1L]
    packages <- c(packages, .qualified_packages(body(f)))
    globals <- codetools::findGlobals(f, merge = FALSE)

    for (nm in setdiff(c(globals$functions, globals$variables), seen)) {
      seen <- c(seen, nm)
      where <- .find_binding(nm, environment(f))
      # unresolved names are usually columns used in data.table expressions
      if (is.null(where)) next
      val <- get(nm, envir = where)

      if (is.function(val)) {
        val_pkg <- .fn_package(val)
        if (!is.na(val_pkg)) {
          packages <- c(packages, val_pkg)
          # a package function stored under another name in the user's environment
          if (!.is_package_env(where)) {
            ref <- .find_in_namespace(val, val_pkg, nm)
            if (!is.na(ref)) aliases[[nm]] <- sprintf("%s::%s", val_pkg, ref)
          }
        } else {
          helpers[[nm]] <- .deparse_fn(val)
          queue <- c(queue, list(val))
        }
      } else if (.is_package_env(where)) {
        packages <- c(packages, sub("^(package|imports):", "", environmentName(where)))
      } else if (.is_storable(val)) {
        values[[nm]] <- val
      } else {
        not_stored <- c(not_stored, nm)
      }
    }
  }

  list(
    values = values,
    helpers = helpers,
    aliases = aliases,
    not_stored = not_stored,
    packages = setdiff(unique(packages), .provenance_base_packages)
  )
}


# Environment in the chain starting at `env` where `name` is bound, or NULL.
#' @keywords internal
.find_binding <- function(name, env) {
  while (!identical(env, emptyenv())) {
    if (exists(name, envir = env, inherits = FALSE)) return(env)
    env <- parent.env(env)
  }
  NULL
}


#' @keywords internal
.is_package_env <- function(env) {
  isNamespace(env) || identical(env, baseenv()) ||
    grepl("^(package|imports):", environmentName(env))
}


# Packages named in pkg::fun or pkg:::fun calls anywhere in an expression.
#' @keywords internal
.qualified_packages <- function(expr) {
  if (!is.call(expr)) return(character(0))
  head <- expr[[1L]]
  if (is.name(head) && as.character(head) %in% c("::", ":::")) {
    return(as.character(expr[[2L]]))
  }
  unlist(lapply(as.list(expr), .qualified_packages), use.names = FALSE)
}


# Plain data (NULL, atomic vectors, lists and data frames made of them) small enough
# to keep in the object.
#' @keywords internal
.is_storable <- function(x) {
  .is_plain_data(x) && as.numeric(utils::object.size(x)) <= .provenance_max_value_bytes
}


#' @keywords internal
.is_plain_data <- function(x) {
  if (is.null(x) || is.atomic(x)) return(TRUE)
  if (is.list(x) && (is.null(attr(x, "class")) || is.data.frame(x))) {
    return(all(vapply(x, .is_plain_data, logical(1))))
  }
  FALSE
}


#' @keywords internal
.split_storable <- function(x) {
  if (length(x) == 0L) return(list(stored = list(), not_stored = character(0)))
  ok <- vapply(x, .is_storable, logical(1))
  list(stored = x[ok], not_stored = names(x)[!ok])
}


#' @keywords internal
.fit_fn_env <- function(fn) {
  env_fn <- attr(fn, "fit_fn_env")
  if (is.null(env_fn)) return(NULL)
  tryCatch(
    lapply(env_fn(), as.character),
    error = function(e) list(error = conditionMessage(e))
  )
}


#' @keywords internal
.package_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NULL)
}


# Commit of a package installed from a git remote (remotes / pak record it).
#' @keywords internal
.package_git_sha <- function(pkg) {
  sha <- tryCatch(utils::packageDescription(pkg)[["RemoteSha"]], error = function(e) NULL)
  if (is.null(sha) || is.na(sha)) NULL else sha
}


#' @keywords internal
.package_imports <- function(pkg) {
  desc <- tryCatch(utils::packageDescription(pkg), error = function(e) NULL)
  if (!is.list(desc)) return(character(0))
  fields <- unlist(desc[intersect(c("Depends", "Imports"), names(desc))], use.names = FALSE)
  if (length(fields) == 0L) return(character(0))
  pkgs <- trimws(sub("\\(.*\\)", "", unlist(strsplit(fields, ","))))
  setdiff(pkgs[nzchar(pkgs)], c("R", .provenance_base_packages))
}


# Named character vector of installed versions, sorted by package name.
#' @keywords internal
.package_versions <- function(pkgs) {
  pkgs <- sort(unique(setdiff(pkgs, .provenance_base_packages)))
  if (length(pkgs) == 0L) return(NULL)
  versions <- vapply(pkgs, function(p) {
    v <- .package_version(p)
    if (is.null(v)) NA_character_ else v
  }, character(1))
  versions <- versions[!is.na(versions)]
  if (length(versions) == 0L) NULL else versions
}


####
# Restoring and checking
####

# Rebuild a function from its record.
#' @keywords internal
.restore_fn <- function(record, allow_source) {
  if (is.null(record$source)) {
    pkg <- record$package
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("package '%s' is needed to re-run this fit but is not installed", pkg))
    }
    fn <- get(record$name, envir = asNamespace(pkg))
    installed <- .package_version(pkg)
    if (!identical(installed, record$version)) {
      .provenance_warning(sprintf(
        "%s::%s was recorded with version %s but %s is installed; results may differ",
        pkg, record$name, record$version, installed
      ))
    } else if (!identical(digest::digest(.deparse_fn(fn), algo = "xxhash64", serialize = FALSE),
                          record$hash)) {
      .provenance_warning(sprintf(
        "the code of %s::%s differs from the recorded one although the version is the same",
        pkg, record$name
      ))
    }
    return(fn)
  }

  if (!allow_source) {
    stop(paste0(
      "this fit was produced by a function stored as source code; ",
      "set allow_source = TRUE to evaluate it (only for objects from a trusted source)"
    ))
  }
  .warn_on_version_change(record$dependencies)

  # A closure built inside a package function resolves the package internals it
  # calls through the namespace; any other function sees the exports of the
  # packages it used, without touching the search path.
  parent <- if (!is.null(record$package)) {
    asNamespace(record$package)
  } else {
    .dependency_env(names(record$dependencies))
  }
  env <- new.env(parent = parent)
  for (nm in names(record$captured)) assign(nm, record$captured[[nm]], envir = env)
  for (nm in names(record$aliases)) {
    parts <- strsplit(record$aliases[[nm]], "::", fixed = TRUE)[[1L]]
    assign(nm, get(parts[[2L]], envir = asNamespace(parts[[1L]])), envir = env)
  }
  for (nm in names(record$helpers)) {
    assign(nm, eval(parse(text = record$helpers[[nm]]), envir = env), envir = env)
  }
  fn <- eval(parse(text = record$source), envir = env)
  if (!is.null(record$info)) attr(fn, "fit_fn_info") <- record$info
  fn
}


# Environment holding the exports of `pkgs`, on top of the global environment.
#' @keywords internal
.dependency_env <- function(pkgs) {
  env <- new.env(parent = globalenv())
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      .provenance_warning(sprintf("package '%s' used by the fit is not installed", p))
      next
    }
    for (nm in getNamespaceExports(p)) {
      # skip S4 class/method tables, which are not callable
      if (startsWith(nm, ".__")) next
      val <- tryCatch(getExportedValue(p, nm), error = function(e) NULL)
      if (!is.null(val)) assign(nm, val, envir = env)
    }
  }
  env
}


#' @keywords internal
.warn_on_version_change <- function(recorded) {
  for (p in names(recorded)) {
    installed <- .package_version(p)
    if (!is.null(installed) && !identical(installed, recorded[[p]])) {
      .provenance_warning(sprintf(
        "package '%s' was version %s when the fit was recorded and is %s now; results may differ",
        p, recorded[[p]], installed
      ))
    }
  }
}


# Re-run one fit and compare it with the stored rows.
#' @keywords internal
.check_one_fit <- function(se, assay, fit_source, tolerance, allow_source) {
  recipe <- get_fit_provenance(se, assay, fit_source)
  notes <- character(0)
  input_assay <- recipe$input_assay
  if (!is.null(input_assay) && input_assay %in% SummarizedExperiment::assayNames(se) &&
      !identical(.assay_hash(se, input_assay), recipe$run$input_hash)) {
    notes <- c(notes, sprintf("the %s assay changed since the fit", input_assay))
  }

  # Keep provenance warnings as notes; warnings of the fit itself were already
  # raised when the fit was first run.
  refit <- tryCatch(
    withCallingHandlers(
      refit_from_provenance(se, assay, fit_source, allow_source = allow_source),
      gDRutils_provenance_warning = function(w) {
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) e
  )

  stored <- .fit_rows(se, assay, fit_source)
  if (inherits(refit, "error")) {
    return(data.table::data.table(
      assay = assay, fit_source = fit_source, reproducible = NA,
      n_rows = NROW(stored), detail = paste(c(conditionMessage(refit), notes), collapse = "; ")
    ))
  }

  diffs <- .compare_rows(stored, .fit_rows(refit, assay, fit_source), tolerance)
  if (!is.null(recipe$summary_assay)) {
    summary_diffs <- .compare_rows(
      .fit_rows(se, recipe$summary_assay, fit_source),
      .fit_rows(refit, recipe$summary_assay, fit_source),
      tolerance
    )
    if (length(summary_diffs) > 0L) diffs <- c(diffs, paste(recipe$summary_assay, summary_diffs))
  }

  data.table::data.table(
    assay = assay,
    fit_source = fit_source,
    reproducible = length(diffs) == 0L,
    n_rows = NROW(stored),
    detail = paste(c(if (length(diffs) == 0L) "results match" else diffs, notes), collapse = "; ")
  )
}


# Rows of one fit source in an assay, in a stable order and with factors as strings.
#' @keywords internal
.fit_rows <- function(se, assay, fit_source) {
  if (!assay %in% SummarizedExperiment::assayNames(se)) return(data.table::data.table())
  dt <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, assay), row.field = "row", column.field = "column"
  ))
  # an assay without a fit_source column (e.g. combination excess) holds one fit
  if ("fit_source" %in% names(dt)) {
    # compute the filter outside [.data.table, where `fit_source` would name the column;
    # a single symbol in i is evaluated in the calling scope
    keep <- as.character(dt[["fit_source"]]) == fit_source
    dt <- dt[keep]
  }
  for (col in names(dt)) {
    if (is.factor(dt[[col]])) data.table::set(dt, j = col, value = as.character(dt[[col]]))
  }
  data.table::setcolorder(dt, sort(names(dt)))
  key_cols <- names(dt)[vapply(dt, is.character, logical(1))]
  if (length(key_cols) > 0L) data.table::setorderv(dt, key_cols)
  dt
}


# Differences between two sets of fit rows, as short messages (empty when equal).
#' @keywords internal
.compare_rows <- function(stored, refit, tolerance) {
  if (!identical(names(stored), names(refit))) {
    return(sprintf("columns differ: %s", toString(sort(union(
      setdiff(names(stored), names(refit)), setdiff(names(refit), names(stored))
    )))))
  }
  if (NROW(stored) != NROW(refit)) {
    return(sprintf("%d rows stored, %d after re-fitting", NROW(stored), NROW(refit)))
  }
  res <- all.equal(stored, refit, tolerance = tolerance, check.attributes = FALSE)
  if (isTRUE(res)) character(0) else res
}


####
# Reading
####

#' @keywords internal
.find_recipe <- function(prov, assay, fit_source) {
  recipe <- prov$fits[[assay]][[fit_source]]
  if (is.null(recipe)) {
    stop(sprintf("no recipe recorded for fit_source '%s' in assay '%s'", fit_source, assay))
  }
  if (!is.null(recipe$summary_of)) recipe <- .find_recipe(prov, recipe$summary_of, fit_source)
  recipe
}


#' @keywords internal
.summarize_provenance <- function(prov, data_source) {
  rows <- list()
  for (a in names(prov$fits)) {
    for (s in names(prov$fits[[a]])) {
      recipe <- prov$fits[[a]][[s]]
      is_summary <- !is.null(recipe$summary_of)
      if (is_summary) recipe <- .find_recipe(prov, recipe$summary_of, s)
      fn <- if (is_summary) recipe$summary_fn else recipe$fn
      rows[[length(rows) + 1L]] <- data.table::data.table(
        data_source = data_source,
        assay = a,
        fit_source = s,
        fn = .format_fn(fn),
        fit_args = if (is_summary) "" else .format_args(recipe[["fit_args"]]),
        model = if (is.null(fn$info$model)) NA_character_ else fn$info$model
      )
    }
  }
  for (src in names(prov$sources)) {
    rows[[length(rows) + 1L]] <- .summarize_provenance(prov$sources[[src]], src)
  }
  out <- data.table::rbindlist(rows, fill = TRUE)
  if (NROW(out) > 0L && all(is.na(out$data_source))) data.table::set(out, j = "data_source", value = NULL)
  out
}


#' @keywords internal
.format_fn <- function(fn) {
  name <- if (is.null(fn$name)) "anonymous function" else fn$name
  if (is.null(fn$source)) {
    sprintf("%s::%s (%s)", fn$package, name, fn$version)
  } else {
    sprintf("%s (source stored)", name)
  }
}


#' @keywords internal
.format_args <- function(args) {
  if (length(args) == 0L) return("")
  toString(vapply(names(args), function(nm) {
    sprintf("%s = %s", nm, paste(deparse(args[[nm]], width.cutoff = 60L), collapse = " "))
  }, character(1)))
}


####
# Helpers
####

# Hash of an assay, to tell later whether a fit is re-run on the same data. It
# hashes the serialized object, so it is comparable within one R version.
#' @keywords internal
.assay_hash <- function(se, assay) {
  digest::digest(SummarizedExperiment::assay(se, assay), algo = "xxhash64")
}


# Call a fit function on one slice. The data goes in by name of a local variable,
# so a failing fit does not print the whole table in its error call.
#' @keywords internal
.call_fit_fn <- function(fn, dt, fit_args) {
  if (length(fit_args) == 0L) return(fn(dt))
  do.call(fn, c(list(quote(dt)), fit_args))
}


# Name of a function argument as the caller wrote it, or NULL for an expression.
#' @keywords internal
.call_name <- function(expr) {
  if (is.name(expr)) as.character(expr) else NULL
}


# Names of the elements of list(a = fn_a, b = fn_b) as the caller wrote them, keyed
# by element name.
#' @keywords internal
.list_call_names <- function(expr) {
  if (!is.call(expr) || !identical(expr[[1L]], as.name("list"))) return(list())
  args <- as.list(expr)[-1L]
  if (is.null(names(args))) return(list())
  lapply(args, .call_name)
}


# A function written in the session is often called through a wrapper argument
# (fit_fn); when it is bound in the global environment, use that name instead.
#' @keywords internal
.global_name <- function(fn, name) {
  if (!is.null(name) && exists(name, envir = globalenv(), inherits = FALSE) &&
      identical(get(name, envir = globalenv()), fn)) {
    return(name)
  }
  for (nm in ls(globalenv())) {
    if (identical(get(nm, envir = globalenv()), fn)) return(nm)
  }
  name
}

#' @keywords internal
.drop_null <- function(x) x[!vapply(x, is.null, logical(1))]


# Warnings raised by the provenance code carry their own class, so
# check_reproducibility() can report them without the warnings of the fit itself.
#' @keywords internal
.provenance_warning <- function(msg) {
  warning(warningCondition(msg, class = "gDRutils_provenance_warning", call = sys.call(-1L)))
}
