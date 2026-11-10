####
# Fit function registry
####
#
# A companion to the fit profile registry in fit_profiles.R. Profiles say how
# a data type is sliced; this says which fit functions exist for it.
#
# apply_fit() takes any function, so nothing here is required to fit anything.
# It exists so that a caller offering a *menu* — a skill, a Shiny panel, a
# vignette — can ask what is available instead of hardcoding a list that goes
# stale. The functions are also split across packages (the reference Hill fit
# lives here, the synergy fits in gDRcore), so a single registry is the only
# place where "what can I fit with" has one answer.
#
# Entries hold `package` and `fn_name` rather than the function itself, and
# resolve on demand. That keeps the registry printable and serialisable, does
# not pin a stale copy of a closure across a devtools::load_all(), and matches
# how refit_from_provenance() restores a function by reference.
#
# Public API:
#   get_fit_fns(data_type)   — table of registered functions
#   get_fit_fn(name)         — resolve one to the function itself
#   register_fit_fn(name, …) — add or update an entry at runtime

.fit_fn_registry_env <- new.env(parent = emptyenv())
.fit_fn_registry_state <- new.env(parent = emptyenv())

# Built-ins are registered on first use, not in .onLoad(), for the same reason
# the profiles are loaded lazily: work done during loadNamespace() is charged to
# every package that imports gDRutils, and a failure there takes all of them
# down rather than only the fitting layer.
#' @keywords internal
.ensure_fit_fns <- function() {
  if (!isTRUE(.fit_fn_registry_state$loaded)) {
    # Set the flag first. register_fit_fn() calls .ensure_fit_fns(), so
    # registering the built-ins would otherwise recurse.
    .fit_fn_registry_state$loaded <- TRUE
    .register_builtin_fit_fns()
  }
  invisible(NULL)
}

# gDRutils' own fit functions. gDRcore registers its own from .onLoad(); see
# the note in get_fit_fns() about providers that have not been loaded yet.
#' @keywords internal
.register_builtin_fit_fns <- function() {
  register_fit_fn(
    "fit_drug_response_metrics",
    package = "gDRutils",
    data_type = "single-agent"
  )
  register_fit_fn(
    "fit_drug_response_metrics_4p",
    package = "gDRutils",
    data_type = "single-agent"
  )
  register_fit_fn(
    "fit_drug_response_metrics_xcx",
    package = "gDRutils",
    data_type = "single-agent"
  )
  invisible(NULL)
}


#' Register or update a fit function
#'
#' Adds a fit function to the registry read by \code{\link{get_fit_fns}}, or
#' overwrites an entry of the same name. The entry is available immediately in
#' the current R session and is not saved anywhere: register it again after
#' restarting R, or from the registering package's \code{.onLoad()}.
#'
#' Registering a function does not change how it is called.
#' \code{\link{apply_fit}} takes the function itself, with or without an entry
#' here.
#'
#' @param name string; unique identifier, conventionally the function name.
#' @param package string; package exporting the function.
#' @param fn_name string; the exported name, when it differs from
#'   \code{name}.
#' @param data_type character vector of fit profile names the function is
#'   written for, e.g. \code{"single-agent"}. A vector rather than a string:
#'   a function may serve more than one profile.
#' @param available function of no arguments returning a flag, or \code{NULL}.
#'   Evaluated by \code{\link{get_fit_fns}} to decide whether the function can
#'   run here — \code{musyc_available} for a fit needing Python, for instance.
#'   \code{NULL} means available whenever the package is installed.
#'
#' @return Invisibly, the registered entry.
#'
#' @examples
#' register_fit_fn(
#'   "fit_drug_response_metrics",
#'   package = "gDRutils",
#'   data_type = "single-agent"
#' )
#' get_fit_fn("fit_drug_response_metrics")
#'
#' @seealso \code{\link{get_fit_fns}}, \code{\link{register_fit_profile}}
#' @keywords fit_registry
#' @export
#'
register_fit_fn <- function(name,
                            package,
                            fn_name = name,
                            data_type,
                            available = NULL) {
  checkmate::assert_string(name)
  checkmate::assert_string(package)
  checkmate::assert_string(fn_name)
  checkmate::assert_character(data_type, min.len = 1L, any.missing = FALSE)
  checkmate::assert_function(available, nargs = 0L, null.ok = TRUE)

  entry <- list(
    name = name,
    package = package,
    fn_name = fn_name,
    data_type = data_type,
    available = available
  )
  .ensure_fit_fns()
  .fit_fn_registry_env[[name]] <- entry
  invisible(entry)
}


#' Get a registered fit function
#'
#' Resolves a registry entry to the function itself.
#'
#' @param name string; the registered identifier.
#'
#' @return The fit function. Stops if the name is not registered, or if its
#'   package is not installed.
#'
#' @examples
#' get_fit_fn("fit_drug_response_metrics")
#'
#' @seealso \code{\link{get_fit_fns}}
#' @keywords fit_registry
#' @export
#'
get_fit_fn <- function(name) {
  checkmate::assert_string(name)
  .ensure_fit_fns()
  entry <- .fit_fn_registry_env[[name]]
  if (is.null(entry)) {
    stop(sprintf(
      "Unknown fit function '%s'. Registered: %s",
      name, toString(ls(.fit_fn_registry_env))
    ), call. = FALSE)
  }
  .resolve_fit_fn(entry, must_work = TRUE)
}


# Resolve an entry to its function. Returns NULL instead of raising when
# must_work is FALSE, so get_fit_fns() can list an entry whose package is not
# installed rather than failing the whole query because of one provider.
#' @keywords internal
.resolve_fit_fn <- function(entry, must_work = FALSE) {
  fn <- tryCatch(
    getExportedValue(entry$package, entry$fn_name),
    error = function(e) NULL
  )
  if (is.null(fn) && must_work) {
    stop(sprintf(
      "'%s' is registered from package '%s', which does not export it here.",
      entry$name, entry$package
    ), call. = FALSE)
  }
  fn
}


#' List the registered fit functions
#'
#' The fit functions known to this session, with the description each one
#' carries, for presenting a choice.
#'
#' A function is listed as unavailable when its package is not installed, or
#' when its \code{available} predicate returns \code{FALSE} — a fit needing
#' Python is registered but not offered on a machine without it.
#'
#' \code{model}, \code{description} and \code{reference} are read from the
#' attributes \code{\link{set_fit_fn_info}} attaches, and are \code{NA} for a
#' function that does not carry them.
#'
#' @section Providers register when they load:
#' An entry exists only once the package registering it has been loaded, so a
#' fit function from another package is missing from this table until then.
#' Call \code{requireNamespace()} on the packages you care about first if you
#' need the full catalogue.
#'
#' @param data_type character vector; keep only functions written for these
#'   profiles. \code{NULL} keeps all.
#'
#' @return \code{data.frame}, one row per registered function, with columns
#'   \code{name}, \code{package}, \code{fn_name}, \code{data_type} (comma
#'   separated), \code{model}, \code{description}, \code{reference} and
#'   \code{available}. Zero rows when nothing matches.
#'
#' @examples
#' get_fit_fns()
#' get_fit_fns(data_type = "single-agent")
#'
#' @seealso \code{\link{register_fit_fn}}, \code{\link{get_fit_profiles}}
#' @keywords fit_registry
#' @export
#'
get_fit_fns <- function(data_type = NULL) {
  checkmate::assert_character(data_type, min.len = 1L, null.ok = TRUE,
                              any.missing = FALSE)
  .ensure_fit_fns()
  entries <- as.list(.fit_fn_registry_env)
  if (!is.null(data_type)) {
    entries <- Filter(
      function(e) length(intersect(e$data_type, data_type)) > 0L,
      entries
    )
  }
  if (length(entries) == 0L) {
    return(.empty_fit_fn_table())
  }

  rows <- lapply(entries[order(names(entries))], function(e) {
    fn <- .resolve_fit_fn(e)
    info <- if (is.null(fn)) NULL else attr(fn, "fit_fn_info")
    data.frame(
      name = e$name,
      package = e$package,
      fn_name = e$fn_name,
      data_type = toString(e$data_type),
      model = .info_field(info, "model"),
      description = .info_field(info, "description"),
      reference = .info_field(info, "reference"),
      available = .fit_fn_available(e, fn),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}


#' @keywords internal
.info_field <- function(info, key) {
  value <- info[[key]]
  if (is.null(value)) NA_character_ else as.character(value)[1]
}


# Unavailable rather than an error: a predicate that cannot answer is not the
# same as a function that cannot run, but neither should be offered.
#' @keywords internal
.fit_fn_available <- function(entry, fn) {
  if (is.null(fn)) {
    return(FALSE)
  }
  if (is.null(entry$available)) {
    return(TRUE)
  }
  isTRUE(tryCatch(entry$available(), error = function(e) FALSE))
}


#' @keywords internal
.empty_fit_fn_table <- function() {
  data.frame(
    name = character(0),
    package = character(0),
    fn_name = character(0),
    data_type = character(0),
    model = character(0),
    description = character(0),
    reference = character(0),
    available = logical(0),
    stringsAsFactors = FALSE
  )
}
