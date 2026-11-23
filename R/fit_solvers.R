####
# Solvers for the Hill fit
####
#
# A solver finds the parameters of the Hill curve, the response falling from x_0 at no
# drug to x_inf at saturating drug, with half of that drop at the concentration ec50 and
# steepness h, for one dose-response series. Everything around it (preprocessing,
# capping, the constant-fit test and the metrics) is fit_drug_response_metrics(), so
# changing the optimiser does not mean rewriting the fit.
#
# A solver takes the concentrations, the responses, the model built by .hill_model() and
# a list of its own settings, and returns the fitted h, x_inf, x_0 and ec50, the residual
# sum of squares, whether it converged and a message. A solver that cannot fit reports
# converged = FALSE rather than throwing.


#' Describe the Hill model a solver has to fit
#'
#' @param x_0 numeric fixed upper asymptote, or \code{NA} to fit it (4-parameter model).
#' @param start,lower,upper named numeric vectors over the free parameters, in the
#'   order of \code{free}.
#' @return list with \code{free}, \code{fixed}, \code{start}, \code{lower}, \code{upper}.
#' @keywords internal
.hill_model <- function(x_0, start, lower, upper) {
  free <- if (is.na(x_0)) c("h", "x_inf", "x_0", "ec50") else c("h", "x_inf", "ec50")
  stopifnot(length(start) == length(free), length(lower) == length(free),
            length(upper) == length(free))
  list(
    free = free,
    fixed = if (is.na(x_0)) list() else list(x_0 = x_0),
    start = stats::setNames(as.numeric(start), free),
    lower = stats::setNames(as.numeric(lower), free),
    upper = stats::setNames(as.numeric(upper), free)
  )
}


#' Hill fit solver using drc
#'
#' The default solver of \code{\link{fit_drug_response_metrics}}: \code{drc::drm()} with
#' \code{drc::LL.3u()} (upper asymptote fixed) or \code{drc::LL.4()} (upper asymptote
#' fitted), within the bounds of the model. With bounds, drc optimises with
#' \code{stats::optim(method = "L-BFGS-B")}.
#'
#' A solver is any function with this signature returning the same list, and can be
#' passed to \code{fit_drug_response_metrics(solver = )} or through
#' \code{apply_fit(fit_args = list(solver = ))}.
#'
#' @param conc numeric vector of concentrations (not logged).
#' @param x numeric vector of responses, the same length as \code{conc}.
#' @param model list describing the model to fit: \code{free} parameter names,
#'   \code{fixed} values (\code{x_0} for the 3-parameter model), and \code{start},
#'   \code{lower}, \code{upper} over the free parameters.
#' @param control named list of \code{drc::drmc()} settings overriding the defaults
#'   (\code{relTol = 1e-4}, \code{errorm = FALSE}, \code{noMessage = TRUE},
#'   \code{rmNA = TRUE}). With the bounds of the model drc optimises with L-BFGS-B, which
#'   ignores \code{relTol}.
#'
#' @return list with \code{params} (named \code{h}, \code{x_inf}, \code{x_0},
#'   \code{ec50}), \code{rss}, \code{converged} and \code{message}.
#'
#' @examples
#' model <- list(free = c("h", "x_inf", "ec50"), fixed = list(x_0 = 1),
#'               start = c(h = 2, x_inf = 0.4, ec50 = 0.1),
#'               lower = c(h = 0.1, x_inf = 0, ec50 = 1e-4),
#'               upper = c(h = 5, x_inf = 1, ec50 = 100))
#' solver_drc(c(0.001, 0.01, 0.1, 1, 10), c(0.95, 0.8, 0.5, 0.2, 0.1), model)
#'
#' @keywords fit_curves
#' @export
solver_drc <- function(conc, x, model, control = list()) {
  checkmate::assert_numeric(conc)
  checkmate::assert_numeric(x, len = length(conc))
  checkmate::assert_list(control, names = "unique")

  four_param <- is.null(model$fixed$x_0)
  fct <- if (four_param) {
    drc::LL.4(names = c("h", "x_inf", "x_0", "ec50"))
  } else {
    drc::LL.3u(upper = model$fixed$x_0, names = c("h", "x_inf", "ec50"))
  }
  ctrl <- utils::modifyList(
    list(relTol = 1e-04, errorm = FALSE, noMessage = TRUE, rmNA = TRUE), control
  )

  # drc >= 4.0-0 writes its na.action argument into the global option and never
  # restores it, so one fit changes how every later model.frame() in the session
  # behaves. Restore it here, and keep the name unqualified: stored as
  # "stats::na.omit" the option breaks stats::lm() with
  # could not find function "stats::na.omit".
  old_na_action <- getOption("na.action")
  on.exit(options(na.action = old_na_action), add = TRUE)

  fit <- tryCatch(
    drc::drm(x ~ conc,
             data = data.table::data.table(x = x, conc = conc),
             logDose = NULL,
             fct = fct,
             start = unname(model$start),
             lowerl = unname(model$lower),
             upperl = unname(model$upper),
             control = do.call(drc::drmc, ctrl),
             na.action = na.omit),
    error = function(e) e
  )
  # With errorm = FALSE a failed optimisation is not an error: drm() warns and returns a
  # plain list with convergence = FALSE, which has no coefficients.
  if (inherits(fit, "error") || !inherits(fit, "drc") || !isTRUE(fit$fit$convergence)) {
    msg <- if (inherits(fit, "error")) conditionMessage(fit) else "drc::drm() did not converge"
    return(list(params = c(h = NA_real_, x_inf = NA_real_, x_0 = NA_real_, ec50 = NA_real_),
                rss = NA_real_, converged = FALSE, message = msg))
  }

  coefs <- stats::coef(fit)
  get_coef <- function(nm) unname(coefs[[paste0(nm, ":(Intercept)")]])
  list(
    params = c(h = get_coef("h"), x_inf = get_coef("x_inf"),
               x_0 = if (four_param) get_coef("x_0") else model$fixed$x_0,
               ec50 = get_coef("ec50")),
    rss = sum(stats::residuals(fit)^2, na.rm = TRUE),
    converged = TRUE,
    message = NA_character_
  )
}


####
# Solver registry
####

.solver_registry <- new.env(parent = emptyenv())
.solver_registry_state <- new.env(parent = emptyenv())

# Built-ins are registered on first use, as the fit functions are (see
# .ensure_fit_fns()): nothing runs during loadNamespace().
#' @keywords internal
.ensure_solvers <- function() {
  if (!isTRUE(.solver_registry_state$loaded)) {
    .solver_registry_state$loaded <- TRUE   # first: register_solver() calls back here
    .register_default_solvers()
  }
  invisible(NULL)
}


#' Register a solver for the Hill fit
#'
#' Records a solver under a name so it can be listed with \code{get_solvers()}. Any
#' function with the \code{\link{solver_drc}} signature can be passed to
#' \code{fit_drug_response_metrics(solver = )} whether or not it is registered; the
#' registry is a menu.
#'
#' @param name string, the solver's name.
#' @param fn the solver function.
#' @param description string, what it does.
#' @param available function of no arguments returning \code{TRUE} when the solver can
#'   run (e.g. its package is installed). Defaults to always available.
#'
#' @return \code{NULL}, invisibly.
#'
#' @examples
#' register_solver("drc", solver_drc, "drc::drm, L-BFGS-B within the model bounds")
#' names(get_solvers())
#'
#' @keywords fit_curves
#' @export
register_solver <- function(name, fn, description = NA_character_,
                            available = function() TRUE) {
  checkmate::assert_string(name)
  checkmate::assert_function(fn, args = c("conc", "x", "model", "control"))
  checkmate::assert_string(description, na.ok = TRUE)
  checkmate::assert_function(available)
  .ensure_solvers()
  assign(name, list(fn = fn, description = description, available = available),
         envir = .solver_registry)
  invisible(NULL)
}


#' List the registered solvers
#'
#' @param available_only logical; if \code{TRUE} (default), only solvers whose
#'   \code{available()} check passes.
#'
#' @return named list of solver functions.
#'
#' @examples
#' names(get_solvers())
#'
#' @keywords fit_curves
#' @export
get_solvers <- function(available_only = TRUE) {
  checkmate::assert_flag(available_only)
  .ensure_solvers()
  out <- lapply(ls(.solver_registry), function(nm) get(nm, envir = .solver_registry))
  names(out) <- ls(.solver_registry)
  if (available_only) {
    out <- Filter(function(s) isTRUE(tryCatch(s$available(), error = function(e) FALSE)), out)
  }
  lapply(out, `[[`, "fn")
}


#' Describe the registered solvers
#'
#' The menu \code{get_solvers()} returns, with what each solver does and whether it
#' can run here, so a caller can offer a choice without knowing the solvers in
#' advance.
#'
#' @param available_only logical; if \code{TRUE}, only solvers whose
#'   \code{available()} check passes. Defaults to \code{FALSE}, so a solver that
#'   cannot run is still listed, and marked.
#'
#' @return \code{data.frame} with one row per solver: \code{name},
#'   \code{description} and \code{available}.
#'
#' @examples
#' get_solver_info()
#'
#' @keywords fit_curves
#' @export
get_solver_info <- function(available_only = FALSE) {
  checkmate::assert_flag(available_only)
  .ensure_solvers()
  nms <- ls(.solver_registry)
  entries <- lapply(nms, function(nm) get(nm, envir = .solver_registry))
  out <- data.frame(
    name = nms,
    description = vapply(entries, function(s) as.character(s$description), character(1)),
    available = vapply(entries, function(s) {
      isTRUE(tryCatch(s$available(), error = function(e) FALSE))
    }, logical(1)),
    stringsAsFactors = FALSE
  )
  if (available_only) {
    out <- out[out$available, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}


#' @keywords internal
.register_default_solvers <- function() {
  register_solver("drc", solver_drc,
                  "drc::drm with LL.3u/LL.4 within the model bounds (stats::optim L-BFGS-B)")
}
