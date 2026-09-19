#' onload function
#'
#' @param libname library name
#' @param pkgname package name
#' @noRd
.onLoad <- function(libname, pkgname) {
  # Fit profiles - load from inst/extdata/fit_profiles.json
  .load_fit_profiles()
}
