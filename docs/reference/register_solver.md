# Register a solver for the Hill fit

Records a solver under a name so it can be listed with
[`get_solvers()`](https://gdrplatform.github.io/gDRstyle/reference/get_solvers.md).
Any function with the
[`solver_drc`](https://gdrplatform.github.io/gDRstyle/reference/solver_drc.md)
signature can be passed to `fit_drug_response_metrics(solver = )`
whether or not it is registered; the registry is a menu.

## Usage

``` r
register_solver(
  name,
  fn,
  description = NA_character_,
  available = function() TRUE
)
```

## Arguments

- name:

  string, the solver's name.

- fn:

  the solver function.

- description:

  string, what it does.

- available:

  function of no arguments returning `TRUE` when the solver can run
  (e.g. its package is installed). Defaults to always available.

## Value

`NULL`, invisibly.

## Examples

``` r
register_solver("drc", solver_drc, "drc::drm, L-BFGS-B within the model bounds")
names(get_solvers())
#> [1] "drc"
```
