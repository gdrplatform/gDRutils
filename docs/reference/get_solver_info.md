# Describe the registered solvers

The menu
[`get_solvers()`](https://gdrplatform.github.io/gDRstyle/reference/get_solvers.md)
returns, with what each solver does and whether it can run here, so a
caller can offer a choice without knowing the solvers in advance.

## Usage

``` r
get_solver_info(available_only = FALSE)
```

## Arguments

- available_only:

  logical; if `TRUE`, only solvers whose `available()` check passes.
  Defaults to `FALSE`, so a solver that cannot run is still listed, and
  marked.

## Value

`data.frame` with one row per solver: `name`, `description` and
`available`.

## Examples

``` r
get_solver_info()
#>   name                                                              description
#> 1  drc drc::drm with LL.3u/LL.4 within the model bounds (stats::optim L-BFGS-B)
#>   available
#> 1      TRUE
```
