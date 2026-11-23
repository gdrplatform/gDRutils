# List the registered solvers

List the registered solvers

## Usage

``` r
get_solvers(available_only = TRUE)
```

## Arguments

- available_only:

  logical; if `TRUE` (default), only solvers whose `available()` check
  passes.

## Value

named list of solver functions.

## Examples

``` r
names(get_solvers())
#> [1] "drc"
```
