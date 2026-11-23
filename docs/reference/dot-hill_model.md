# Describe the Hill model a solver has to fit

Describe the Hill model a solver has to fit

## Usage

``` r
.hill_model(x_0, start, lower, upper)
```

## Arguments

- x_0:

  numeric fixed upper asymptote, or `NA` to fit it (4-parameter model).

- start, lower, upper:

  named numeric vectors over the free parameters, in the order of
  `free`.

## Value

list with `free`, `fixed`, `start`, `lower`, `upper`.
