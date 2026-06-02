# Geometric mean

Auxiliary function for calculating geometric mean with possibility to
handle -Inf

## Usage

``` r
geometric_mean(x, fixed = TRUE, maxlog10Concentration = 1)
```

## Arguments

- x:

  numeric vector

- fixed:

  flag should be add fix for -Inf

- maxlog10Concentration:

  numeric value needed to calculate minimal value

## Value

numeric vector

## Examples

``` r
geometric_mean(c(2, 8))
#> [1] 4
```
