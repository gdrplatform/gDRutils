# get_non_empty_assays

get non empty assays

## Usage

``` r
get_non_empty_assays(mae)
```

## Arguments

- mae:

  MultiAssayExperiment object

## Value

charvec with non-empty experiments

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
get_non_empty_assays(mae)
#> [1] "single-agent"
```
