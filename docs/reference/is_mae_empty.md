# is_mae_empty

check if all mae experiments are empty

## Usage

``` r
is_mae_empty(mae)
```

## Arguments

- mae:

  MultiAssayExperiment object

## Value

logical

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
is_mae_empty(mae)
#> [1] FALSE
```
