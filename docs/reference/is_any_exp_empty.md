# is_any_exp_empty

check if any experiment is empty

## Usage

``` r
is_any_exp_empty(mae)
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
is_any_exp_empty(mae)
#> [1] FALSE
```
