# is_exp_empty

check if experiment (SE) is empty

## Usage

``` r
is_exp_empty(exp)
```

## Arguments

- exp:

  SummarizedExperiment object.

## Value

logical

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
is_exp_empty(se)
#> [1] FALSE
```
