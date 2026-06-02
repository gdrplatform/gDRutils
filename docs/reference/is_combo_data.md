# Checks if `se` is combo dataset.

Checks if `se` is combo dataset.

## Usage

``` r
is_combo_data(se)
```

## Arguments

- se:

  SummarizedExperiment

## Value

logical

## Examples

``` r
se <- get_synthetic_data("finalMAE_combo_matrix.qs2")[[1]]
is_combo_data(se)
#> [1] TRUE
se <- get_synthetic_data("finalMAE_combo_matrix.qs2")[[2]]
is_combo_data(se)
#> [1] FALSE
se <- get_synthetic_data("finalMAE_small.qs2")[[1]]
is_combo_data(se)
#> [1] FALSE
```
