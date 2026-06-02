# Identify unique metadata fields from a list of `SummarizedExperiment`s

Identify unique metadata fields from a list of `SummarizedExperiment`s

## Usage

``` r
identify_unique_se_metadata_fields(SElist)
```

## Arguments

- SElist:

  named list of `SummarizedExperiment`s

## Value

character vector of unique names of metadata

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
SElist <- list(
  se,
  se
)
identify_unique_se_metadata_fields(SElist)
#> [1] "identifiers"         "experiment_metadata" "Keys"               
#> [4] "fit_parameters"      ".internal"          
```
