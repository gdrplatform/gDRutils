# Standardize MAE by switching from custom identifiers into gDR-default

Standardize MAE by switching from custom identifiers into gDR-default

## Usage

``` r
standardize_mae(mae, use_default = TRUE)
```

## Arguments

- mae:

  a MultiAssayExperiment object with drug-response data generate by gDR
  pipeline

- use_default:

  boolean indicating whether or not to use default identifiers for
  standardization

## Value

mae a MultiAssayExperiment with default gDR identifiers

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
S4Vectors::metadata(mae[[1]])$identifiers$drug <- "druug"
standardize_mae(mae)
#> Warning: overwriting existing metadata entry: 'identifiers'
#> A MultiAssayExperiment object of 1 listed
#>  experiment with a user-defined name and respective class.
#>  Containing an ExperimentList class object of length 1:
#>  [1] single-agent: SummarizedExperiment with 10 rows and 10 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
