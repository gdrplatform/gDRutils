# Get synthetic data from gDRtestData package

Get synthetic data from gDRtestData package

## Usage

``` r
get_synthetic_data(qs)
```

## Arguments

- qs:

  dataset name or qs2 filename (e.g. `"small"` or
  `"finalMAE_small.qs2"`)

## Value

loaded data

## Examples

``` r
get_synthetic_data("finalMAE_small.qs2")
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
