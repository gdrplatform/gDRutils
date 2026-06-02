# Validate MultiAssayExperiment object

Function validates correctness of SE included in MAE by checking
multiple cases:

- detection of duplicated rowData/colData,

- incompatibility of rownames/colnames,

- occurrence of necessary assays,

- detection of mismatch of CLIDs inside colData and colnames (different
  order),

- correctness of metadata names.

## Usage

``` r
validate_MAE(mae)
```

## Arguments

- mae:

  MultiAssayExperiment object produced by the gDR pipeline

## Value

`NULL` invisibly if the MultiAssayExperiment is valid. Throws an error
if the MultiAssayExperiment is not valid.

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
validate_MAE(mae)
```
