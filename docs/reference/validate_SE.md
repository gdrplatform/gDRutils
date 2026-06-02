# Validate SummarizedExperiment object

Function validates correctness of SE by checking multiple cases:

- detection of duplicated rowData/colData,

- incompatibility of rownames/colnames,

- occurrence of necessary assays,

- detection of mismatch of CLIDs inside colData and colnames (different
  order),

- correctness of metadata names.

## Usage

``` r
validate_SE(se, expect_single_agent = FALSE)
```

## Arguments

- se:

  SummarizedExperiment object produced by the gDR pipeline

- expect_single_agent:

  a logical indicating if the function should check whether the
  SummarizedExperiment is single-agent data

## Value

`NULL` invisibly if the SummarizedExperiment is valid. Throws an error
if the SummarizedExperiment is not valid.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
validate_SE(se)
```
