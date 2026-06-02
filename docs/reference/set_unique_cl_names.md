# Set Unique Parental Identifiers

This function sets the `CellLineName` field in `colData` to be unique by
appending the `clid` in parentheses for duplicates.

## Usage

``` r
set_unique_cl_names(se)
```

## Arguments

- se:

  A SummarizedExperiment object.

## Value

A SummarizedExperiment object with unique `CellLineName` in `colData`.

## Examples

``` r
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(1:4, ncol = 2)),
  colData = S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2"))
)
se <- set_unique_cl_names(se)
```
