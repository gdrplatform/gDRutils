# Set Unique Drug Names

This function sets the `DrugName`, `DrugName_2`, and `DrugName_3` fields
in `rowData` to be unique by appending the corresponding `Gnumber`,
`Gnumber_2`, and `Gnumber_3` in parentheses for duplicates.

## Usage

``` r
set_unique_drug_names(se)
```

## Arguments

- se:

  A SummarizedExperiment object.

## Value

A SummarizedExperiment object with unique `DrugName` fields in
`rowData`.

## Examples

``` r
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(1:9, ncol = 3)),
  rowData = S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA", "DrugB"),
  Gnumber = c("G1", "G2", "G5"),
  DrugName_2 = c("DrugC", "DrugC", "DrugD"),
  Gnumber_2 = c("G3", "G4", "G5")
))
se <- set_unique_drug_names(se)
```
