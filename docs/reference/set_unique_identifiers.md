# Set Unique Identifiers in MultiAssayExperiment

This function sets the `CellLineName` in `colData` and `DrugName` fields
in `rowData` to be unique for each `SummarizedExperiment` in a
`MultiAssayExperiment`.

## Usage

``` r
set_unique_identifiers(mae)
```

## Arguments

- mae:

  A MultiAssayExperiment object.

## Value

A MultiAssayExperiment object with unique identifiers.

## Examples

``` r
se1 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(1:4, ncol = 2)),
  colData = S4Vectors::DataFrame(parental_identifier = c("ID1", "ID1"), clid = c("C1", "C2")),
  rowData = S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA"), Gnumber = c("G1", "G2"))
)
rownames(SummarizedExperiment::colData(se1)) <- c("cl1", "cl2")
rownames(SummarizedExperiment::rowData(se1)) <- c("g1", "g")
se2 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(5:8, ncol = 2)),
  colData = S4Vectors::DataFrame(parental_identifier = c("ID2", "ID2"), clid = c("C3", "C4")),
  rowData = S4Vectors::DataFrame(DrugName = c("DrugB", "DrugB"), Gnumber = c("G3", "G4"))
)
rownames(SummarizedExperiment::colData(se2)) <- c("cl3", "cl4")
rownames(SummarizedExperiment::rowData(se2)) <- c("g3", "g4")
mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(se1 = se1, se2 = se2))
mae <- set_unique_identifiers(mae)
```
