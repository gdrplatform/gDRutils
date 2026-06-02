# Demote a metadata field in the `rowData` or `colData` of a `SummarizedExperiment` object to a nested field of a `BumpyMatrix` assay.

Demote a metadata field in the `rowData` or `colData` of a
`SummarizedExperiment` object to a nested field of a `BumpyMatrix`
assay.

## Usage

``` r
demote_fields(se, fields)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

- fields:

  Character vector of metadata fields to demote as nested columns.

## Value

A `SummarizedExperiment` object with new dimensions resulting from
demoting given `fields` to nested columns.

## Details

Revert this operation using `promote_fields`.

## See also

promote_fields

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
se <- promote_fields(se, "ReadoutValue", 2)
#> Warning: dropping assay 'Controls' from the final SE
#> Warning: dropping assay 'Normalized' from the final SE
#> Warning: dropping assay 'Averaged' from the final SE
#> Warning: dropping assay 'Metrics' from the final SE
demote_fields(se, "ReadoutValue")
#> class: SummarizedExperiment 
#> dim: 10 10 
#> metadata(5): identifiers experiment_metadata Keys fit_parameters
#>   .internal
#> assays(1): RawTreated
#> rownames(10): 1 14 ... 111 127
#> rowData names(4): Gnumber DrugName drug_moa Duration
#> colnames(10): 1 271 ... 2161 2431
#> colData names(4): clid CellLineName Tissue ReferenceDivisionTime
```
