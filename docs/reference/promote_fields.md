# Promote a nested field to be represented as a metadata field of the `SummarizedExperiment` as either the `rowData` or `colData`.

Promote a nested field to be represented as a metadata field of the
`SummarizedExperiment` as either the `rowData` or `colData`.

## Usage

``` r
promote_fields(se, fields, MARGIN = c(1, 2))
```

## Arguments

- se:

  `SummarizedExperiment` object.

- fields:

  Character vector of nested fields in a `BumpyMatrix` object to promote
  to metadata fields of a `se`.

- MARGIN:

  Numeric of values `1` or `2` indicating whether to promote fields to
  rows or columns respectively.

## Value

A `SummarizedExperiment` object with new dimensions resulting from
promoting given `fields`.

## Details

Revert this operation using `demote_fields`.

## See also

demote_fields

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
se <- promote_fields(se, "ReadoutValue", 2)
#> Warning: dropping assay 'Controls' from the final SE
#> Warning: dropping assay 'Normalized' from the final SE
#> Warning: dropping assay 'Averaged' from the final SE
#> Warning: dropping assay 'Metrics' from the final SE
```
