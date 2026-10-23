# Build a SummarizedExperiment from a long table

Turns a long table, one row per measurement, into an SE with a single
BumpyMatrix assay. Meant for assays that do not come out of the gDR
pipeline, e.g. compounds x targets in a biochemical screen; register a
profile with
[`register_fit_profile`](https://gdrplatform.github.io/gDRstyle/reference/register_fit_profile.md)
and pass the result to
[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md).

## Usage

``` r
convert_dt_to_se(dt, row_cols, col_cols, assay_name)
```

## Arguments

- dt:

  data.table with one row per measurement.

- row_cols:

  character vector of columns identifying an SE row.

- col_cols:

  character vector of columns identifying an SE column.

- assay_name:

  string; name of the assay.

## Value

`SummarizedExperiment` with `row_cols` in `rowData`, `col_cols` in
`colData` and the remaining columns nested in assay `assay_name`.

## Details

Several `row_cols` or `col_cols` are pasted with `"_"` into the row or
column names. Combinations absent from `dt` come out as empty cells.

## Examples

``` r
dt <- data.table::CJ(compound = c("cA", "cB"), target = "t1", Concentration = c(0.1, 1))
dt$Activity <- c(90, 40, 80, 20)
convert_dt_to_se(dt, "compound", "target", "Biochemical")
#> class: SummarizedExperiment 
#> dim: 2 1 
#> metadata(0):
#> assays(1): Biochemical
#> rownames(2): cA cB
#> rowData names(1): compound
#> colnames(1): t1
#> colData names(1): target
```
