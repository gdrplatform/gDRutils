# Flatten a table

Flatten a stacked table into a wide format.

## Usage

``` r
flatten(tbl, groups, wide_cols, sep = "_")
```

## Arguments

- tbl:

  table to flatten.

- groups:

  character vector of column names representing uniquifying groups in
  expansion.

- wide_cols:

  character vector of column names to flatten.

- sep:

  string representing separator between `wide_cols` columns, used in
  column renaming. Defaults to `"_"`.

## Value

table of flattened data as defined by `wide_cols`.

## Details

flattened columns will be named with original column names prefixed by
`wide_cols` columns, concatenated together and separated by `sep`.

A common use case for this function is when a flattened version of the
`"Metrics"` assay is desired.

## See also

convert_se_assay_to_dt

## Examples

``` r
n <- 4
m <- 5
grid <- expand.grid(normalization_type = c("GR", "RV"),
  source = c("GDS", "GDR"))
repgrid <- data.table::rbindlist(rep(list(grid), m))
repgrid$wide <- seq(m * n)
repgrid$id <- rep(LETTERS[1:m], each = n)

groups <- colnames(grid)
wide_cols <- c("wide")

flatten(repgrid, groups = groups, wide_cols = wide_cols)
#>    GR_GDS_wide     id RV_GDS_wide GR_GDR_wide RV_GDR_wide
#>          <int> <char>       <int>       <int>       <int>
#> 1:           1      A           2           3           4
#> 2:           5      B           6           7           8
#> 3:           9      C          10          11          12
#> 4:          13      D          14          15          16
#> 5:          17      E          18          19          20
```
