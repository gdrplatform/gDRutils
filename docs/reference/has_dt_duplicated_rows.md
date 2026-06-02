# check if data.table contains duplicated data

An auxiliary function that checks for duplicates in the data.table (or
its subset)

## Usage

``` r
has_dt_duplicated_rows(dt, col_names = NULL)
```

## Arguments

- dt:

  data.table

- col_names:

  charvec with columns to be used for subsetting

## Value

logical flag indicating if a dt contains duplicated rows or not

## Examples

``` r
dt <- data.table::data.table(a = c(1, 2, 3), b = c(3, 2, 2))
has_dt_duplicated_rows(dt, "b")
#> [1] TRUE
```
