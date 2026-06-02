# Helper function to find duplicated rows

Helper function to find duplicated rows

## Usage

``` r
get_duplicated_rows(x, col_names = NULL, output = "index")
```

## Arguments

- x:

  DataFrame or data.table

- col_names:

  character vector, columns in which duplication are searched for

- output:

  string with the output format to be returned - one of "index" (index
  of duplicates) or "data" (subset of input data with duplicates)

## Value

integer vector or data.table with duplicated rows

## Examples

``` r
dt <- data.table::data.table(a = c(1, 2, 3), b = c(3, 2, 2))
get_duplicated_rows(dt, "b")
#> [1] 2 3
get_duplicated_rows(dt, "b", output = "data")
#>        a     b
#>    <num> <num>
#> 1:     2     2
#> 2:     3     2
```
