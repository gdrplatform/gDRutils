# Split big table

Helper function for saving big tables in an Excel file. Excel has a
sheet size limit, if the table is too large it will not be possible to
save such a file. This function allows you to split the table into
smaller parts so that saving can be possible

## Usage

``` r
split_big_table_for_xlsx(dt_list, max_row = 1e+06, max_col = 16000)
```

## Arguments

- dt_list:

  list of data.tables. Each data.table will be checked and split if meet
  the criteria

- max_row:

  integer defining the maximum number of rows in one sheet, the rows
  will be divided into portions of this size. Default value, 1 000 000,
  is based on excel limit - 1 048 576 with extra safety margin

- max_col:

  integer defining the maximum number of columns in one sheet, the
  columns will be divided into portions of this size. Default value, 16
  000, is based on excel limit - 16 384 with extra safety margin

## Value

list of data.tables

## Examples

``` r
too_large_dt <- list(data.table::data.table(matrix(seq_len(300)), nrow = 10))
split_big_table_for_xlsx(too_large_dt, max_row = 250)
#> $`_1`
#>         V1  nrow
#>      <int> <num>
#>   1:     1    10
#>   2:     2    10
#>   3:     3    10
#>   4:     4    10
#>   5:     5    10
#>  ---            
#> 246:   246    10
#> 247:   247    10
#> 248:   248    10
#> 249:   249    10
#> 250:   250    10
#> 
#> $`_2`
#>        V1  nrow
#>     <int> <num>
#>  1:   251    10
#>  2:   252    10
#>  3:   253    10
#>  4:   254    10
#>  5:   255    10
#>  6:   256    10
#>  7:   257    10
#>  8:   258    10
#>  9:   259    10
#> 10:   260    10
#> 11:   261    10
#> 12:   262    10
#> 13:   263    10
#> 14:   264    10
#> 15:   265    10
#> 16:   266    10
#> 17:   267    10
#> 18:   268    10
#> 19:   269    10
#> 20:   270    10
#> 21:   271    10
#> 22:   272    10
#> 23:   273    10
#> 24:   274    10
#> 25:   275    10
#> 26:   276    10
#> 27:   277    10
#> 28:   278    10
#> 29:   279    10
#> 30:   280    10
#> 31:   281    10
#> 32:   282    10
#> 33:   283    10
#> 34:   284    10
#> 35:   285    10
#> 36:   286    10
#> 37:   287    10
#> 38:   288    10
#> 39:   289    10
#> 40:   290    10
#> 41:   291    10
#> 42:   292    10
#> 43:   293    10
#> 44:   294    10
#> 45:   295    10
#> 46:   296    10
#> 47:   297    10
#> 48:   298    10
#> 49:   299    10
#> 50:   300    10
#>        V1  nrow
#>     <int> <num>
#> 
```
