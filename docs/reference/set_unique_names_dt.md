# Set unique primary identifiers in the data.frame-like objects

This function sets the primary field in the data.frame-like objects to
be unique by appending the secondary field in parentheses for
duplicates.

## Usage

``` r
set_unique_names_dt(dt, primary_name, secondary_name, sep = " ")
```

## Arguments

- dt:

  data.table, data.frame or DFrame with data

- primary_name:

  string with the name of the primary field

- secondary_name:

  string with the name of the secondary field

- sep:

  string with separator added before suffix

## Value

fixed input table with unique primary field in the table

## Examples

``` r
col_data <- S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2"))
col_data <- set_unique_names_dt(col_data, primary_name = "CellLineName", secondary_name = "clid")
```
