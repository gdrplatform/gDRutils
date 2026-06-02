# Set unique primary cell line identifiers in the table

This function sets the primary cell line field in data.frame-like object
to be unique by appending the secondary cell line field in parentheses
for duplicates.

## Usage

``` r
set_unique_cl_names_dt(
  dt,
  primary_name = get_env_identifiers("cellline_name"),
  secondary_name = get_env_identifiers("cellline"),
  sep = " "
)
```

## Arguments

- dt:

  data.table, data.frame or DFrame with the data

- primary_name:

  string with the name of the primary cell line field

- secondary_name:

  string with the name of the secondary cell line field

- sep:

  string with separator added before suffix

## Value

fixed input table with unique primary cell line field in dt

## Examples

``` r
col_data <- S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2"))
col_data <- set_unique_cl_names_dt(col_data)
```
