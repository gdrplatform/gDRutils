# Set unique primary drug identifiers in the table

This function sets the primary drug field(s) in data.frame-like object
to be unique by appending the secondary drug field(s) in parentheses for
duplicates. By default `DrugName`, `DrugName_2`, and `DrugName_3` are
primary drug fields, while `Gnumber`, `Gnumber_2`, and `Gnumber_3` are
their respective secondary drug fields.

## Usage

``` r
set_unique_drug_names_dt(
  dt,
  primary_names = unlist(get_env_identifiers()[(c("drug_name", "drug_name2",
    "drug_name3"))]),
  secondary_names = unlist(get_env_identifiers()[(c("drug", "drug2", "drug3"))]),
  sep = " "
)
```

## Arguments

- dt:

  data.table, data.frame or DFrame with the data

- primary_names:

  charvec with the names of the primary drug field(s)

- secondary_names:

  charvec with the name of the secondary drug field(s)

- sep:

  string with separator added before suffix

## Value

fixed input table with unique primary drug field in dt

## Examples

``` r
row_data <- S4Vectors::DataFrame(
  DrugName = c("DrugA", "DrugA", "DrugB"),
  Gnumber = c("G1", "G2", "G5"),
  DrugName_2 = c("DrugC", "DrugC", "DrugD"),
  Gnumber_2 = c("G3", "G4", "G5")
)
row_data <- set_unique_drug_names_dt(row_data)
```
