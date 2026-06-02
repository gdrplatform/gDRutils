# Update drug name with additional variables

Concatenates the values of specified additional variables to the
existing drug identifier columns in a data.table, using the variables
defined in `get_env_identifiers`.

## Usage

``` r
update_drug_name(dt, additional_vars)
```

## Arguments

- dt:

  A data.table containing drug-response information, including drug
  identifier columns (e.g., `DrugName`, `Gnumber`) and the
  `additional_vars`.

- additional_vars:

  Character vector of column names (variables) to merge into the drug
  identifier columns.

## Value

A copy of the input data.table `dt` with the relevant drug identifier
columns updated to include the additional variable information in the
format: `Identifier (variable = value)`.

## Examples

``` r
# Assuming get_env_identifiers() returns c("DrugName", "Gnumber") for drug identifiers
dt <- data.table::data.table(
  DrugName = c("DrugA", "DrugA", "DrugB"),
  Gnumber = c("G1", "G1", "G2"),
  Var1 = c(NA, "X", NA),
  Var2 = c(NA, "Y", "Z")
)
additional_vars <- c("Var1", "Var2")
 dt_updated <- update_drug_name(dt, additional_vars)
# Would update DrugName and Gnumber
dt_updated
#>                       DrugName                  Gnumber   Var1   Var2
#>                         <char>                   <char> <char> <char>
#> 1:                       DrugA                       G1   <NA>   <NA>
#> 2: DrugA (Var1 = X) (Var2 = Y) G1 (Var1 = X) (Var2 = Y)      X      Y
#> 3:            DrugB (Var2 = Z)            G2 (Var2 = Z)   <NA>      Z
```
