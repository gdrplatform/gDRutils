# Identify and return additional variables in list of dt

Identify and return additional variables in list of dt

## Usage

``` r
get_additional_variables(dt_list, unique = FALSE, prettified = FALSE)
```

## Arguments

- dt_list:

  list of data.table or data.table containing additional variables

- unique:

  logical flag indicating if all variables should be returned or only
  those containing more than one unique value

- prettified:

  Flag indicating if the provided identifiers in the dt are prettified

## Value

vector of variable names with additional variables

## Examples

``` r
dt <- data.table::data.table(
  Gnumber = seq_len(10),
  Concentration = runif(10),
  Ligand = c(rep(0.5, 5), rep(0, 5))
)
get_additional_variables(dt)
#> [1] "Ligand"
```
