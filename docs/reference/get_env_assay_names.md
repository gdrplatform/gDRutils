# get default assay names for the specified filters, i.e. set of assay types, assay groups and assay data types

get default assay names for the specified filters, i.e. set of assay
types, assay groups and assay data types

## Usage

``` r
get_env_assay_names(
  type = NULL,
  group = NULL,
  data_type = NULL,
  prettify = FALSE,
  simplify = TRUE
)
```

## Arguments

- type:

  charvec of assay types

- group:

  charvec of assay groups

- data_type:

  charvec assay of data types

- prettify:

  logical flag, prettify the assay name?

- simplify:

  logical flag, simplify the output? will return single string instead
  of named vector with single element useful when function is expected
  to return single element/assay only

## Value

charvec

## Author

Arkadiusz Gładki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
get_env_assay_names()
#>            raw        control     normalized       averaged        metrics 
#>   "rawTreated"     "Controls"   "Normalized"     "Averaged"      "Metrics" 
#>         excess         scores   isobolograms 
#>       "excess"       "scores" "isobolograms" 
```
