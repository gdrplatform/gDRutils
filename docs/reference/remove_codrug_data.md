# Remove Codrug Data

Remove Codrug Data

## Usage

``` r
remove_codrug_data(
  data,
  prettify_identifiers = TRUE,
  codrug_identifiers = c("drug_name2", "concentration2")
)
```

## Arguments

- data:

  data.table with input data

- prettify_identifiers:

  logical flag specifying if identifiers are expected to be prettified

- codrug_identifiers:

  character vector with identifiers for the codrug columns

## Value

data.table without combination columns

## Examples

``` r
dt <-
  data.table::data.table(
    "Drug Name" = letters[seq_len(3)],
    "Concentration" = seq_len(3),
    "Drug Name 2" = letters[4:6],
    "Concentration 2" = 4:6
  )
dt
#>    Drug Name Concentration Drug Name 2 Concentration 2
#>       <char>         <int>      <char>           <int>
#> 1:         a             1           d               4
#> 2:         b             2           e               5
#> 3:         c             3           f               6
remove_codrug_data(dt)
#>    Drug Name Concentration
#>       <char>         <int>
#> 1:         a             1
#> 2:         b             2
#> 3:         c             3
```
