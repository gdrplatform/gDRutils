# Has Valid Codrug Data

Has Valid Codrug Data

## Usage

``` r
has_valid_codrug_data(
  data,
  prettify_identifiers = TRUE,
  codrug_name_identifier = "drug_name2",
  codrug_conc_identifier = "concentration2"
)
```

## Arguments

- data:

  data.table with input data

- prettify_identifiers:

  logical flag specifying if identifiers are expected to be prettified

- codrug_name_identifier:

  string with the identifiers for the codrug drug_name column

- codrug_conc_identifier:

  string with the identifiers for the codrug concentration column

## Value

logical flag

## Examples

``` r
dt <-
  data.table::data.table(
    "Drug Name" = letters[seq_len(3)],
    "Concentration" = seq_len(3),
    "Drug Name 2" = letters[4:6],
    "Concentration 2" = 4:6
  )
has_valid_codrug_data(dt)
#> [1] TRUE

dt$`Concentration 2` <- NULL
has_valid_codrug_data(dt)
#> [1] FALSE
```
