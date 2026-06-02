# Has Single Codrug Data

Has Single Codrug Data

## Usage

``` r
has_single_codrug_data(
  cols,
  prettify_identifiers = TRUE,
  codrug_identifiers = c("drug_name2", "concentration2")
)
```

## Arguments

- cols:

  character vector with the columns of the input data

- prettify_identifiers:

  logical flag specifying if identifiers are expected to be prettified

- codrug_identifiers:

  character vector with identifiers for the codrug columns

## Value

logical flag

## Examples

``` r
has_single_codrug_data("Drug Name")
#> [1] FALSE
has_single_codrug_data(c("Drug Name", "Cell Lines"))
#> [1] FALSE
has_single_codrug_data(c("Drug Name 2", "Concentration 2"))
#> [1] TRUE
has_single_codrug_data(
  get_prettified_identifiers(
    c("concentration2", "drug_name2"),
    simplify = FALSE
  )
)
#> [1] TRUE
```
