# Standardize concentration values.

Standardize concentration values.

## Usage

``` r
.standardize_conc(conc)
```

## Arguments

- conc:

  numeric vector of the concentrations

## Value

vector of standardized concentrations

## Details

If no `conc` are passed, `NULL` is returned.

## Examples

``` r
concs <- 10 ^ (seq(-1, 1, 0.9))
.standardize_conc(concs)
#> [1] 0.100 0.794 6.310
```
