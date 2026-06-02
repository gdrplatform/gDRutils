# Set fit parameters for an invalid fit.

Set fit parameters for an invalid fit.

## Usage

``` r
.set_invalid_fit_params(out, norm_values)
```

## Arguments

- out:

  Named list of fit parameters.

- norm_values:

  Numeric vector used to estimate an `xc50` value.

## Value

Modified named list of fit parameters.

## Examples

``` r
.set_invalid_fit_params(list(), norm_values = rep(0.3, 6))
#> $fit_type
#> [1] "DRCInvalidFitResult"
#> 
#> $r2
#> [1] NA
#> 
#> $xc50
#> [1] -Inf
#> 
```
