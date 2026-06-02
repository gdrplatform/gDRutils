# Set fit parameters for a constant fit.

Replace values for flat fits: ec50 = 0, h = 0.0001 and xc50 = +/- Inf

## Usage

``` r
set_constant_fit_params(out, mean_norm_value)
```

## Arguments

- out:

  Named list of fit parameters.

- mean_norm_value:

  Numeric value that be used to set all parameters that can be
  calculated from the mean.

## Value

Modified named list of fit parameters.

## Examples

``` r
na <- list(x_0 = NA)
set_constant_fit_params(na, mean_norm_value = 0.6)
#> $x_0
#> [1] 0.6
#> 
#> $fit_type
#> [1] "DRCConstantFitResult"
#> 
#> $ec50
#> [1] 0
#> 
#> $h
#> [1] 1e-04
#> 
#> $r2
#> [1] 0
#> 
#> $xc50
#> [1] Inf
#> 
#> $x_mean
#> [1] 0.6
#> 
#> $x_inf
#> [1] 0.6
#> 
#> $x_AOC
#> [1] 0.4
#> 
#> $x_AOC_range
#> [1] 0.4
#> 
```
