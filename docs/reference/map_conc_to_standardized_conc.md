# Create a mapping of concentrations to standardized concentrations.

Create a mapping of concentrations to standardized concentrations.

## Usage

``` r
map_conc_to_standardized_conc(conc1, conc2)
```

## Arguments

- conc1:

  numeric vector of the concentrations for drug 1.

- conc2:

  numeric vector of the concentrations for drug 2.

## Value

data.table of 2 columns named `"concs"` and `"rconcs"` containing the
original concentrations and their closest matched standardized
concentrations respectively. and their new standardized concentrations.

## Details

The concentrations are standardized in that they will contain regularly
spaced dilutions and close values will be rounded.

## Examples

``` r
ratio <- 0.5
conc1 <- c(0, 10 ^ (seq(-3, 1, ratio)))

shorter_range <- conc1[-1]
noise <- runif(length(shorter_range), 1e-12, 1e-11)
conc2 <- shorter_range + noise

map_conc_to_standardized_conc(conc1, conc2)
#>            concs   rconcs
#>            <num>    <num>
#>  1:  0.000000000  0.00000
#>  2:  0.001000000  0.00100
#>  3:  0.003162278  0.00316
#>  4:  0.010000000  0.01000
#>  5:  0.031622777  0.03160
#>  6:  0.100000000  0.10000
#>  7:  0.316227766  0.31600
#>  8:  1.000000000  1.00000
#>  9:  3.162277660  3.16000
#> 10: 10.000000000 10.00000
#> 11:  0.001000000  0.00100
#> 12:  0.003162278  0.00316
#> 13:  0.010000000  0.01000
#> 14:  0.031622777  0.03160
#> 15:  0.100000000  0.10000
#> 16:  0.316227766  0.31600
#> 17:  1.000000000  1.00000
#> 18:  3.162277660  3.16000
#> 19: 10.000000000 10.00000
```
