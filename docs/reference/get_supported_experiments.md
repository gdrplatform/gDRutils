# get_supported_experiments

get supported experiments

## Usage

``` r
get_supported_experiments(type = NULL)
```

## Arguments

- type:

  String indicating the type of experiment

## Value

charvec with supported experiment name(s)

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
get_supported_experiments()
#> [1] "single-agent" "combination"  "co-dilution"  "time-course" 
```
