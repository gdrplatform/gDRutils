# get_experiment_groups

get experiment groups

## Usage

``` r
get_experiment_groups(type = NULL)
```

## Arguments

- type:

  String indicating the name of an assay group. Defaults to all
  experiment groups.

## Value

list with experiment groups or string (if type not NULL)

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
get_experiment_groups()
#> $`single-agent`
#>   single-agent    co-dilution 
#> "single-agent"  "co-dilution" 
#> 
#> $combination
#> [1] "combination"
#> 
#> $`time-course`
#> [1] "time-course"
#> 
```
