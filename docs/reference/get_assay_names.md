# get assay names of the given se/dataset fetch the data from the se if provided as metadata use predefined values from `get_env_assay_names` otherwise

get assay names of the given se/dataset fetch the data from the se if
provided as metadata use predefined values from `get_env_assay_names`
otherwise

## Usage

``` r
get_assay_names(se = NULL, ...)
```

## Arguments

- se:

  SummarizedExperiment or NULL

- ...:

  Additional arguments to pass to `get_env_assay_names`.

## Value

charvec

## Author

Arkadiusz Gładki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
get_assay_names()
#>            raw        control     normalized       averaged        metrics 
#>   "rawTreated"     "Controls"   "Normalized"     "Averaged"      "Metrics" 
#>         excess         scores   isobolograms 
#>       "excess"       "scores" "isobolograms" 
```
