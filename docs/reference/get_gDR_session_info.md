# get gDR package and their version installed in the environment

get gDR package and their version installed in the environment

## Usage

``` r
get_gDR_session_info(pattern = "^gDR")
```

## Arguments

- pattern:

  string with the pattern to grep R packages from the list of installed
  packages

## Value

data.table with gDR packages and their versions

## Examples

``` r
get_gDR_session_info()
#>        Package Version
#>         <char>  <char>
#> 1:    gDRutils  1.11.3
#> 2:    gDRstyle  1.10.0
#> 3: gDRtestData  1.10.0
```
