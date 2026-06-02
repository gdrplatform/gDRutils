# get names of combo assays

get names of combo assays

## Usage

``` r
get_combo_assay_names(se = NULL, ...)
```

## Arguments

- se:

  SummarizedExperiment or NULL

- ...:

  Additional arguments to pass to `get_assay_names`.

## Value

charvec of combo assay names.

## Author

Arkadiusz Gładki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
get_combo_assay_names()
#>         excess         scores   isobolograms 
#>       "excess"       "scores" "isobolograms" 
```
