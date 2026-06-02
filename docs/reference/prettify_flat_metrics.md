# Prettify metric names in flat 'Metrics' assay

Map existing column names of a flattened 'Metrics' assay to prettified
names.

## Usage

``` r
prettify_flat_metrics(
  x,
  human_readable = FALSE,
  normalization_type = c("GR", "RV")
)
```

## Arguments

- x:

  character vector of names to prettify.

- human_readable:

  boolean indicating whether or not to return column names in human
  readable format. Defaults to `FALSE`.

- normalization_type:

  character vector with a specified normalization type. Defaults to
  `c("GR", "RV")`.

## Value

character vector of prettified names.

## Details

A common use case for this function is to prettify column names from a
flattened version of the `"Metrics"` assay. Mode
`"human_readable" = TRUE` is often used for prettification in the
context of front-end applications, whereas `"human_readable" = FALSE` is
often used for prettification in the context of the command line.

## Examples

``` r
x <- c("CellLineName", "Tissue", "Primary Tissue", "GR_gDR_x_mean", "GR_gDR_xc50", "RV_GDS_x_mean")
prettify_flat_metrics(x, human_readable = FALSE)
#> [1] "CellLineName"   "Tissue"         "Primary Tissue" "GR_mean"       
#> [5] "GR50"           "GDS_RV_mean"   
```
