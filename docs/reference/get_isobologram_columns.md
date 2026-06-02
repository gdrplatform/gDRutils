# Get isobologram column names

Get isobologram column names

## Usage

``` r
get_isobologram_columns(k = NULL, prettify = TRUE)
```

## Arguments

- k:

  key

- prettify:

  change to upper case and add underscore, iso_level –\> Iso_Level

## Value

character vector of isobologram column names for combination data

## Examples

``` r
get_isobologram_columns()
#> [1] "Iso_Level"        "Pos_x"            "Pos_x_Ref"        "Pos_y"           
#> [5] "Pos_y_Ref"        "Log10_Ratio_Conc" "Log2_CI"         
get_isobologram_columns("iso_level", prettify = TRUE)
#> [1] "Iso_Level"
```
