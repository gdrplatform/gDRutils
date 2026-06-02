# get combo assay names based on the field name

get combo assay names based on the field name

## Usage

``` r
convert_combo_field_to_assay(field)
```

## Arguments

- field:

  String containing name of the field for which the assay name should be
  returned

## Value

charvec

## Examples

``` r
convert_combo_field_to_assay("hsa_score")
#> [1] "scores"
```
