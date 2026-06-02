# check if assay data contains duplicated data

An auxiliary function that checks for duplicates in the assay data

## Usage

``` r
has_assay_dt_duplicated_rows(dt)
```

## Arguments

- dt:

  data.table with assay data

## Value

logical flag indicating if a dt contains duplicated rows or not

## Examples

``` r
sdata <- get_synthetic_data("finalMAE_small.qs2")
smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
has_assay_dt_duplicated_rows(smetrics_data)
#> [1] FALSE
```
