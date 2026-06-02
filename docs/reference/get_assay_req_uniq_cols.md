# get columns in the assay data required to have unique data

get columns in the assay data required to have unique (non-duplicated)
data

## Usage

``` r
get_assay_req_uniq_cols(dt)
```

## Arguments

- dt:

  data.table with assay data

## Value

charvec with columns required to have unique data

## Examples

``` r
sdata <- get_synthetic_data("finalMAE_small.qs2")
smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
get_assay_req_uniq_cols(smetrics_data)
#> [1] "Gnumber"            "Duration"           "DrugName"          
#> [4] "clid"               "CellLineName"       "normalization_type"
```
