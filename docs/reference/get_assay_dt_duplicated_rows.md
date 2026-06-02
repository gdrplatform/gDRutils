# Helper function to find duplicated rows in assay data

Helper function to find duplicated rows in assay data

## Usage

``` r
get_assay_dt_duplicated_rows(dt, output = "index")
```

## Arguments

- dt:

  data.table

- output:

  string with the output format to be returned

## Value

integer vector or data.table with duplicated rows

## Examples

``` r
sdata <- get_synthetic_data("finalMAE_small.qs2")
smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
get_assay_dt_duplicated_rows(smetrics_data, output = "data")
#> Empty data.table (0 rows and 28 cols): rId,cId,x_mean,x_AOC,x_AOC_range,xc50...
get_assay_dt_duplicated_rows(smetrics_data)
#> integer(0)
```
