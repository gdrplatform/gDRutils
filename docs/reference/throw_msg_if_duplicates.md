# throw message if assay data.table contains duplicated rows

An auxiliary function that checks for duplicated rows in assay
data.table, In case of duplicates it throws a message. The messsage
function is by default [`stop()`](https://rdrr.io/r/base/stop.html) The
message function can be customized with `msg_f` parameter

## Usage

``` r
throw_msg_if_duplicates(
  dt,
  assay_name = "unknown",
  msg_f = stop,
  preview_max_numb = 4
)
```

## Arguments

- dt:

  data.table with assay data

- assay_name:

  string with the name of the assay

- msg_f:

  function to be used to throw the message

- preview_max_numb:

  number of rows to preview if duplicates found

## Examples

``` r
sdata <- get_synthetic_data("finalMAE_small.qs2")
smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
throw_msg_if_duplicates(smetrics_data, assay_name = "Metrics", msg_f = futile.logger::flog.info)
```
