# Average biological replicates on the data table side.

Average biological replicates on the data table side.

## Usage

``` r
average_biological_replicates_dt(
  dt,
  var,
  prettified = FALSE,
  fixed = TRUE,
  geometric_average_fields = get_header("metric_average_fields")$geometric_mean,
  fit_type_average_fields = get_header("metric_average_fields")$fit_type,
  blacklisted_fields = get_header("metric_average_fields")$blacklisted,
  add_sd = FALSE,
  copy_dt = TRUE
)
```

## Arguments

- dt:

  data.table with Metric data

- var:

  String representing additional metadata of replicates

- prettified:

  Flag indicating if the provided identifiers in the dt are prettified

- fixed:

  Flag indicating whether to add a fix for -Inf in the geometric mean.

- geometric_average_fields:

  Character vector of column names in `dt` to take the geometric average
  of.

- fit_type_average_fields:

  Character vector of column names in `dt` that should be treated as a
  column with fit type data

- blacklisted_fields:

  Character vector of column names in `dt` that should be skipped in
  averaging

- add_sd:

  Flag indicating whether to add standard deviation and count columns.

- copy_dt:

  Flag indicating whether to copy the input data.table before modifying.
  Set to `FALSE` when the caller already provides a copy.

## Value

data.table without replicates

## Examples

``` r
dt <- data.table::data.table(a = c(seq_len(10), 1),
b = c(rep("drugA", 10), rep("drugB", 1)))
average_biological_replicates_dt(dt, var = "a")
#>         b
#>    <char>
#> 1:  drugA
#> 2:  drugB
```
