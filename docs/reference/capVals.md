# Cap metric values

Convenience function to apply caps to outlying metric values.

## Usage

``` r
capVals(x)
```

## Arguments

- x:

  `data.table` containing growth metrics extracted from a
  `SummarizedExperiment`

## Value

A data table with capped values.

## Details

The following metrics are capped at the respective values:

- `E max`: 0 - 1.1

- `GR max`: -1 - 1.1

- `RV AOC within set range`: over -0.1

- `GR AOC within set range`: over of -0.1

- `GR50`: 1e-4 to 30

- `IC50`: 1e-4 to 30

- `EC50`: 1e-4 to 30 (change 0 to NA beforehand)

## See also

`convert_se_assay_to_dt`,
[`oob`](https://scales.r-lib.org/reference/oob.html)

## Examples

``` r
dt <- data.table::data.table(
  `E Max` = c(-0.1, 0, 0.5, 1.2),
  `GR Max` = c(-1.1, -1, 0.5, 1.2),
  `RV AOC within set range` = c(-0.2, -0.1, 0, 3),
  `GR AOC within set range` = c(-0.2, -0.1, 0, 3),
  `GR50` = c(0, 1e-7, 10, 34),
  `IC50` = c(0, 1e-7, 10, 34),
  `EC50` = c(0, 1e-7, 10, 34),
  check.names = FALSE
)
dt
#>    E Max GR Max RV AOC within set range GR AOC within set range    GR50    IC50
#>    <num>  <num>                   <num>                   <num>   <num>   <num>
#> 1:  -0.1   -1.1                    -0.2                    -0.2 0.0e+00 0.0e+00
#> 2:   0.0   -1.0                    -0.1                    -0.1 1.0e-07 1.0e-07
#> 3:   0.5    0.5                     0.0                     0.0 1.0e+01 1.0e+01
#> 4:   1.2    1.2                     3.0                     3.0 3.4e+01 3.4e+01
#>       EC50
#>      <num>
#> 1: 0.0e+00
#> 2: 1.0e-07
#> 3: 1.0e+01
#> 4: 3.4e+01
dt1 <- capVals(dt)
dt1
#>    E Max GR Max RV AOC within set range GR AOC within set range  GR50  IC50
#>    <num>  <num>                   <num>                   <num> <num> <num>
#> 1:   0.0   -1.0                    -0.1                    -0.1 1e-04 1e-04
#> 2:   0.0   -1.0                    -0.1                    -0.1 1e-04 1e-04
#> 3:   0.5    0.5                     0.0                     0.0 1e+01 1e+01
#> 4:   1.1    1.1                     3.0                     3.0 3e+01 3e+01
#>     EC50
#>    <num>
#> 1:    NA
#> 2: 1e-04
#> 3: 1e+01
#> 4: 3e+01
```
