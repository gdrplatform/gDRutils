# Cap XC50 value.

Set IC50/GR50 value to `Inf` or `-Inf` based on upper and lower limits.

## Usage

``` r
cap_xc50(xc50, max_conc, min_conc = NA, capping_fold = 5)
```

## Arguments

- xc50:

  Numeric value of the IC50/GR50 to cap.

- max_conc:

  Numeric value of the highest concentration in a dose series used to
  calculate the `xc50`.

- min_conc:

  Numeric value of the lowest concentration in a dose series used to
  calculate the `xc50`. If `NA` (default), using `max_conc/1e5` instead.

- capping_fold:

  Integer value of the fold number to use for capping. Defaults to `5`.

## Value

Capped IC50/GR50 value.

## Details

Note: `xc50` and `max_conc` should share the same units. Ideally, the
`lower_cap` should be based on the lowest tested concentration. However,
since we don't record that, it is set 5 orders of magnitude below the
highest dose.

## Examples

``` r
cap_xc50(xc50 = 1, max_conc = 2)
#> [1] 1
cap_xc50(xc50 = 2, max_conc = 5, min_conc = 1)
#> [1] 2
cap_xc50(xc50 = 26, max_conc = 5, capping_fold = 5)
#> [1] Inf
```
