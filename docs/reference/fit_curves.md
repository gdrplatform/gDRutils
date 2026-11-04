# Fit curves

Fit GR and RV curves from a data.table.

## Usage

``` r
fit_curves(
  df_,
  series_identifiers,
  e_0 = 1,
  GR_0 = 1,
  n_point_cutoff = 4,
  range_conc = c(0.005, 5),
  force_fit = FALSE,
  pcutoff = 0.05,
  cap = 0.1,
  normalization_type = c("GR", "RV")
)
```

## Arguments

- df\_:

  data.table containing data to fit. See details.

- series_identifiers:

  character vector of the column names in `data.table` whose combination
  represents a unique series for which to fit curves.

- e_0:

  numeric value representing the `x_0` value for the RV curve. Defaults
  to `1`.

- GR_0:

  numeric value representing the `x_0` value for the GR curve. Defaults
  to `1`.

- n_point_cutoff:

  integer of how many points should be considered the minimum required
  to try to fit a curve. Defaults to `4`.

- range_conc:

  numeric vector of length 2 indicating the lower and upper
  concentration ranges. Defaults to `c(5e-3, 5)`. See details.

- force_fit:

  boolean indicating whether or not to force a constant fit. Defaults to
  `FALSE`.

- pcutoff:

  numeric of pvalue significance threshold above or equal to which to
  use a constant fit. Defaults to `0.05`.

- cap:

  numeric value capping `norm_values` to stay below (`x_0` + cap).
  Defaults to `0.1`.

- normalization_type:

  character vector of types of curves to fit. Defaults to
  `c("GR", "RV")`.

## Value

data.table of fit parameters as specified by the `normalization_type`.

## Details

The `df_` expects the following columns:

- RelativeViability normalized relative viability values (if
  `normalization_type` includes `"RV"`)

- GRvalue normalized GR values (if `normalization_type` includes `"GR"`)

The `range_conc` is used to calculate the `x_AOC_range` statistic. The
purpose of this statistic is to enable comparison across different
experiments with slightly different concentration ranges.

## Examples

``` r
df_ <- data.table::data.table(Concentration = c(0.001, 0.00316227766016838,
0.01, 0.0316227766016838),
x_std = c(0.1, 0.1, 0.1, 0.1), normalization_types = c("RV", "RV", "RV", "RV"),
x = c(0.9999964000144, 0.999964001439942, 0.999640143942423, 0.996414342629482))

fit_curves(df_, "Concentration", normalization_type = "RV")
#>       x_mean        x_AOC x_AOC_range  xc50     x_max      ec50     x_inf   x_0
#>        <num>        <num>       <num> <num>     <num>     <num>     <num> <num>
#> 1: 0.9994663 0.0005336782  0.07214082   Inf 0.9964143 0.2030077 0.8446057     1
#>           h    r2      p_value          rss maxlog10Concentration N_conc
#>       <num> <num>        <num>        <num>                 <num>  <int>
#> 1: 2.014503     1 3.799064e-08 3.425792e-13                  -1.5      4
#>    x_sd_avg               fit_type normalization_type fit_source
#>       <num>                 <char>             <char>     <char>
#> 1:      0.1 DRC3pHillFitModelFixS0                 RV        gDR
```
