# Predict a concentration for a given efficacy with fit parameters.

Predict a concentration for a given efficacy with fit parameters.

## Usage

``` r
predict_conc_from_efficacy(efficacy, x_inf, x_0, ec50, h)
```

## Arguments

- efficacy:

  Numeric vector representing efficacies to predict concentrations for.

- x_inf:

  Numeric vector representing the asymptotic value of the sigmoidal fit
  to the dose-response data as concentration goes to infinity.

- x_0:

  Numeric vector representing the asymptotic metric value corresponding
  to a concentration of 0 for the primary drug.

- ec50:

  Numeric vector representing the drug concentration at half-maximal
  effect.

- h:

  Numeric vector representing the hill coefficient of the fitted curve,
  which reflects how steep

## Value

Numeric vector representing predicted concentrations from given
efficacies and fit parameters.

## Details

The inverse of this function is `predict_efficacy_from_conc`.

## See also

predict_efficacy_from_conc .calculate_x50

## Examples

``` r
predict_conc_from_efficacy(efficacy = c(1, 1.5), x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
#> [1] 0 0
```
