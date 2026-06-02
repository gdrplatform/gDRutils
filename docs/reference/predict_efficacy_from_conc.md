# Predict efficacy values given fit parameters and a concentration.

Predict efficacy values given fit parameters and a concentration.

## Usage

``` r
predict_efficacy_from_conc(c, x_inf, x_0, ec50, h)
```

## Arguments

- c:

  Numeric vector representing concentrations to predict efficacies for.

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
  which reflects how steep the dose-response curve is.

## Value

Numeric vector representing predicted efficacies from given
concentrations and fit parameters.

## Details

The inverse of this function is `predict_conc_from_efficacy`.

## See also

predict_conc_from_efficacy

## Examples

``` r
predict_efficacy_from_conc(c = 1, x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
#> [1] 0.28
```
