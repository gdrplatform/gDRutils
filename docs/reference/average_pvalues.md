# Average p-values using Fisher's method Combines a vector of p-values into a single representative p-value. It implements Fisher's method, where the test statistic is calculated as \$\$X\_{2k}^2 = -2 \sum\_{i=1}^{k} \ln(p_i)\$\$. This statistic follows a chi-squared distribution with 2k degrees of freedom (where k is the number of p-values), from which the combined p-value is derived.

Average p-values using Fisher's method Combines a vector of p-values
into a single representative p-value. It implements Fisher's method,
where the test statistic is calculated as \$\$X\_{2k}^2 = -2
\sum\_{i=1}^{k} \ln(p_i)\$\$. This statistic follows a chi-squared
distribution with 2k degrees of freedom (where k is the number of
p-values), from which the combined p-value is derived.

## Usage

``` r
average_pvalues(p_values)
```

## Arguments

- p_values:

  A numeric vector of p-values. Values are expected to be between 0
  and 1. The function assumes at least one non-NA value is provided.

## Value

A single, combined p-value as a numeric value.
