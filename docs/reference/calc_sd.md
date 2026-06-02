# Calculate Standard Deviation or Return Zero

This function calculates the standard deviation of a numeric vector. If
the vector has a length of 1 and it is numeric, it returns 0.

## Usage

``` r
calc_sd(x)
```

## Arguments

- x:

  A numeric vector.

## Value

The standard deviation of the vector if its length is greater than 1 or
it is not numeric, otherwise 0.

## Examples

``` r
calc_sd(c(1, 2, 3, 4, 5)) # Should return the standard deviation
#> [1] 1.581139
calc_sd(c(1)) # Should return 0
#> [1] 0
calc_sd(numeric(0)) # Should return NA
#> [1] NA
calc_sd(c("a", "b", "c")) # Should return NA
#> Warning: NAs introduced by coercion
#> [1] NA
```
