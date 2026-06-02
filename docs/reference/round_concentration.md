# Round concentration to ndigit significant digits

Round concentration to ndigit significant digits

## Usage

``` r
round_concentration(x, ndigit = 3)
```

## Arguments

- x:

  value to be rounded.

- ndigit:

  number of significant digits (default = 4).

## Value

rounded x

## Examples

``` r
round_concentration(x = c(0.00175,0.00324,0.0091), ndigit = 1)
#> [1] 0.002 0.003 0.010
```
