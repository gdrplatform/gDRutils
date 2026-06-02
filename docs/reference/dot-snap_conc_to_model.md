# Snap a concentration to the nearest available model concentration

Finds the value in a vector of available concentrations that is closest
(on a log scale) to a user-specified concentration. This is an internal
helper function.

## Usage

``` r
.snap_conc_to_model(user_conc, available_concs)
```

## Arguments

- user_conc:

  A single numeric value for the desired concentration.

- available_concs:

  A numeric vector of concentrations for which a model exists.

## Value

A single numeric value from 'available_concs'.
