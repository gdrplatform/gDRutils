# Validate dimnames

Assure that dimnames of two objects are the same

## Usage

``` r
validate_dimnames(obj, obj2, skip_empty = TRUE)
```

## Arguments

- obj:

  first object with dimnames to compare

- obj2:

  second object with dimnames to compare

- skip_empty:

  a logical indicating whether to skip comparison if a given dimname is
  empty in the case of both objects

## Value

`NULL`
