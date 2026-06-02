# assert choices

assert choices

## Usage

``` r
assert_choices(x, choices, ...)
```

## Arguments

- x:

  charvec expected subset

- choices:

  charvec reference set

- ...:

  Additional arguments to pass to
  [`checkmate::test_choice`](https://mllg.github.io/checkmate/reference/checkChoice.html)

## Value

`NULL`

## Examples

``` r
assert_choices("x", c("x","y"))
```
