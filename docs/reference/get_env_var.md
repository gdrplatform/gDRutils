# safe wrapper of Sys.getenv()

So far the helper is needed to handle env vars containing `:` for which
the backslash is automatically added in some contexts and R could not
get the original value for these env vars.

## Usage

``` r
get_env_var(x, ...)
```

## Arguments

- x:

  string with the name of the environmental variable

- ...:

  additional params for Sys.getenev

## Value

sanitized value of the env variable

## Examples

``` r
get_env_var("HOME")
#> [1] "/home/runner"
```
