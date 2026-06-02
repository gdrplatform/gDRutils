# add arbitrary S3 class to an object

Modify and object's `class` attribute.

## Usage

``` r
addClass(x, newClass)
```

## Arguments

- x:

  an object

- newClass:

  character string; class to be added

## Value

The same object with an added S3 class.

## Details

This is a simple convenience function that an item to the `class`
attribute of an object so that it can be dispatched to a proper S3
method. This is purely for code clarity, so that individual methods do
not clutter the definitions of higher order functions.

## Examples

``` r
addClass(data.table::data.table(), "someClass")
#> Null data.table (0 rows and 0 cols)
```
