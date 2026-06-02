# modify assay with additional data

Reduce dimensionality of an assay by dropping extra data or combining
variables.

## Usage

``` r
modifyData(x, ...)

# S3 method for class 'drug_name2'
modifyData(x, option, keep, ...)

# S3 method for class 'data_source'
modifyData(x, option, keep, ...)

# Default S3 method
modifyData(x, option, keep, ...)
```

## Arguments

- x:

  a `data.table` containing an assay

- ...:

  additional arguments passed to methods

- option:

  character string specifying the action to be taken, see `Details`

- keep:

  character string specifying the value of the active variable that will
  be kept

## Value

modified object

## Details

If an essay extracted from a `SummarizedExperiment` contains additional
information, i.e. factors beyond `DrugName` and `CellLineName`, that
information will be treated in one of three ways, depending on the value
of `option`:

- `drop`: Some information will be discarded and only one value of the
  additional variable (chosen by the user) will be kept.

- `toDrug`: The information is pasted together with the primary drug
  name. All observations are kept.

- `toCellLine`: As above, but pasting is done with cell line name.

Depending on the type of additional information, the exact details will
differ. This is handled in the app by adding special classes to the data
tables and dispatching to S3 methods.

Following modification, the additional columns are discarded.

## Methods (by class)

- `modifyData(drug_name2)`: includes the name and concentration of the
  second drug

- `modifyData(data_source)`: includes the data source

- `modifyData(default)`: includes the name of other additional variables

## Examples

``` r
dt <- data.table::data.table(a = as.character(1:10), b = "data")
dt <- addClass(dt, "a")
modifyData(dt, "average", keep = "b")
#>         b
#>    <char>
#> 1:   data
```
