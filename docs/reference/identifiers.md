# Get, set, or reset identifiers for one or all identifier field(s)

Get, set, or reset the expected identifier(s) for one or all identifier
field(s). Identifiers are used by the gDR processing functions to
identify which columns in a `data.table` correspond to certain expected
fields. Functions of the family `*et_identifier` will look for
identifiers from the environment while functions of the family
`*et_SE_identifiers` will look for identifiers in the `metadata` slot of
a `SummarizedExperiment` object. See details for expected identifiers
and their definitions.

## Usage

``` r
get_env_identifiers(k = NULL, simplify = TRUE)

get_prettified_identifiers(k = NULL, simplify = TRUE)

set_env_identifier(k, v)

reset_env_identifiers()
```

## Arguments

- k:

  String corresponding to identifier name.

- simplify:

  Boolean indicating whether output should be simplified.

- v:

  Character vector corresponding to the value for given identifier `k`.

## Value

For any `set`ting or `reset`ting functionality, a `NULL` invisibly. For
`get_env_identifiers` a character vector of identifiers for field `k`.
For functions called with no arguments, the entire available identifier
list is returned.

list or charvec depends on unify param

list or charvec depends on unify param

`NULL`

`NULL`

## Details

Identifiers supported by the gDR suite include:

- "barcode": String of column name containing barcode metadata

- "cellline": String of column name containing unique, machine-readable
  cell line identifiers

- "cellline_name": String of column name containing human-friendly cell
  line names

- "cellline_tissue": String of column name containing metadata on cell
  line tissue type

- "cellline_ref_div_time": String of column name containing reference
  division time for cell lines

- "cellline_parental_identifier": String of column name containing
  unique, machine-readable parental cell line identifiers. Used in the
  case of derived or engineered cell lines.

- "drug": String of column name containing unique, machine-readable drug
  identifiers

- "drug_name": String of column name containing human-friendly drug
  names

- "drug_moa": String of column name containing metadata for drug mode of
  action

- "duration": String of column name containing metadata on duration that
  cells were treated (in hours)

- "template": String of column name containing template metadata

- "untreated_tag": Character vector of entries that identify control,
  untreated wells

- "well_position": Character vector of column names containing metadata
  on well positions on a plate

## Examples

``` r
get_env_identifiers("duration") # "Duration"
#> [1] "Duration"
```
