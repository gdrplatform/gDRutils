# Convert colData to JSON

Convert colData to JSON format for elasticsearch indexing.

## Usage

``` r
convert_colData_to_json(
  cdata,
  identifiers,
  req_cols = c("cellline", "cellline_name", "cellline_tissue", "cellline_ref_div_time")
)
```

## Arguments

- cdata:

  data.table of `colData`.

- identifiers:

  charvec with identifiers

- req_cols:

  charvec required columns

## Value

JSON string capturing the `cdata`.

## Details

Standardizes the `cdata` to common schema fields and tidies formatting
to be condusive to joining with other JSON responses.

## Examples

``` r
cdata <- data.table::data.table(
  mycellline = letters,
  mycelllinename = letters,
  mycelllinetissue = letters,
  cellline_ref_div_time = "cellline_ref_div_time")
identifiers <- list(cellline = "mycellline",
                    cellline_name = "mycelllinename",
                    cellline_ref_div_time = "cellline_ref_div_time",
                    cellline_tissue = "mycelllinetissue")
convert_colData_to_json(cdata, identifiers)
#> [1] "\"cellline\":[\"a\",\"b\",\"c\",\"d\",\"e\",\"f\",\"g\",\"h\",\"i\",\"j\",\"k\",\"l\",\"m\",\"n\",\"o\",\"p\",\"q\",\"r\",\"s\",\"t\",\"u\",\"v\",\"w\",\"x\",\"y\",\"z\"],\"cellline_name\":[\"a\",\"b\",\"c\",\"d\",\"e\",\"f\",\"g\",\"h\",\"i\",\"j\",\"k\",\"l\",\"m\",\"n\",\"o\",\"p\",\"q\",\"r\",\"s\",\"t\",\"u\",\"v\",\"w\",\"x\",\"y\",\"z\"],\"cellline_tissue\":[\"a\",\"b\",\"c\",\"d\",\"e\",\"f\",\"g\",\"h\",\"i\",\"j\",\"k\",\"l\",\"m\",\"n\",\"o\",\"p\",\"q\",\"r\",\"s\",\"t\",\"u\",\"v\",\"w\",\"x\",\"y\",\"z\"],\"cellline_ref_div_time\":[\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\",\"cellline_ref_div_time\"], \"misc_coldata\": {}"
```
