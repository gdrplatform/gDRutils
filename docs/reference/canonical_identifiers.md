# Canonical form of a set of identifiers, and whether two sets agree

`canonical_identifiers()` completes `ids` with the session's value of
every key it lacks (`fill`), recording those keys in
`attr(, "assumed")`: a key an object does not store is read from the
session, so it is compared as the session's value too. Keys that are not
identifiers (e.g. from an older object) are left out and named in
`attr(, "unknown")`. Given `data`, it also resolves an identifier that
lists alternative columns (e.g. `barcode = c("Barcode", "Plate")`) to
the one present, so two objects built the same way store the same
identifiers.

## Usage

``` r
canonical_identifiers(ids, data = NULL, fill = get_env_identifiers())

same_identifiers(a, b)
```

## Arguments

- ids, a, b:

  named lists of identifiers.

- data:

  optional `data.frame`, or character vector of column names, the
  identifiers describe.

- fill:

  named list giving the value of a key `ids` lacks. Defaults to the
  session's identifiers.

## Value

`canonical_identifiers()`: a named list over every identifier key.
`same_identifiers()`: a flag.

## Details

`same_identifiers()` is `TRUE` when two sets name the same columns. A
key resolved to one column agrees with a list of alternatives that
contains it, and lists of alternatives or tags (`barcode`, `template`,
`untreated_tag`) are compared as sets; `well_position` keeps its order.
The keys that differ are in `attr(, "differ")`.

## Examples

``` r
ids <- canonical_identifiers(list(barcode = c("Barcode", "Plate")), data = "Plate")
ids$barcode # "Plate"
#> [1] "Plate"
attr(ids, "assumed")
#>  [1] "duration"                     "cellline"                    
#>  [3] "cellline_name"                "cellline_tissue"             
#>  [5] "cellline_ref_div_time"        "cellline_parental_identifier"
#>  [7] "cellline_subtype"             "drug"                        
#>  [9] "drug_name"                    "drug_moa"                    
#> [11] "untreated_tag"                "masked_tag"                  
#> [13] "well_position"                "concentration"               
#> [15] "template"                     "drug2"                       
#> [17] "drug_name2"                   "drug_moa2"                   
#> [19] "concentration2"               "drug3"                       
#> [21] "drug_name3"                   "drug_moa3"                   
#> [23] "concentration3"               "data_source"                 
#> [25] "replicate"                    "normalization_type"          
same_identifiers(list(barcode = "Plate"), get_default_identifiers())
#> [1] TRUE
#> attr(,"differ")
#> character(0)
```
