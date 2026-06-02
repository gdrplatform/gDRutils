# Check that specified identifier values exist in the data.

Check that specified identifier values exist in the data and error
otherwise.

## Usage

``` r
validate_identifiers(
  df,
  identifiers = NULL,
  req_ids = NULL,
  exp_one_ids = NULL
)
```

## Arguments

- df:

  data.table with `colnames`.

- identifiers:

  Named list of identifiers where the `names` are standardized
  identifier names. If not passed, defaults to
  [`get_env_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.md).

- req_ids:

  Character vector of standardized identifier names required to pass
  identifier validation.

- exp_one_ids:

  Character vector of standardized identifiers names where only one
  identifier value is expected. If not passed, defaults to
  [`get_expect_one_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/get_expect_one_identifiers.md).

## Value

Named list of identifiers modified to pass validation against the input
data. Errors with explanatory message if validation cannot pass with the
given identifiers and data.

## Details

Note that this does NOT set the identifiers anywhere (i.e. environment
or `SummarizedExperiment` object). If identifiers do not validate, will
throw error as side effect.

## Examples

``` r
validate_identifiers(
  S4Vectors::DataFrame("Barcode" = NA, "Duration" = NA, "Template" = NA, "clid" = NA),
  req_ids = "barcode"
)
#> $duration
#> [1] "Duration"
#> 
#> $cellline
#> [1] "clid"
#> 
#> $cellline_name
#> [1] "CellLineName"
#> 
#> $cellline_tissue
#> [1] "Tissue"
#> 
#> $cellline_ref_div_time
#> [1] "ReferenceDivisionTime"
#> 
#> $cellline_parental_identifier
#> [1] "parental_identifier"
#> 
#> $cellline_subtype
#> [1] "subtype"
#> 
#> $drug
#> [1] "Gnumber"
#> 
#> $drug_name
#> [1] "DrugName"
#> 
#> $drug_moa
#> [1] "drug_moa"
#> 
#> $untreated_tag
#> [1] "vehicle"   "untreated"
#> 
#> $masked_tag
#> [1] "masked"
#> 
#> $well_position
#> [1] "WellRow"    "WellColumn"
#> 
#> $concentration
#> [1] "Concentration"
#> 
#> $template
#> [1] "Template"
#> 
#> $barcode
#> [1] "Barcode"
#> 
#> $drug2
#> [1] "Gnumber_2"
#> 
#> $drug_name2
#> [1] "DrugName_2"
#> 
#> $drug_moa2
#> [1] "drug_moa_2"
#> 
#> $concentration2
#> [1] "Concentration_2"
#> 
#> $drug3
#> [1] "Gnumber_3"
#> 
#> $drug_name3
#> [1] "DrugName_3"
#> 
#> $drug_moa3
#> [1] "drug_moa_3"
#> 
#> $concentration3
#> [1] "Concentration_3"
#> 
#> $data_source
#> [1] "data_source"
#> 
#> $replicate
#> [1] "Replicate"
#> 
#> $normalization_type
#> [1] "normalization_type"
#> 
```
