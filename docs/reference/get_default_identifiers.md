# Get gDR default identifiers required for downstream analysis.

Get gDR default identifiers required for downstream analysis.

## Usage

``` r
get_default_identifiers()
```

## Value

charvec

## Examples

``` r
get_default_identifiers()
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
#> [1] "Template"  "Treatment"
#> 
#> $barcode
#> [1] "Barcode" "Plate"  
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
