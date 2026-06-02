# Get descriptions for identifiers

Get descriptions for identifiers

## Usage

``` r
get_identifiers_dt(k = NULL, get_description = FALSE, get_example = FALSE)
```

## Arguments

- k:

  identifier key, string

- get_description:

  return descriptions only, boolean

- get_example:

  return examples only, boolean

## Value

named list

## Examples

``` r
get_identifiers_dt()
#> $duration
#> $duration$description
#> [1] "Duration of treatment (in h)"
#> 
#> $duration$example
#> [1] "120"
#> 
#> 
#> $clid
#> $clid$description
#> [1] "Cell Line identifier"
#> 
#> $clid$example
#> [1] "CL00011"
#> 
#> 
#> $cellline_name
#> $cellline_name$description
#> [1] "Cell Line common name"
#> 
#> $cellline_name$example
#> [1] "MCF-10A"
#> 
#> 
#> $cellline_tissue
#> $cellline_tissue$description
#> [1] "Tissue name"
#> 
#> $cellline_tissue$example
#> [1] "Lung"
#> 
#> 
#> $cellline_ref_div_time
#> $cellline_ref_div_time$description
#> [1] "Cell Line reference division time (in h)"
#> 
#> $cellline_ref_div_time$example
#> [1] "50"
#> 
#> 
#> $cellline_parental_identifier
#> $cellline_parental_identifier$description
#> [1] "Cell Line parental identifier"
#> 
#> $cellline_parental_identifier$example
#> [1] "MCF-10A"
#> 
#> 
#> $cellline_subtype
#> $cellline_subtype$description
#> [1] "Cell Line subtype"
#> 
#> $cellline_subtype$example
#> [1] "TNBC"
#> 
#> 
#> $drug
#> $drug$description
#> [1] "Drug identifier"
#> 
#> $drug$example
#> [1] "G012345"
#> 
#> 
#> $drug_name
#> $drug_name$description
#> [1] "Drug common name"
#> 
#> $drug_name$example
#> [1] "Paclitaxel"
#> 
#> 
#> $drug_moa
#> $drug_moa$description
#> [1] "Drug mechanism of action"
#> 
#> $drug_moa$example
#> [1] "Tubulin"
#> 
#> 
#> $untreated_tag
#> $untreated_tag$description
#> [1] "Tags for untreated cellines"
#> 
#> $untreated_tag$example
#> [1] "Untreated"
#> 
#> 
#> $masked_tag
#> $masked_tag$description
#> [1] "Tag for masked celllines"
#> 
#> $masked_tag$example
#> [1] TRUE
#> 
#> 
#> $well_position
#> $well_position$description
#> [1] "Position in the wells"
#> 
#> $well_position$example
#> [1] "[A, 1]"
#> 
#> 
#> $concentration
#> $concentration$description
#> [1] "Concentration of the drug used (uM)"
#> 
#> $concentration$exammple
#> [1] "3.333"
#> 
#> 
#> $template
#> $template$description
#> [1] "Treatment layout template"
#> 
#> $template$example
#> [1] "Template_untreated.xlsx"
#> 
#> 
#> $barcode
#> $barcode$description
#> [1] "Plate identifier"
#> 
#> $barcode$example
#> [1] "1234C65465876_001"
#> 
#> 
#> $drug_2
#> $drug_2$description
#> [1] "Drug 2 identifier"
#> 
#> $drug_2$example
#> [1] "G123456"
#> 
#> 
#> $drug_name2
#> $drug_name2$description
#> [1] "Drug 2 common name"
#> 
#> $drug_name2$example
#> [1] "Erlotinib"
#> 
#> 
#> $drug_moa2
#> $drug_moa2$description
#> [1] "Drug 2 mechanism of action"
#> 
#> $drug_moa2$example
#> [1] "EGFR"
#> 
#> 
#> $concentration2
#> $concentration2$description
#> [1] "Concentration of Drug 2 used (uM)"
#> 
#> $concentration2$example
#> [1] "1.111"
#> 
#> 
#> $drug_3
#> $drug_3$description
#> [1] "Drug 3 identifier"
#> 
#> $drug_3$example
#> [1] "G345678"
#> 
#> 
#> $drug_name3
#> $drug_name3$description
#> [1] "Drug 3 common name"
#> 
#> $drug_name3$example
#> [1] "Lapatinib"
#> 
#> 
#> $drug_moa3
#> $drug_moa3$description
#> [1] "Drug 3 mechanism of action"
#> 
#> $drug_moa3$example
#> [1] "HER2"
#> 
#> 
#> $concentration3
#> $concentration3$description
#> [1] "Concentration of Drug 3 used"
#> 
#> $concentration3$example
#> [1] "10.0"
#> 
#> 
#> $data_source
#> $data_source$description
#> [1] "Source of the data"
#> 
#> $data_source$example
#> [1] "gDR"
#> 
#> 
```
