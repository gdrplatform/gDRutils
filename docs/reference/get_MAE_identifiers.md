# get_MAE_identifiers

get the identifiers of all SE's in the MAE

## Usage

``` r
get_MAE_identifiers(mae)
```

## Arguments

- mae:

  MultiAssayExperiment

## Value

named list with identifiers for each SE

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
get_MAE_identifiers(mae)
#> $`single-agent`
#> $`single-agent`$duration
#> [1] "Duration"
#> 
#> $`single-agent`$cellline
#> [1] "clid"
#> 
#> $`single-agent`$cellline_name
#> [1] "CellLineName"
#> 
#> $`single-agent`$cellline_tissue
#> [1] "Tissue"
#> 
#> $`single-agent`$cellline_ref_div_time
#> [1] "ReferenceDivisionTime"
#> 
#> $`single-agent`$cellline_parental_identifier
#> [1] "parental_identifier"
#> 
#> $`single-agent`$cellline_subtype
#> [1] "subtype"
#> 
#> $`single-agent`$drug
#> [1] "Gnumber"
#> 
#> $`single-agent`$drug_name
#> [1] "DrugName"
#> 
#> $`single-agent`$drug_moa
#> [1] "drug_moa"
#> 
#> $`single-agent`$untreated_tag
#> [1] "vehicle"   "untreated"
#> 
#> $`single-agent`$masked_tag
#> [1] "masked"
#> 
#> $`single-agent`$well_position
#> [1] "WellRow"    "WellColumn"
#> 
#> $`single-agent`$concentration
#> [1] "Concentration"
#> 
#> $`single-agent`$template
#> [1] "Template"  "Treatment"
#> 
#> $`single-agent`$barcode
#> [1] "Barcode" "Plate"  
#> 
#> $`single-agent`$drug2
#> [1] "Gnumber_2"
#> 
#> $`single-agent`$drug_name2
#> [1] "DrugName_2"
#> 
#> $`single-agent`$drug_moa2
#> [1] "drug_moa_2"
#> 
#> $`single-agent`$concentration2
#> [1] "Concentration_2"
#> 
#> $`single-agent`$drug3
#> [1] "Gnumber_3"
#> 
#> $`single-agent`$drug_name3
#> [1] "DrugName_3"
#> 
#> $`single-agent`$drug_moa3
#> [1] "drug_moa_3"
#> 
#> $`single-agent`$concentration3
#> [1] "Concentration_3"
#> 
#> $`single-agent`$data_source
#> [1] "data_source"
#> 
#> $`single-agent`$replicate
#> [1] "Replicate"
#> 
#> $`single-agent`$normalization_type
#> [1] "normalization_type"
#> 
#> 
```
