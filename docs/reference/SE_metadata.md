# Get and set metadata for parameters on a SummarizedExperiment object.

Set metadata for the fitting parameters that define the Metrics assay in
SummarizedExperiment object metadata.

## Usage

``` r
set_SE_fit_parameters(se, value)

set_SE_processing_metadata(se, value)

set_SE_keys(se, value)

set_SE_experiment_metadata(se, value, append = TRUE)

set_SE_experiment_raw_data(se, value)

get_SE_fit_parameters(se)

get_SE_processing_metadata(se)

get_SE_experiment_raw_data(se)

get_SE_experiment_metadata(se)

get_SE_keys(se, key_type = NULL)

get_SE_identifiers(se, id_type = NULL, simplify = TRUE)

set_SE_identifiers(se, value)
```

## Arguments

- se:

  a SummarizedExperiment object for which to add fit parameter metadata.

- value:

  named list of metadata for fit parameters.

- append:

  Boolean indicating whether to append the new metadata value to the
  existing entry.

- key_type:

  string of a specific key type (i.e. 'nested_keys', etc.).

- id_type:

  string of a specific id type (i.e. 'duration', 'cellline_name', etc.).

- simplify:

  Boolean indicating whether output should be simplified.

## Value

`se` with added metadata.

## Details

For `*et_SE_processing_metadata`, get/set metadata for the processing
info that defines the date_processed and packages versions in
SummarizedExperiment object metadata. For `*et_SE_fit_parameters`,
get/set metadata for fit parameters used to construct the Metrics assay
in a SummarizedExperiment object.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
get_SE_fit_parameters(se)
#> $n_point_cutoff
#> [1] 4
#> 
#> $range_conc
#> [1] 0.005 5.000
#> 
#> $force_fit
#> [1] FALSE
#> 
#> $pcutoff
#> [1] 0.05
#> 
#> $cap
#> [1] 0.1
#> 

mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
meta <- get_SE_processing_metadata(se)

mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
get_SE_experiment_raw_data(se)
#> NULL

mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
get_SE_experiment_metadata(se)
#> DataFrame with 0 rows and 0 columns

mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
get_SE_identifiers(se)
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
