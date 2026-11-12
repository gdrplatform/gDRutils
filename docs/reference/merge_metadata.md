# Merge metadata

Merge metadata

## Usage

``` r
merge_metadata(SElist, metadata_fields)
```

## Arguments

- SElist:

  named list of `SummarizedExperiment`s

- metadata_fields:

  vector of metadata names that will be merged

## Value

list of merged metadata

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
listSE <- list(
  se,
  se
)
metadata_fields <- identify_unique_se_metadata_fields(listSE)
merge_metadata(listSE, metadata_fields)
#> $identifiers
#> NULL
#> 
#> $experiment_metadata
#> $experiment_metadata$sources
#> list()
#> 
#> 
#> $Keys
#> NULL
#> 
#> $provenance
#> $provenance$schema_version
#> [1] 1
#> 
#> $provenance$fits
#> list()
#> 
#> $provenance$steps
#> list()
#> 
#> $provenance$sources
#> $provenance$sources[[1]]
#> $provenance$sources[[1]]$schema_version
#> [1] 1
#> 
#> $provenance$sources[[1]]$fits
#> $provenance$sources[[1]]$fits$Metrics
#> $provenance$sources[[1]]$fits$Metrics$gDR
#> $provenance$sources[[1]]$fits$Metrics$gDR$called_by
#> [1] "apply_fit"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$name
#> [1] "fit_drug_response_metrics"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$package
#> [1] "gDRutils"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$version
#> [1] "1.11.22"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$hash
#> [1] "ff16d78c11f4cdea"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$dependencies
#>          BumpyMatrix            checkmate            codetools 
#>             "1.18.0"              "2.3.4"             "0.2.20" 
#>           data.table               digest                  drc 
#>             "1.18.4"             "0.6.39"              "3.0.1" 
#>             gDRutils             jsonlite         jsonvalidate 
#>            "1.11.22"              "2.0.0"              "1.5.0" 
#> MultiAssayExperiment                  qs2            S4Vectors 
#>             "1.36.2"              "0.2.2"             "0.48.1" 
#>              stringr SummarizedExperiment 
#>              "1.6.0"             "1.40.0" 
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$info
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$info$model
#> [1] "3-parameter log-logistic (drc::LL.3u), upper asymptote fixed at 1"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$info$description
#> [1] "Hill curve of the response against concentration, fitted per normalization type. Falls back to a constant fit when fewer than n_point_cutoff concentrations are available or, unless force_fit = TRUE, when the F-test p-value is at least pcutoff."
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fn$info$reference
#> [1] "Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021"
#> 
#> 
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fit_args
#> $provenance$sources[[1]]$fits$Metrics$gDR$fit_args$capping_fold
#> [1] 5
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fit_args$range_conc
#> [1] 0.005 5.000
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fit_args$pcutoff
#> [1] 0.05
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fit_args$n_point_cutoff
#> [1] 4
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fit_args$force_fit
#> [1] FALSE
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$fit_args$cap
#> [1] 0.1
#> 
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$input_assay
#> [1] "Averaged"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$output_assay
#> [1] "Metrics"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$slicing_cols
#> [1] "normalization_type"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$slicing_values
#> [1] "GR" "RV"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$merge
#> [1] "replace"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$run
#> $provenance$sources[[1]]$fits$Metrics$gDR$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$run$gDRutils
#> [1] "1.11.22"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$run$input_hash
#> [1] "077d8eabaafbd390"
#> 
#> $provenance$sources[[1]]$fits$Metrics$gDR$run$n_rows
#> [1] 200
#> 
#> 
#> 
#> 
#> 
#> $provenance$sources[[1]]$steps
#> $provenance$sources[[1]]$steps$create
#> $provenance$sources[[1]]$steps$create$settings
#> $provenance$sources[[1]]$steps$create$settings$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[1]]$steps$create$settings$readout
#> [1] "ReadoutValue"
#> 
#> $provenance$sources[[1]]$steps$create$settings$nested_identifiers
#> [1] "Concentration"
#> 
#> $provenance$sources[[1]]$steps$create$settings$nested_confounders
#> [1] "Barcode"
#> 
#> $provenance$sources[[1]]$steps$create$settings$override_untrt_controls
#> NULL
#> 
#> $provenance$sources[[1]]$steps$create$settings$masked
#> $provenance$sources[[1]]$steps$create$settings$masked$n_wells
#> [1] 0
#> 
#> 
#> 
#> $provenance$sources[[1]]$steps$create$run
#> $provenance$sources[[1]]$steps$create$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[1]]$steps$create$run$gDRutils
#> [1] "1.11.22"
#> 
#> 
#> 
#> $provenance$sources[[1]]$steps$normalize
#> $provenance$sources[[1]]$steps$normalize$functions
#> $provenance$sources[[1]]$steps$normalize$functions$control_mean_fxn
#> $provenance$sources[[1]]$steps$normalize$functions$control_mean_fxn$name
#> [1] "control_mean_fxn"
#> 
#> $provenance$sources[[1]]$steps$normalize$functions$control_mean_fxn$package
#> [1] "gDRcore"
#> 
#> $provenance$sources[[1]]$steps$normalize$functions$control_mean_fxn$version
#> [1] "1.11.26"
#> 
#> $provenance$sources[[1]]$steps$normalize$functions$control_mean_fxn$hash
#> [1] "33fb9ece036126b7"
#> 
#> $provenance$sources[[1]]$steps$normalize$functions$control_mean_fxn$source
#> [1] "function (x) \n{\n    mean(x, trim = 0.25)\n}"
#> 
#> 
#> 
#> $provenance$sources[[1]]$steps$normalize$settings
#> $provenance$sources[[1]]$steps$normalize$settings$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$nested_identifiers
#> [1] "Concentration"
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$nested_confounders
#> [1] "Barcode"
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$ndigit_rounding
#> [1] 4
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$control_assay
#> [1] "Controls"
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$raw_treated_assay
#> [1] "RawTreated"
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$normalized_assay
#> [1] "Normalized"
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$gr
#> $provenance$sources[[1]]$steps$normalize$settings$gr$n_pairs
#>         day0 ref_div_time         none 
#>            0          100            0 
#> 
#> $provenance$sources[[1]]$steps$normalize$settings$gr$ref_div_time_cell_lines
#>  [1] "cellline_BA" "cellline_CA" "cellline_DA" "cellline_EA" "cellline_FA"
#>  [6] "cellline_GB" "cellline_HB" "cellline_IB" "cellline_JB" "cellline_KB"
#> 
#> 
#> 
#> $provenance$sources[[1]]$steps$normalize$run
#> $provenance$sources[[1]]$steps$normalize$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[1]]$steps$normalize$run$gDRutils
#> [1] "1.11.22"
#> 
#> 
#> 
#> $provenance$sources[[1]]$steps$average
#> $provenance$sources[[1]]$steps$average$functions
#> $provenance$sources[[1]]$steps$average$functions$average_FUN
#> $provenance$sources[[1]]$steps$average$functions$average_FUN$name
#> [1] "average_FUN"
#> 
#> $provenance$sources[[1]]$steps$average$functions$average_FUN$package
#> [1] "gDRcore"
#> 
#> $provenance$sources[[1]]$steps$average$functions$average_FUN$version
#> [1] "1.11.26"
#> 
#> $provenance$sources[[1]]$steps$average$functions$average_FUN$hash
#> [1] "6f45617f8192fe7c"
#> 
#> $provenance$sources[[1]]$steps$average$functions$average_FUN$dependencies
#>         BiocParallel          BumpyMatrix            checkmate 
#>             "1.44.0"             "1.18.0"              "2.3.4" 
#>           data.table                  drc        futile.logger 
#>             "1.18.4"              "3.0.1"              "1.4.9" 
#>              gDRcore             gDRutils             jsonlite 
#>            "1.11.26"            "1.11.22"              "2.0.0" 
#> MultiAssayExperiment                purrr            S4Vectors 
#>             "1.36.2"              "1.2.2"             "0.48.1" 
#>              stringr SummarizedExperiment 
#>              "1.6.0"             "1.40.0" 
#> 
#> 
#> 
#> $provenance$sources[[1]]$steps$average$settings
#> $provenance$sources[[1]]$steps$average$settings$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[1]]$steps$average$settings$series_identifiers
#> [1] "Concentration"
#> 
#> $provenance$sources[[1]]$steps$average$settings$trt_keys
#> [1] "CellLineName" "clid"         "drug_moa"     "DrugName"     "Duration"    
#> [6] "Gnumber"      "Tissue"      
#> 
#> $provenance$sources[[1]]$steps$average$settings$normalized_assay
#> [1] "Normalized"
#> 
#> $provenance$sources[[1]]$steps$average$settings$averaged_assay
#> [1] "Averaged"
#> 
#> 
#> $provenance$sources[[1]]$steps$average$run
#> $provenance$sources[[1]]$steps$average$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[1]]$steps$average$run$gDRutils
#> [1] "1.11.22"
#> 
#> 
#> 
#> 
#> 
#> $provenance$sources[[2]]
#> $provenance$sources[[2]]$schema_version
#> [1] 1
#> 
#> $provenance$sources[[2]]$fits
#> $provenance$sources[[2]]$fits$Metrics
#> $provenance$sources[[2]]$fits$Metrics$gDR
#> $provenance$sources[[2]]$fits$Metrics$gDR$called_by
#> [1] "apply_fit"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$name
#> [1] "fit_drug_response_metrics"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$package
#> [1] "gDRutils"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$version
#> [1] "1.11.22"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$hash
#> [1] "ff16d78c11f4cdea"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$dependencies
#>          BumpyMatrix            checkmate            codetools 
#>             "1.18.0"              "2.3.4"             "0.2.20" 
#>           data.table               digest                  drc 
#>             "1.18.4"             "0.6.39"              "3.0.1" 
#>             gDRutils             jsonlite         jsonvalidate 
#>            "1.11.22"              "2.0.0"              "1.5.0" 
#> MultiAssayExperiment                  qs2            S4Vectors 
#>             "1.36.2"              "0.2.2"             "0.48.1" 
#>              stringr SummarizedExperiment 
#>              "1.6.0"             "1.40.0" 
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$info
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$info$model
#> [1] "3-parameter log-logistic (drc::LL.3u), upper asymptote fixed at 1"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$info$description
#> [1] "Hill curve of the response against concentration, fitted per normalization type. Falls back to a constant fit when fewer than n_point_cutoff concentrations are available or, unless force_fit = TRUE, when the F-test p-value is at least pcutoff."
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fn$info$reference
#> [1] "Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021"
#> 
#> 
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fit_args
#> $provenance$sources[[2]]$fits$Metrics$gDR$fit_args$capping_fold
#> [1] 5
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fit_args$range_conc
#> [1] 0.005 5.000
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fit_args$pcutoff
#> [1] 0.05
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fit_args$n_point_cutoff
#> [1] 4
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fit_args$force_fit
#> [1] FALSE
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$fit_args$cap
#> [1] 0.1
#> 
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$input_assay
#> [1] "Averaged"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$output_assay
#> [1] "Metrics"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$slicing_cols
#> [1] "normalization_type"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$slicing_values
#> [1] "GR" "RV"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$merge
#> [1] "replace"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$run
#> $provenance$sources[[2]]$fits$Metrics$gDR$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$run$gDRutils
#> [1] "1.11.22"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$run$input_hash
#> [1] "077d8eabaafbd390"
#> 
#> $provenance$sources[[2]]$fits$Metrics$gDR$run$n_rows
#> [1] 200
#> 
#> 
#> 
#> 
#> 
#> $provenance$sources[[2]]$steps
#> $provenance$sources[[2]]$steps$create
#> $provenance$sources[[2]]$steps$create$settings
#> $provenance$sources[[2]]$steps$create$settings$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[2]]$steps$create$settings$readout
#> [1] "ReadoutValue"
#> 
#> $provenance$sources[[2]]$steps$create$settings$nested_identifiers
#> [1] "Concentration"
#> 
#> $provenance$sources[[2]]$steps$create$settings$nested_confounders
#> [1] "Barcode"
#> 
#> $provenance$sources[[2]]$steps$create$settings$override_untrt_controls
#> NULL
#> 
#> $provenance$sources[[2]]$steps$create$settings$masked
#> $provenance$sources[[2]]$steps$create$settings$masked$n_wells
#> [1] 0
#> 
#> 
#> 
#> $provenance$sources[[2]]$steps$create$run
#> $provenance$sources[[2]]$steps$create$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[2]]$steps$create$run$gDRutils
#> [1] "1.11.22"
#> 
#> 
#> 
#> $provenance$sources[[2]]$steps$normalize
#> $provenance$sources[[2]]$steps$normalize$functions
#> $provenance$sources[[2]]$steps$normalize$functions$control_mean_fxn
#> $provenance$sources[[2]]$steps$normalize$functions$control_mean_fxn$name
#> [1] "control_mean_fxn"
#> 
#> $provenance$sources[[2]]$steps$normalize$functions$control_mean_fxn$package
#> [1] "gDRcore"
#> 
#> $provenance$sources[[2]]$steps$normalize$functions$control_mean_fxn$version
#> [1] "1.11.26"
#> 
#> $provenance$sources[[2]]$steps$normalize$functions$control_mean_fxn$hash
#> [1] "33fb9ece036126b7"
#> 
#> $provenance$sources[[2]]$steps$normalize$functions$control_mean_fxn$source
#> [1] "function (x) \n{\n    mean(x, trim = 0.25)\n}"
#> 
#> 
#> 
#> $provenance$sources[[2]]$steps$normalize$settings
#> $provenance$sources[[2]]$steps$normalize$settings$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$nested_identifiers
#> [1] "Concentration"
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$nested_confounders
#> [1] "Barcode"
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$ndigit_rounding
#> [1] 4
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$control_assay
#> [1] "Controls"
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$raw_treated_assay
#> [1] "RawTreated"
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$normalized_assay
#> [1] "Normalized"
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$gr
#> $provenance$sources[[2]]$steps$normalize$settings$gr$n_pairs
#>         day0 ref_div_time         none 
#>            0          100            0 
#> 
#> $provenance$sources[[2]]$steps$normalize$settings$gr$ref_div_time_cell_lines
#>  [1] "cellline_BA" "cellline_CA" "cellline_DA" "cellline_EA" "cellline_FA"
#>  [6] "cellline_GB" "cellline_HB" "cellline_IB" "cellline_JB" "cellline_KB"
#> 
#> 
#> 
#> $provenance$sources[[2]]$steps$normalize$run
#> $provenance$sources[[2]]$steps$normalize$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[2]]$steps$normalize$run$gDRutils
#> [1] "1.11.22"
#> 
#> 
#> 
#> $provenance$sources[[2]]$steps$average
#> $provenance$sources[[2]]$steps$average$functions
#> $provenance$sources[[2]]$steps$average$functions$average_FUN
#> $provenance$sources[[2]]$steps$average$functions$average_FUN$name
#> [1] "average_FUN"
#> 
#> $provenance$sources[[2]]$steps$average$functions$average_FUN$package
#> [1] "gDRcore"
#> 
#> $provenance$sources[[2]]$steps$average$functions$average_FUN$version
#> [1] "1.11.26"
#> 
#> $provenance$sources[[2]]$steps$average$functions$average_FUN$hash
#> [1] "6f45617f8192fe7c"
#> 
#> $provenance$sources[[2]]$steps$average$functions$average_FUN$dependencies
#>         BiocParallel          BumpyMatrix            checkmate 
#>             "1.44.0"             "1.18.0"              "2.3.4" 
#>           data.table                  drc        futile.logger 
#>             "1.18.4"              "3.0.1"              "1.4.9" 
#>              gDRcore             gDRutils             jsonlite 
#>            "1.11.26"            "1.11.22"              "2.0.0" 
#> MultiAssayExperiment                purrr            S4Vectors 
#>             "1.36.2"              "1.2.2"             "0.48.1" 
#>              stringr SummarizedExperiment 
#>              "1.6.0"             "1.40.0" 
#> 
#> 
#> 
#> $provenance$sources[[2]]$steps$average$settings
#> $provenance$sources[[2]]$steps$average$settings$data_type
#> [1] "single-agent"
#> 
#> $provenance$sources[[2]]$steps$average$settings$series_identifiers
#> [1] "Concentration"
#> 
#> $provenance$sources[[2]]$steps$average$settings$trt_keys
#> [1] "CellLineName" "clid"         "drug_moa"     "DrugName"     "Duration"    
#> [6] "Gnumber"      "Tissue"      
#> 
#> $provenance$sources[[2]]$steps$average$settings$normalized_assay
#> [1] "Normalized"
#> 
#> $provenance$sources[[2]]$steps$average$settings$averaged_assay
#> [1] "Averaged"
#> 
#> 
#> $provenance$sources[[2]]$steps$average$run
#> $provenance$sources[[2]]$steps$average$run$R
#> [1] "4.5.1"
#> 
#> $provenance$sources[[2]]$steps$average$run$gDRutils
#> [1] "1.11.22"
#> 
#> 
#> 
#> 
#> 
#> 
#> 
#> $fit_parameters
#> NULL
#> 
#> $.internal
#> $.internal$date_processed
#> [1] "2026-10-26"
#> 
#> $.internal$versions
#> $.internal$versions$R
#> [1] "4.5.1"
#> 
#> $.internal$versions$gDRcore
#> [1] "1.11.26"
#> 
#> $.internal$versions$gDRutils
#> [1] "1.11.22"
#> 
#> 
#> 
```
