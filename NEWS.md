## gDRutils 1.3.6 - 2024-07-17
* update `define_matrix_grid_positions`

## gDRutils 1.3.5 - 2024-07-12
* move `get_combo_col_settings` and `get_iso_colors` to `gDRplots` package

## gDRutils 1.3.4 - 2024-07-08
* add residual sum of square and p-value to Metrics assay

## gDRutils 1.3.3 - 2024-07-03
* add vignette section about prettifying logic

## gDRutils 1.3.2 - 2024-06-24
* fixed issue in vignette 

## gDRutils 1.3.1 - 2024-05-27
* synchronize Bioconductor and GitHub versioning

## gDRutils 1.1.17 - 2024-05-22
* fix check notes

## gDRutils 1.1.16 - 2024-05-22
* move `convert_se_assay_to_custom_dt`, `capVals` and `get_settings_from_json` from `gDRcomponents` package

## gDRutils 1.1.15 - 2024-05-17
* simplify logic of supported experiments

## gDRutils 1.1.14 - 2024-05-16
* move `define_matrix_grid_positions` and `round_concentration` functions from `gDRcore` package

## gDRutils 1.1.13 - 2024-05-08
* fix typo

## gDRutils 1.1.12 - 2024-04-30
* move validator functions from `gDRcomponents` to `gDRutils`

## gDRutils 1.1.11 - 2024-04-15
* add `get_testdata_combo` and `get_testdata_codilution`

## gDRutils 1.1.10 - 2024-03-07
* clean up the package

## gDRutils 1.1.9 - 2024-03-07
* simplify keywords

## gDRutils 1.1.8 - 2024-02-28
* add fit_source to header list

## gDRutils 1.1.7 - 2024-02-26
* improve pkgdown site
  * improved references
  * valid NEWS.md

## gDRutils 1.1.6 - 2024-02-22
* restore tooltips in table

## gDRutils 1.1.5 - 2024-02-14
* make documentation compatible with pkdgdown

## gDRutils 1.1.4 - 2024-01-30
* rename `matrix` into `combination`

## gDRutils 1.1.3 - 2024-01-22
* add new description fields

## gDRutils 1.1.2 - 2023-12-01
* fix bug with refining rowData
* extend the list of headers

## gDRutils 1.1.1 - 2023-11-22
* sync master with devel branch
* update schema to support NA in reference division time
* add minor fix in code styling
* add new function `gemoetric_mean`
* transform values into numeric in `predict_efficacy_from_conc` function 
* add "Treatment" as template identifier

## gDRutils 1.1.0 - 2023-10-24
* release Bioc 3.18

## gDRutils 1.0.0 - 2023-10-24
* prerelease Bioc 3.18

## gDRutils 0.99.34 - 2023-10-18
* adjust NEWS to Bioc format

## gDRutils 0.99.33 - 2023-10-09
* add support for flattening averaged assays

## gDRutils 0.99.32 - 2023-09-22
* fix bug in the case of conc=0 for evaluating efficacy

## gDRutils 0.99.31 - 2023-09-19
* add `wide_structure` param to `convert_mae_assay_to_dt`

## gDRutils 0.99.30 - 2023-09-08
* updated `experimentalist` description in schema

## gDRutils 0.99.29 - 2023-09-05
* add `Replicate` as a new identifier

## gDRutils 0.99.28 - 2023-09-05
* improve the logic of standardize_MAE to keep SE-specific metadata and be able to revert standardization

## gDRutils 0.99.27 - 2023-08-01
* keep unchanged names in DataFrame

## gDRutils 0.99.26 - 2023-08-01
* tidy code

## gDRutils 0.99.25 - 2023-06-27
* add assert for missing rownames

## gDRutils 0.99.24 - 2023-06-22
* replaced RDS with qs

## gDRutils 0.99.23 - 2023-06-20
* fix check in R 4.3

## gDRutils 0.99.22 - 2023-06-12
* switch from `merge` to `[[`

## gDRutils 0.99.21 - 2023-06-12
* replace `order` with `data.table::setorder`
* add support for custom identifiers in merge_SE

## gDRutils 0.99.20 - 2023-06-07
* switch from `aggregate` to `data.table`

## gDRutils 0.99.19 - 2023-06-06
* replaced reshape2 functions by functions from data.table

## gDRutils 0.99.18 - 2023-05-22
* format the vignette with BiocStyle

## gDRutils 0.99.17 - 2023-05-22
* fix related with data.table
* remove `.get_treated_conditions` and `.get_untreated_conditions`

## gDRutils 0.99.16 - 2023-05-18
* add support for merging combination-data assays

## gDRutils 0.99.15 - 2023-05-12
* update after unifying normalization types

## gDRutils 0.99.14 - 2023-05-12
* fix lintr

## gDRutils 0.99.13 - 2023-05-09
* removed `cotreatment` entry from EXPERIMENT_GROUPS

## gDRutils 0.99.12 - 2023-05-09
* fix bug in `convert_mae_assay_to_dt`

## gDRutils 0.99.11 - 2023-05-08
* refactor code with single ampersand in if statements

## gDRutils 0.99.10 - 2023-04-28
* change order of untreated tags

## gDRutils 0.99.9 - 2023-04-24
* changed data.frame to data.table

## gDRutils 0.99.8 - 2023-04-20
* switch to OSI license

## gDRutils 0.99.7 - 2023-04-20
* clean-up vignette

## gDRutils 0.99.6 - 2023-04-18
* extend the logic of `apply_bumpy_function`

## gDRutils 0.99.5 - 2023-04-17
* add R 4.2 as a dependency

## gDRutils 0.99.4 - 2023-04-14
* fix examples

## gDRutils 0.99.3 - 2023-04-13
* make linter happy

## gDRutils 0.99.2 - 2023-04-12
* add licence

## gDRutils 0.99.1 - 2023-04-07
* update maintainer

## gDRutils 0.99.0 - 2023-03-28
* downgrade version to make it Bioconductor compatible

## gDRutils 0.1.3.22 - 2023-03-13
* tidy code

## gDRutils 0.1.3.21 - 2023-03-09
* better handling of corner cases for the single-agent fitting

## gDRutils 0.1.3.20 - 2023-03-08
* add support for splitting normalization data types

## gDRutils 0.1.3.19 - 2023-03-08
* restore necessary functions

## gDRutils 0.1.3.18 - 2023-03-01
* add examples of identifiers

## gDRutils 0.1.3.17 - 2023-03-01
* remove obsolete code

## gDRutils 0.1.3.16 - 2023-02-22
* tidy code

## gDRutils 0.1.3.15 - 2023-02-10
* add tests for `convert_se_ref_assay_to_dt`

## gDRutils 0.1.3.14 - 2023-01-10
* add helpers for colData/rowData refinement

## gDRutils 0.1.3.13 - 2022-12-22
* fix warnings in unit tests

## gDRutils 0.1.3.12 - 2022-12-20
* `R CMD check` returns `0 errors ✓ | 0 warnings ✓ | 0 notes ✓` 

## gDRutils 0.1.3.11 - 2022-12-15
* default apply_bumpy_function parallelize to TRUE 

## gDRutils 0.1.3.10 - 2022-12-08
* add apply_bumpy_function support

## gDRutils 0.1.3.9 - 2022-12-07
* add co-dilution to single-agent group

## gDRutils 0.1.3.8 - 2022-11-30
* bugfix in validate_SE

## gDRutils 0.1.3.7 - 2022-10-18
* add update_env_idfs_from_mae function

## gDRutils 0.1.3.6 - 2022-10-04
* add helper for dealing with idfs' synonyms

## gDRutils 0.1.3.5 - 2022-09-21
* promote or demote fields in a BumpyMatrix object and perform summarization on assays.

## gDRutils 0.1.3.4 - 2022-08-11
* major improvements in JSON validation/conversion logic 
  * provide info about JSON schemes via env variables
  * convert and validate MAE summary next to the SE experiments

## gDRutils 0.1.3.3 - 2022-07-25
* move json validation/convertion logic from gDRelastic

## gDRutils 0.1.3.2 - 2022-07-11
* standardize MAE using default gDR identifiers

## gDRutils 0.1.3.1 - 2022-06-29
* remove adding integer identifiers at the end of colnames/rownames

## gDRutils 0.1.3.0 - 2022-06-02
* release 1.3.0

## gDRutils 0.1.0.48 - 2022-05-27
* correct recognition of empty SE

## gDRutils 0.1.0.47 - 2022-05-25
* remove redundant validation of rownames in SE

## gDRutils 0.1.0.46 - 2022-05-03
* update prettify function
* fix hardcoded identifiers in validated_SE

## gDRutils 0.1.0.45 - 2022-04-29
* avoid using `grep` for getting cotreatment identifiers

## gDRutils 0.1.0.44 - 2022-04-26
* switched from unnamed to named vector of experiment groups for `single-agent`
* fix the logic in validating single-agent experiments

## gDRutils 0.1.0.43 - 2022-04-13
* set r2 value to NA for invalid and 0 for constant fits

## gDRutils 0.1.0.42 - 2022-04-08
* added identifier descriptions

## gDRutils 0.1.0.41 - 2022-04-08
* fix wrong order of elements in rownames in SE

## gDRutils 0.1.0.40 - 2022-04-07
* add cap_ic50 function

## gDRutils 0.1.0.39 - 2022-03-31
* extend possible `Barcode` identifiers

## gDRutils 0.1.0.38 - 2022-03-30
* fix hardcoded identifiers in validate SE

## gDRutils 0.1.0.37 - 2022-03-28
* add a space between two-word cotreatment identifiers

## gDRutils 0.1.0.36 - 2022-03-24
* remove all R CMD check warnings

## gDRutils 0.1.0.35 - 2022-03-22
* change prettify functions to not substitute metadata

## gDRutils 0.1.0.34 - 2022-03-21
* add helper function for MAE/experiments

## gDRutils 0.1.0.33 - 2022-03-19
* move constant fit warning

## gDRutils 0.1.0.32 - 2022-03-18
* add .calculate_complement 

## gDRutils 0.1.0.31 - 2022-03-16
* move duplication warning 

## gDRutils 0.1.0.30 - 2022-03-14
* add support for identifier validation

## gDRutils 0.1.0.29 - 2022-03-03
* add getter and setter for experiment_raw_data

## gDRutils 0.1.0.28 - 2022-02-16
* refactor identifier `drugname` to `drug_name`

## gDRutils 0.1.0.27 - 2022-02-10
* add `MAEpply` function

## gDRutils 0.1.0.26 - 2022-02-01
* update documentation

## gDRutils 0.1.0.25 - 2022-01-31
* support SE-  init for different BioC versions

## gDRutils 0.1.0.24 - 2022-01-25
* standardize/improve CI

## gDRutils 0.1.0.23 - 2022-01-25
* switch unit tests from SE to MAE from `gDRtestData`

## gDRutils 0.1.0.22 - 2022-01-07
* update assert in `validate_MAE`

## gDRutils 0.1.0.21 - 2022-01-07
* update SE-related functions of gDRutils to support MAE

## gDRutils 0.1.0.20 - 2021-12-21
* move combo-related functions from gDRviz
* add unit tests for combo-related functions

## gDRutils 0.1.0.19 - 2021-11-04
* do not create a nested list of identifiers during mergin SE

## gDRutils 0.1.0.18 - 2021-11-01
* fix issued with new SummarizedExperiment

## gDRutils 0.1.0.17 - 2021-10-22
* refactor: generalize prediction help functions for fits

## gDRutils 0.1.0.16 - 2021-09-27
* update validate_SE as per combos

## gDRutils 0.1.0.15 - 2021-09-22
* feat: use 'GRvalue' and 'RelativeViability' as normalization_types in 'fit_curves'

## gDRutils 0.1.0.14 - 2021-09-15
* fix bug with getting_SE_identifiers for untreated_tag

## gDRutils 0.1.0.13 - 2021-09-13
* fix obsolete arguments in `reset_env_identifiers`

## gDRutils 0.1.0.12 - 2021-09-07
* fix bug with getting identifiers for untreated_tag

## gDRutils 0.1.0.11 - 2021-08-26
* move NA logic for elements of BumpyMatrix with no data to fit_curve

## gDRutils 0.1.0.10 - 2021-08-25
* remove additional ordering line in df_to_bm_assay.R

## gDRutils 0.1.0.9 - 2021-08-19
* update list of available identifiers - 2nd and 3rd drug, data_source
* update the logic for get_identifier and .get_id
* add prettified_identifier- 

## gDRutils 0.1.0.8 - 2021-08-13
* fix bug with wrong order of rows and cols in `df_to_bm_assay`

## gDRutils 0.1.0.7 - 2021-07-06
* remove hard version equality for pkg deps 

## gDRutils 0.1.0.6 - 2021-06-29
* add barcode as identifiers and store identifiers within split_SE_components

## gDRutils 0.1.0.5 - 2021-06-25
* update the logic for CI/CD - repos fetching is now handled with technical user

## gDRutils 0.1.0.4 - 2021-06-18
* add concentration and template as additional identifiers

## gDRutils 0.1.0.3 - 2021-06-21
* refactor logisticFit with error handling

## gDRutils 0.1.0.2 - 2021-06-18
* remove deprecated functions
* switch from getMetadata to split_SE_components
* refactor df_to_bm_assay

## gDRutils 0.1.0.1 - 2021-06-02
* upgrade validate_SE by checking if rowData and colData do not have empty strings
* isolate flatten function

## gDRutils 0.1.0.0 - 2021-06-02
* release 1.0.0

## gDRutils 0.0.0.49 - 2021-06-02
* fix bug with merging flatten assays

## gDRutils 0.0.0.48 - 2021-05-24
* refactor validate_SE function

## gDRutils 0.0.0.47 - 2021-05-18
* add `merge_SE` function

## gDRutils 0.0.0.46 - 2021-05-18
* fix typo in SE validator function

## gDRutils 0.0.0.45 - 2021-05-18
* add SE validator function

## gDRutils 0.0.0.44 - 2021-05-13
* refactor `prettify_flat_metrics` function

## gDRutils 0.0.0.43 - 2021-04-30
* remove `Metrics_rownames` during flattening data.frame/data.table

## gDRutils 0.0.0.42 - 2021-04-29
* add `prettify_flat_metrics` function

## gDRutils 0.0.0.41 - 2021-04-29
* bugfix flattening data.tables

## gDRutils 0.0.0.40 - 2021-04-29
* bugfix data.table merge in `convert_se_assay_ref_to_dt`

## gDRutils 0.0.0.39 - 2021-04-27
* add support for flattening data.tables

## gDRutils 0.0.0.38 - 2021-04-21
* switch from `processing_metadata` to `.internal`

## gDRutils 0.0.0.37 - 2021-04-20
* add support for getting and setting processing info metadata

## gDRutils 0.0.0.36 - 2021-04-19
* sort BumpyMatrix created from dt

## gDRutils 0.0.0.35 - 2021-04-14
* revert metric arguments in 'convert_se_assay_to_dt'

## gDRutils 0.0.0.34 - 2021-04-09
* add support for getting and setting fit parameter metadata

## gDRutils 0.0.0.33 - 2021-04-07
* standardize metric header names

## gDRutils 0.0.0.32 - 2021-04-07
* add more options for returned data - with 'Metrics' assay in the case of `convert_se_assay_to_dt`
* add 'convert_se_ref_assay_to_dt' function

## gDRutils 0.0.0.31 - 2021-03-30
* add support for getting and setting metadata on the SummarizedExperiment object

## gDRutils 0.0.0.30 - 2021-03-30
* add support for getting and setting metadata on the SummarizedExperiment object
* add numeric type assertions and tests for logistic_4parameters

## gDRutils 0.0.0.29 - 2021-03-23
* remove positional naming dependence on assay_to_dt for - merge_metrics = TRUE argument

## gDRutils 0.0.0.28 - 2021-03-16
* minor refactor to assay_to_dt
* deprecated support for assay_to_dt- include_controls = TRUE argument - not backwards compatible

## gDRutils 0.0.0.27 - 2021-03-09
* update .estimate_xc50 and add tests

## gDRutils 0.0.0.26 - 2021-03-04
* fix bug with wrong class for S3 methods in convert_assay_data_to_dt
* update the documentation

## gDRutils 0.0.0.25 - 2021-03-02
* fix bug with missing class in assert for matrix

## gDRutils 0.0.0.24 - 2021-02-10
* modify fit_curves to take in flexible curve_type- s
* clean-up assay_to_dt

## gDRutils 0.0.0.23 - 2021-02-04
* added linter

## gDRutils 0.0.0.22 - 2021-02-03
* move assay_to_dt from gDR to gDRutils
* refactor assay_to_dt to support two assay types
  * list of DFrame- s
  * BumpyMatrix objects
* move .get_treated_conditions and .get_untreated_conditions from gDR to gDRutils 

## gDRutils 0.0.0.21 - 2021-01-26
* allow for DFrame as input

## gDRutils 0.0.0.20 - 2021-01-26
* 's/BumpyMatrix::splitToBumpyMatrix/BumpyMatrix::splitAsBumpyMatrix/'

## gDRutils 0.0.0.19 - 2021-01-20
* refactor RVGRfits and rename to fit_curves
* make logisticFit function independent of curve_type

## gDRutils 0.0.0.18 - 2021-01-19
* move df_to_assay.R from gDRcore to gDRutils
* move df_to_bm_assay.R from gDRcore to gDRutils
* update identifiers_list.R

## gDRutils 0.0.0.15 - 2020-12-14
* refactor get_identifiers and get_headers to be settable and cached

## gDRutils 0.0.0.14 - 2020-12-14
* totally, finally, and unscrupulously remove dplyr package from gDRutils
* replace IC by RV

## gDRutils 0.0.0.13 - 2020-10-12
* add CI

## gDRutils 0.0.0.12 - 2020-10-09
* fix assay_to_df

## gDRutils 0.0.0.11 - 2020-10-09
* update namespaces

## gDRutils 0.0.0.10 - 2020-10-05
* minor refactor

## gDRutils 0.0.0.9 - 2020-09-14
* add small fixes for assays with empty DataFrame

## gDRutils 0.0.0.8 - 2020-09-02
* update logic as per new db model

## gDRutils 0.0.0.7 - 2020-08-07
* update variable names as per new db model

## gDRutils 0.0.0.6 - 2020-07-02
* import pipes from magrittr

## gDRutils 0.0.0.4 - 2020-06-10
* including the masked field to be able to remove the masked data from averages
