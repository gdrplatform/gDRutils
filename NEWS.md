## 1.1.6 (2024-02-22)
- add fit_source to header list

## 1.1.5 (2024-02-14)
- make documentation compatible with pkdgdown

## 1.1.4 (2024-01-30)
- rename `matrix` into `combination`

## 1.1.3 (2024-01-22)
- add new description fields

## 1.1.2 (2023-12-01)
- fix bug with refining rowData
- extend the list of headers

## 1.1.1 (2023-11-22)
- sync master with devel branch
- update schema to support NA in reference division time
- added minor fix in code styling
- added new function `gemoetric_mean`
- transform values into numeric in `predict_efficacy_from_conc` function 
- add "Treatment" as template identifier

## 1.1.0 (2023-10-24)
- release Bioc 3.18

## 1.0.0 (2023-10-24)
- prerelease Bioc 3.18

## 0.99.34 (2023-10-18)
- adjust NEWS to Bioc format

## 0.99.33 (2023-10-09)
- add support for flattening averaged assays

## 0.99.32 (2023-09-22)
- fix bug in the case of conc=0 for evaluating efficacy

## 0.99.31 (2023-09-19)
- add `wide_structure` param to `convert_mae_assay_to_dt`

## 0.99.30 (2023-09-08)
- updated `experimentalist` description in schema

## 0.99.29 (2023-09-05)
- add `Replicate` as a new identifier

## 0.99.28 (2023-09-05)
- improve the logic of standardize_MAE to keep SE-specific metadata and be able to revert standardization

## 0.99.27 (2023-08-01)
- keep unchanged names in DataFrame

## 0.99.26 (2023-08-01)
- tidy code

## 0.99.25 (2023-06-27)
- add assert for missing rownames

## 0.99.24 (2023-06-22)
- replaced RDS with qs

## 0.99.23 (2023-06-20)
- fix check in R 4.3

## 0.99.22 (2023-06-12)
- switch from `merge` to `[[`

## 0.99.21 (2023-06-12)
- replace `order` with `data.table::setorder`
- add support for custom identifiers in merge_SE

## 0.99.20 (2023-06-07)
- switch from `aggregate` to `data.table`

## 0.99.19 (2023-06-06)
- replaced reshape2 functions by functions from data.table

## 0.99.18 (2023-05-22)
- format the vignette with BiocStyle

## 0.99.17 (2023-05-22)
- fix related with data.table
- remove `.get_treated_conditions` and `.get_untreated_conditions`

## 0.99.16 (2023-05-18)
- add support for merging combination-data assays

## 0.99.15 (2023-05-12)
- update after unifying normalization types

## 0.99.14 (2023-05-12)
- fix lintr

## 0.99.13 (2023-05-09)
- removed `cotreatment` entry from EXPERIMENT_GROUPS

## 0.99.12 (2023-05-09)
- fix bug in `convert_mae_assay_to_dt`

## 0.99.11 (2023-05-08)
- refactor code with single ampersand in if statements

## 0.99.10 (2023-04-28)
- change order of untreated tags

## 0.99.9 (2023-04-24)
- changed data.frame to data.table

## 0.99.8 (2023-04-20)
- switch to OSI license

## 0.99.7 (2023-04-20)
- clean-up vignette

## 0.99.6 (2023-04-18)
- extend the logic of `apply_bumpy_function`

## 0.99.5 (2023-04-17)
- add R 4.2 as a dependency

## 0.99.4 (2023-04-14)
- fix examples

## 0.99.3 (2023-04-13)
- make linter happy

## 0.99.2 (2023-04-12)
- add licence

## 0.99.1 (2023-04-07)
- update maintainer

## 0.99.0 (2023-03-28)
- downgrade version to make it Bioconductor compatible

## 0.1.3.22 (2023-03-13)
- tidy code

## 0.1.3.21 (2023-03-09)
- better handling of corner cases for the single-agent fitting

## 0.1.3.20 (2023-03-08)
- add support for splitting normalization data types

## 0.1.3.19 (2023-03-08)
- restore necessary functions

## 0.1.3.18 (2023-03-01)
- add examples of identifiers

## 0.1.3.17 (2023-03-01)
- remove obsolete code

## 0.1.3.16 (2023-02-22)
- tidy code

## 0.1.3.15 (2023-02-10)
- add tests for `convert_se_ref_assay_to_dt`

## 0.1.3.14 (2023-01-10)
- add helpers for colData/rowData refinement

## 0.1.3.13 (2022-12-22)
- fix warnings in unit tests

## 0.1.3.12 (2022-12-20)
- `R CMD check` returns `0 errors ✓ | 0 warnings ✓ | 0 notes ✓` 

## 0.1.3.11 (2022-12-15)
- default apply_bumpy_function parallelize to TRUE 

## 0.1.3.10 (2022-12-08)
- add apply_bumpy_function support

## 0.1.3.9 (2022-12-07)
- add co-dilution to single-agent group

## 0.1.3.8 (2022-11-30)
- bugfix in validate_SE

## 0.1.3.7 (2022-10-18)
- add update_env_idfs_from_mae function

## 0.1.3.6 (2022-10-04)
- add helper for dealing with idfs' synonyms

## 0.1.3.5 (2022-09-21)
- promote or demote fields in a BumpyMatrix object and perform summarization on assays.

## 0.1.3.4 (2022-08-11)
- major improvements in JSON validation/conversion logic 
  * provide info about JSON schemes via env variables
  * convert and validate MAE summary next to the SE experiments

## 0.1.3.3 (2022-07-25)
- move json validation/convertion logic from gDRelastic

## 0.1.3.2 (2022-07-11)
- standardize MAE using default gDR identifiers

## 0.1.3.1 (2022-06-29)
- remove adding integer identifiers at the end of colnames/rownames

## 0.1.3.0 (2022-06-02)
- release 1.3.0

## 0.1.0.48 (2022-05-27)
- correct recognition of empty SE

## 0.1.0.47 (2022-05-25)
- remove redundant validation of rownames in SE

## 0.1.0.46 (2022-05-03)
- update prettify function
- fix hardcoded identifiers in validated_SE

## 0.1.0.45 (2022-04-29)
- avoid using `grep` for getting cotreatment identifiers

## 0.1.0.44 (2022-04-26)
- switched from unnamed to named vector of experiment groups for `single-agent`
- fix the logic in validating single-agent experiments

## 0.1.0.43 (2022-04-13)
- set r2 value to NA for invalid and 0 for constant fits

## 0.1.0.42 (2022-04-08)
- added identifier descriptions

## 0.1.0.41 (2022-04-08)
- fix wrong order of elements in rownames in SE

## 0.1.0.40 (2022-04-07)
- add cap_ic50 function

## 0.1.0.39 (2022-03-31)
- extend possible `Barcode` identifiers

## 0.1.0.38 (2022-03-30)
- fix hardcoded identifiers in validate SE

## 0.1.0.37 (2022-03-28)
- added a space between two-word cotreatment identifiers

## 0.1.0.36 (2022-03-24)
- remove all R CMD check warnings

## 0.1.0.35 (2022-03-22)
- change prettify functions to not substitute metadata

## 0.1.0.34 (2022-03-21)
- add helper function for MAE/experiments

## 0.1.0.33 (2022-03-19)
- move constant fit warning

## 0.1.0.32 (2022-03-18)
- add .calculate_complement 

## 0.1.0.31 (2022-03-16)
- move duplication warning 

## 0.1.0.30 (2022-03-14)
- add support for identifier validation

## 0.1.0.29 (2022-03-03)
- add getter and setter for experiment_raw_data

## 0.1.0.28 (2022-02-16)
- refactor identifier `drugname` to `drug_name`

## 0.1.0.27 (2022-02-10)
- add `MAEpply` function

## 0.1.0.26 (2022-02-01)
- update documentation

## 0.1.0.25 (2022-01-31)
- support SE() init for different BioC versions

## 0.1.0.24 (2022-01-25)
- standardize/improve CI

## 0.1.0.23 (2022-01-25)
- switch unit tests from SE to MAE from `gDRtestData`

## 0.1.0.22 (2022-01-07)
- update assert in `validate_MAE`

## 0.1.0.21 (2022-01-07)
- update SE-related functions of gDRutils to support MAE

## 0.1.0.20 (2021-12-21)
- move combo-related functions from gDRviz
- add unit tests for combo-related functions

## 0.1.0.19 (2021-11-04)
- do not create a nested list of identifiers during mergin SE

## 0.1.0.18 (2021-11-01)
- fix issued with new SummarizedExperiment

## 0.1.0.17 (2021-10-22)
- refactor: generalize prediction help functions for fits

## 0.1.0.16 (2021-09-27)
- update validate_SE as per combos

## 0.1.0.15 (2021-09-22)
- feat: use 'GRvalue' and 'RelativeViability' as normalization_types in 'fit_curves'

## 0.1.0.14 (2021-09-15)
- fix bug with getting_SE_identifiers for untreated_tag

## 0.1.0.13 (2021-09-13)
- fix obsolete arguments in `reset_env_identifiers`

## 0.1.0.12 (2021-09-07)
- fix bug with getting identifiers for untreated_tag

## 0.1.0.11 (2021-08-26)
- move NA logic for elements of BumpyMatrix with no data to fit_curve

## 0.1.0.10 (2021-08-25)
- remove additional ordering line in df_to_bm_assay.R

## 0.1.0.9 (2021-08-19)
- update list of available identifiers (2nd and 3rd drug, data_source)
- update the logic for get_identifier and .get_id
- add prettified_identifier()

## 0.1.0.8 (2021-08-13)
- fix bug with wrong order of rows and cols in `df_to_bm_assay`

## 0.1.0.7 (2021-07-06)
- remove hard version equality for pkg deps 

## 0.1.0.6 (2021-06-29)
- add barcode as identifiers and store identifiers within split_SE_components

## 0.1.0.5 (2021-06-25)
- update the logic for CI/CD (repos fetching is now handled with technical user)

## 0.1.0.4 (2021-06-18)
- add concentration and template as additional identifiers

## 0.1.0.3 (2021-06-21)
- refactor logisticFit with error handling

## 0.1.0.2 (2021-06-18)
- remove deprecated functions
- switch from getMetadata to split_SE_components
- refactor df_to_bm_assay

## 0.1.0.1 (2021-06-02)
- upgrade validate_SE by checking if rowData and colData do not have empty strings
- isolate flatten function

## 0.1.0.0 (2021-06-02)
- release 1.0.0

## 0.0.0.49 (2021-06-02)
- fix bug with merging flatten assays

## 0.0.0.48 (2021-05-24)
- refactor validate_SE function

## 0.0.0.47 (2021-05-18)
- add `merge_SE` function

## 0.0.0.46 (2021-05-18)
- fix typo in SE validator function

## 0.0.0.45 (2021-05-18)
- add SE validator function

## 0.0.0.44 (2021-05-13)
- refactor `prettify_flat_metrics` function

## 0.0.0.43 (2021-04-30)
- remove `Metrics_rownames` during flattening data.frame/data.table

## 0.0.0.42 (2021-04-29)
- add `prettify_flat_metrics` function

## 0.0.0.41 (2021-04-29)
- bugfix flattening data.tables

## 0.0.0.40 (2021-04-29)
- bugfix data.table merge in `convert_se_assay_ref_to_dt`

## 0.0.0.39 (2021-04-27)
- add support for flattening data.tables

## 0.0.0.38 (2021-04-21)
- switch from `processing_metadata` to `.internal`

## 0.0.0.37 (2021-04-20)
- add support for getting and setting processing info metadata

## 0.0.0.36 (2021-04-19)
- sort BumpyMatrix created from dt

## 0.0.0.35 (2021-04-14)
- revert metric arguments in 'convert_se_assay_to_dt'

## 0.0.0.34 (2021-04-09)
- add support for getting and setting fit parameter metadata

## 0.0.0.33 (2021-04-07)
- standardize metric header names

## 0.0.0.32 (2021-04-07)
- add more options for returned data (with 'Metrics' assay) in the case of `convert_se_assay_to_dt`
- add 'convert_se_ref_assay_to_dt' function

## 0.0.0.31 (2021-03-30)
- add support for getting and setting metadata on the SummarizedExperiment object

## 0.0.0.30 (2021-03-30)
- add support for getting and setting metadata on the SummarizedExperiment object
- add numeric type assertions and tests for logistic_4parameters

## 0.0.0.29 (2021-03-23)
- remove positional naming dependence on assay_to_dt for (merge_metrics = TRUE) argument

## 0.0.0.28 (2021-03-16)
- minor refactor to assay_to_dt
- deprecated support for assay_to_dt(include_controls = TRUE) argument (not backwards compatible)

## 0.0.0.27 (2021-03-09)
- update .estimate_xc50 and add tests

## 0.0.0.26 (2021-03-04)
- fix bug with wrong class for S3 methods in convert_assay_data_to_dt
- update the documentation

## 0.0.0.25 (2021-03-02)
- fix bug with missing class in assert for matrix

## 0.0.0.24 (2021-02-10)
- modify fit_curves to take in flexible curve_type(s)
- clean-up assay_to_dt

## 0.0.0.23 (2021-02-04)
- added linter

## 0.0.0.22 (2021-02-03)
- move assay_to_dt from gDR to gDRutils
- refactor assay_to_dt to support two assay types
  * list of DFrame(s)
  * BumpyMatrix objects
- move .get_treated_conditions and .get_untreated_conditions from gDR to gDRutils 

## 0.0.0.21 (2021-01-26)
- allow for DFrame as input

## 0.0.0.20 (2021-01-26)
- 's/BumpyMatrix::splitToBumpyMatrix/BumpyMatrix::splitAsBumpyMatrix/'

## 0.0.0.19 (2021-01-20)
- refactor RVGRfits and rename to fit_curves
- make logisticFit function independent of curve_type

## 0.0.0.18 (2021-01-19)
- move df_to_assay.R from gDRcore to gDRutils
- move df_to_bm_assay.R from gDRcore to gDRutils
- update identifiers_list.R

## 0.0.0.15 (2020-12-14)
- refactor get_identifiers and get_headers to be settable and cached

## 0.0.0.14 (2020-12-14)
- totally, finally, and unscrupulously remove dplyr package from gDRutils
- replace IC by RV

## 0.0.0.13 (2020-10-12)
- add CI

## 0.0.0.12 (2020-10-09)
- fix assay_to_df

## 0.0.0.11 (2020-10-09)
- update namespaces

## 0.0.0.10 (2020-10-05)
- minor refactor

## 0.0.0.9 (2020-09-14)
- add small fixes for assays with empty DataFrame

## 0.0.0.8 (2020-09-02)
- update logic as per new db model

## 0.0.0.7 (2020-08-07)
- update variable names as per new db model

## 0.0.0.6 (2020-07-02)
- import pipes from magrittr

## 0.0.0.4 (2020-06-10)
- including the masked field to be able to remove the masked data from averages
