# Change to v.1.3.22
- Restore cap_ic50

# Change to v.1.3.21
- Better handling of corner cases for the single-agent fitting

# Change to v.1.3.20
- Add support for splitting normalization data types

# Change to v.1.3.19
- Restore necessary functions

# Change to v.1.3.18
- Add examples of identifiers

# Change to v.1.3.17
- Remove obsolete code

# Change to v.1.3.16
- Tidy code

# Change to v.1.3.15
- Add tests for convert_se_ref_assay_to_dt

# Change to v.1.3.14
- Add helpers for colData/rowData refinement

# Change to v.1.3.13
- Fixed warnings in unit tests

# Change to v.1.3.12
- `R CMD check` returns `0 errors ✓ | 0 warnings ✓ | 0 notes ✓` 

# Change to v.1.3.11
- default apply_bumpy_function parallelize to TRUE 

# Change to v.1.3.10
- add apply_bumpy_function support

# Change to v.1.3.9
- add co-dilution to single-agent group

# Change to v.1.3.8
- Bugfix in validate_SE

# Change to v.1.3.7
- add update_env_idfs_from_mae function

# Change to v.1.3.6
- add helper for dealing with idfs' synonyms

# Change to v.1.3.5
- Promote or demote fields in a BumpyMatrix object and perform summarization on assays.

# Change to v.1.3.4
- Major improvents in JSON validation/conversion logic 
  * provide info about JSON schemas via env variables
  * convert and validate MAE summary next to the SE experiments

# Change to v.1.3.3
- Move json validation/convertion logic from gDRelastic

# Change to v.1.3.2
- Standardize MAE using default gDR identifiers

# Change to v.1.3.1
- Remove adding integer identifiers at the end of colnames/rownames

# Change to v.1.3.0
- Release

# Change to v.1.0.48
- correct recognition of empty SE

# Change to v.1.0.47
- remove redundant validation of rownames in SE

# Change to v.1.0.46
- update prettify function
- fix hardcoded identifiers in validated_SE

# Change to v.1.0.45
- avoid using `grep` for getting cotreatment identifiers

# Change to v.1.0.44
- switched from unnamed to named vector of experiment groups for `single-agent`
- fixed the logic in validating single-agent experiments

# Change to v.1.0.43
- set r2 value to NA for invalid and 0 for constant fits

# Change to v.1.0.42
- added identifier descriptions

# Change to v.1.0.41
- fix wrong order of elements in rownames in SE

# Change to v.1.0.40
- add cap_ic50 function

# Change to v.1.0.39
- extend possible `Barcode` identifiers

# Change to v.1.0.38
- fix hardcoded identifiers in validate SE

# Change to v.1.0.37
- added a space between two-word cotreatment identifiers

# Change to v.1.0.36
- remove all R CMD check warnings

# Change to v.1.0.35
- change prettify functions to not substitute metadata

# Change to v.1.0.34
- add helper function for MAE/experiments

# Change to v.1.0.33
- move constant fit warning

# Change to v.1.0.32
- add .calculate_complement 

# Change to v.1.0.31
- move duplication warning 

# Change to v.1.0.30
- add support for identifier validation

# Change to v.1.0.29
- add getter and setter for experiment_raw_data

# Change to v.1.0.28
- refactor identifier `drugname` to `drug_name`

# Change to v.1.0.27
- add `MAEpply` function

# Change to v.1.0.26
- update documentation

# Change to v.1.0.25
- support SE() init for different BioC versions

# Change to v.1.0.24
- standardize/improve CI

# Change to v.1.0.23
- switch unit tests from SE to MAE from `gDRtestData`

# Change to v.1.0.21
- update assert in `validate_MAE`

# Change to v.1.0.21
- update SE-related functions of gDRutils to support MAE

# Change to v.1.0.20
- move combo-related functions from gDRviz
- add unit tests for combo-related functions

# Change to v.1.0.19
- do not create a nested list of identifiers during mergin SE

# Changes to v.1.0.18
- fixed issued with new SummarizedExperiment

# Changes to v.1.0.17
- refactor: generalize prediction help functions for fits

# Changes to v.1.0.16
- update validate_SE as per combos

# Changes to v.1.0.15
- feat: use 'GRvalue' and 'RelativeViability' as normalization_types in 'fit_curves'

# Changes to v.1.0.14
- fix bug with getting_SE_identifiers for untreated_tag

# Changes to v.1.0.13
- fix obsolete arguments in `reset_env_identifiers`

# Changes to v.1.0.12
- fix bug with getting identifiers for untreated_tag

# Changes to v.1.0.11
- move NA logic for elements of BumpyMatrix with no data to fit_curve

# Change to v.1.0.10
- remove additional ordering line in df_to_bm_assay.R

# Changes to v.1.0.9
- update list of available identifiers (2nd and 3rd drug, data_source)
- update the logic for get_identifier and .get_id
- add prettified_identifier()

# Changes to v.1.0.8
- fix bug with wrong order of rows and cols in `df_to_bm_assay`

# Changes to v.1.0.7
- remove hard version equality for pkg deps 

# Changes to v.1.0.6
- add barcode as identifiers and store identifiers within split_SE_components

# Changes to v.1.0.5
- update the logic for CI/CD (repos fetching is now handled with technical user)

# Changes to v.1.0.4
- add concentration and template as additional identifiers

# Changes to v.1.0.3
- refactor logisticFit with error handling

# Changes to v.1.0.2
- remove deprecated functions
- switch from getMetadata to split_SE_components
- refactor df_to_bm_assay

# Changes to v.1.0.1
- upgrade validate_SE by checking if rowData and colData do not have empty strings
- isolate flatten function

# Changes to v.1.0.0
- release 1.0.0

# Changes to v.0.0.49
- fix bug with merging flatten assays

# Changes to v.0.0.48
- refactor validate_SE function

# Changes to v.0.0.47
- add `merge_SE` function

# Changes to v.0.0.46
- fix typo in SE validator function

# Changes to v.0.0.45
- add SE validator function

# Changes to v.0.0.44
- refactor `prettify_flat_metrics` function

# Changes to v.0.0.43
- remove `Metrics_rownames` during flattening data.frame/data.table

# Changes to v.0.0.42
- add `prettify_flat_metrics` function

# Changes to v.0.0.41
- bugfix flattening data.tables

# Changes to v.0.0.40
- bugfix data.table merge in `convert_se_assay_ref_to_dt`

# Changes to v.0.0.39
- add support for flattening data.tables

# Changes to v.0.0.38
- switch from `processing_metadata` to `.internal`

# Changes to v.0.0.37
- add support for getting and setting processing info metadata

# Changes to v.0.0.36
- sort BumpyMatrix created from dt

# Changes to v.0.0.35
- revert metric arguments in 'convert_se_assay_to_dt'

# Changes to v.0.0.34
- add support for getting and setting fit parameter metadata

# Changes to v.0.0.33
- standardize metric header names

# Changes to v.0.0.32
- add more options for returned data (with 'Metrics' assay) in the case of 'convert_se_assay_to_dt'
- add 'convert_se_ref_assay_to_dt' function

# Changes to v.0.0.31
- add support for getting and setting metadata on the SummarizedExperiment object

# Changes to v.0.0.30
- add support for getting and setting metadata on the SummarizedExperiment object
- add numeric type assertions and tests for logistic_4parameters

# Changes to v.0.0.29
- remove positional naming dependence on assay_to_dt for (merge_metrics = TRUE) argument

# Changes to v.0.0.28
- minor refactor to assay_to_dt
- deprecated support for assay_to_dt(include_controls = TRUE) argument
(not backwards compatible)

# Changes to v.0.0.27
- update .estimate_xc50 and add tests

# Changes to v.0.0.26
- fix bug with wrong class for S3 methods in convert_assay_data_to_dt
- update the documentation

# Changes to v.0.0.25
- fix bug with missing class in assert for matrix

# Changes to v.0.0.24
- modify fit_curves to take in flexible curve_type(s)
- clean-up assay_to_dt

# Changes to v.0.0.23
- added linter

# Changes to v.0.0.22
- move assay_to_dt from gDR to gDRutils
- refactor assay_to_dt to support two assay types
  * list of DFrame(s)
  * BumpyMatrix objects
- move .get_treated_conditions and .get_untreated_conditions from gDR to
    gDRutils 

# Changes to v.0.0.21
- allow for DFrame as input

# Changes to v.0.0.20
- 's/BumpyMatrix::splitToBumpyMatrix/BumpyMatrix::splitAsBumpyMatrix/'

# Changes to v.0.0.19
- refactor RVGRfits and rename to fit_curves
- make logisticFit function independent of curve_type

# Changes to v.0.0.18
- move df_to_assay.R from gDRcore to gDRutils
- move df_to_bm_assay.R from gDRcore to gDRutils
- update identifiers_list.R

# Changes to v.0.0.15
- refactor get_identifiers and get_headers to be settable and cached

# Changes to v.0.0.14
- totally, finally, and unscrupulously remove dplyr package from gDRutils
- replace IC by RV

# Changes to v.0.0.13
- add CI

# Changes to v.0.0.12
- fix assay_to_df

# Changes to v.0.0.11
- update namespaces

# Changes to v.0.0.10
- minor refactor

# Changes to v.0.0.9
- add small fixes for assays with empty DataFrame

# Changes to v.0.0.8
- update logic as per new db model

# Changes to v.0.0.7
- update variable names as per new db model

# Changes to v.0.0.6 - 2020.07.02
- import pipes from magrittr

# Changes to v.0.0.4 - 2020.06.10
- including the masked field to be able to remove the masked data from averages

