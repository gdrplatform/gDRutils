# Changes to v.1.0.1
- upgrade validate_SE by checking if rowData and colData do not have empty strings

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

