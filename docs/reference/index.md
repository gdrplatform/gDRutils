# Package index

## SE operators

- [`set_SE_fit_parameters()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`set_SE_processing_metadata()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`set_SE_keys()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`set_SE_experiment_metadata()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`set_SE_experiment_raw_data()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`get_SE_fit_parameters()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`get_SE_processing_metadata()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`get_SE_experiment_raw_data()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`get_SE_experiment_metadata()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`get_SE_keys()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`get_SE_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  [`set_SE_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
  : Get and set metadata for parameters on a SummarizedExperiment
  object.

- [`aggregate_assay()`](https://gdrplatform.github.io/gDRstyle/reference/aggregate_assay.md)
  :

  Aggregate a `BumpyMatrix` assay by a given aggregation function.

- [`demote_fields()`](https://gdrplatform.github.io/gDRstyle/reference/demote_fields.md)
  :

  Demote a metadata field in the `rowData` or `colData` of a
  `SummarizedExperiment` object to a nested field of a `BumpyMatrix`
  assay.

- [`get_MAE_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/get_MAE_identifiers.md)
  : get_MAE_identifiers

- [`identify_unique_se_metadata_fields()`](https://gdrplatform.github.io/gDRstyle/reference/identify_unique_se_metadata_fields.md)
  :

  Identify unique metadata fields from a list of `SummarizedExperiment`s

- [`merge_MAE()`](https://gdrplatform.github.io/gDRstyle/reference/merge_MAE.md)
  : Merge multiple MultiAssayExperiment objects

- [`merge_SE()`](https://gdrplatform.github.io/gDRstyle/reference/merge_SE.md)
  : Merge multiple Summarized Experiments

- [`merge_assay()`](https://gdrplatform.github.io/gDRstyle/reference/merge_assay.md)
  : Merge assay data

- [`merge_metadata()`](https://gdrplatform.github.io/gDRstyle/reference/merge_metadata.md)
  : Merge metadata

- [`promote_fields()`](https://gdrplatform.github.io/gDRstyle/reference/promote_fields.md)
  :

  Promote a nested field to be represented as a metadata field of the
  `SummarizedExperiment` as either the `rowData` or `colData`.

- [`split_SE_components()`](https://gdrplatform.github.io/gDRstyle/reference/split_SE_components.md)
  : split_SE_components

## Converters

- [`convert_mae_assay_to_dt()`](https://gdrplatform.github.io/gDRstyle/reference/convert_mae_assay_to_dt.md)
  : Convert a MultiAssayExperiment assay to a long data.table
- [`convert_se_assay_to_custom_dt()`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_custom_dt.md)
  : Convert a SummarizedExperiment assay to a long data.table and
  conduct some post processing steps
- [`convert_se_assay_to_dt()`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.md)
  : Convert a SummarizedExperiment assay to a long data.table
- [`df_to_bm_assay()`](https://gdrplatform.github.io/gDRstyle/reference/df_to_bm_assay.md)
  : df_to_bm_assay
- [`flatten()`](https://gdrplatform.github.io/gDRstyle/reference/flatten.md)
  : Flatten a table

## Identifiers helpers/operators

- [`get_default_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/get_default_identifiers.md)
  : Get gDR default identifiers required for downstream analysis.
- [`get_expect_one_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/get_expect_one_identifiers.md)
  : Get identifiers that expect only one value for each identifier.
- [`get_identifiers_dt()`](https://gdrplatform.github.io/gDRstyle/reference/get_identifiers_dt.md)
  : Get descriptions for identifiers
- [`get_idfs_synonyms()`](https://gdrplatform.github.io/gDRstyle/reference/get_idfs_synonyms.md)
  : Get gDR synonyms for the identifiers
- [`get_required_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/get_required_identifiers.md)
  : Get identifiers required for downstream analysis.
- [`get_header()`](https://gdrplatform.github.io/gDRstyle/reference/headers.md)
  : Get or reset headers for one or all header field(s) respectively
- [`get_env_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.md)
  [`get_prettified_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.md)
  [`set_env_identifier()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.md)
  [`reset_env_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.md)
  : Get, set, or reset identifiers for one or all identifier field(s)
- [`prettify_flat_metrics()`](https://gdrplatform.github.io/gDRstyle/reference/prettify_flat_metrics.md)
  : Prettify metric names in flat 'Metrics' assay
- [`update_env_idfs_from_mae()`](https://gdrplatform.github.io/gDRstyle/reference/update_env_idfs_from_mae.md)
  : Update environment identifiers from MAE object identifiers
- [`update_idfs_synonyms()`](https://gdrplatform.github.io/gDRstyle/reference/update_idfs_synonyms.md)
  : Update gDR synonyms for the identifier
- [`validate_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/validate_identifiers.md)
  : Check that specified identifier values exist in the data.

## Experiment helpers

- [`get_experiment_groups()`](https://gdrplatform.github.io/gDRstyle/reference/get_experiment_groups.md)
  : get_experiment_groups
- [`get_supported_experiments()`](https://gdrplatform.github.io/gDRstyle/reference/get_supported_experiments.md)
  : get_supported_experiments
- [`validate_MAE()`](https://gdrplatform.github.io/gDRstyle/reference/validate_MAE.md)
  : Validate MultiAssayExperiment object
- [`validate_SE()`](https://gdrplatform.github.io/gDRstyle/reference/validate_SE.md)
  : Validate SummarizedExperiment object
- [`validate_dimnames()`](https://gdrplatform.github.io/gDRstyle/reference/validate_dimnames.md)
  : Validate dimnames
- [`validate_se_assay_name()`](https://gdrplatform.github.io/gDRstyle/reference/validate_se_assay_name.md)
  : Check whether or not an assay exists in a SummarizedExperiment
  object.

## Assay names

- [`get_assay_names()`](https://gdrplatform.github.io/gDRstyle/reference/get_assay_names.md)
  :

  get assay names of the given se/dataset fetch the data from the se if
  provided as metadata use predefined values from `get_env_assay_names`
  otherwise

- [`get_combo_assay_names()`](https://gdrplatform.github.io/gDRstyle/reference/get_combo_assay_names.md)
  : get names of combo assays

- [`get_combo_base_assay_names()`](https://gdrplatform.github.io/gDRstyle/reference/get_combo_base_assay_names.md)
  : get names of combo base assays

- [`get_combo_score_assay_names()`](https://gdrplatform.github.io/gDRstyle/reference/get_combo_score_assay_names.md)
  : get names of combo score assays

- [`get_env_assay_names()`](https://gdrplatform.github.io/gDRstyle/reference/get_env_assay_names.md)
  : get default assay names for the specified filters, i.e. set of assay
  types, assay groups and assay data types

## Combination data

- [`convert_combo_data_to_dt()`](https://gdrplatform.github.io/gDRstyle/reference/convert_combo_data_to_dt.md)
  : convert combo assays from SummarizedExperiments to the list of
  data.tables

- [`convert_combo_field_to_assay()`](https://gdrplatform.github.io/gDRstyle/reference/convert_combo_field_to_assay.md)
  : get combo assay names based on the field name

- [`define_matrix_grid_positions()`](https://gdrplatform.github.io/gDRstyle/reference/define_matrix_grid_positions.md)
  : Define matrix grid positions

- [`get_additional_variables()`](https://gdrplatform.github.io/gDRstyle/reference/get_additional_variables.md)
  : Identify and return additional variables in list of dt

- [`get_combo_excess_field_names()`](https://gdrplatform.github.io/gDRstyle/reference/get_combo_excess_field_names.md)
  : get names of combo excess fields

- [`get_combo_score_field_names()`](https://gdrplatform.github.io/gDRstyle/reference/get_combo_score_field_names.md)
  : get names of combo score fields

- [`has_single_codrug_data()`](https://gdrplatform.github.io/gDRstyle/reference/has_single_codrug_data.md)
  : Has Single Codrug Data

- [`has_valid_codrug_data()`](https://gdrplatform.github.io/gDRstyle/reference/has_valid_codrug_data.md)
  : Has Valid Codrug Data

- [`is_combo_data()`](https://gdrplatform.github.io/gDRstyle/reference/is_combo_data.md)
  :

  Checks if `se` is combo dataset.

- [`remove_codrug_data()`](https://gdrplatform.github.io/gDRstyle/reference/remove_codrug_data.md)
  : Remove Codrug Data

- [`round_concentration()`](https://gdrplatform.github.io/gDRstyle/reference/round_concentration.md)
  : Round concentration to ndigit significant digits

## Fit curves

- [`cap_xc50()`](https://gdrplatform.github.io/gDRstyle/reference/cap_xc50.md)
  : Cap XC50 value.
- [`.set_invalid_fit_params()`](https://gdrplatform.github.io/gDRstyle/reference/dot-set_invalid_fit_params.md)
  : Set fit parameters for an invalid fit.
- [`fit_curves()`](https://gdrplatform.github.io/gDRstyle/reference/fit_curves.md)
  : Fit curves
- [`logisticFit()`](https://gdrplatform.github.io/gDRstyle/reference/logisticFit.md)
  : Logistic fit
- [`predict_conc_from_efficacy()`](https://gdrplatform.github.io/gDRstyle/reference/predict_conc_from_efficacy.md)
  : Predict a concentration for a given efficacy with fit parameters.
- [`predict_efficacy_from_conc()`](https://gdrplatform.github.io/gDRstyle/reference/predict_efficacy_from_conc.md)
  : Predict efficacy values given fit parameters and a concentration.
- [`predict_smooth_from_combo()`](https://gdrplatform.github.io/gDRstyle/reference/predict_smooth_from_combo.md)
  : Predict a smoothed response for a drug combination
- [`set_constant_fit_params()`](https://gdrplatform.github.io/gDRstyle/reference/set_constant_fit_params.md)
  : Set fit parameters for a constant fit.

## JSON conversion

- [`convert_colData_to_json()`](https://gdrplatform.github.io/gDRstyle/reference/convert_colData_to_json.md)
  : Convert colData to JSON
- [`convert_mae_to_json()`](https://gdrplatform.github.io/gDRstyle/reference/convert_mae_to_json.md)
  : Create JSON document.
- [`convert_metadata_to_json()`](https://gdrplatform.github.io/gDRstyle/reference/convert_metadata_to_json.md)
  : Convert experiment metadata to JSON format for elasticsearch
  indexing.
- [`convert_rowData_to_json()`](https://gdrplatform.github.io/gDRstyle/reference/convert_rowData_to_json.md)
  : Convert rowData to JSON
- [`convert_se_to_json()`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_to_json.md)
  : Convert a SummarizedExperiment object to a JSON document.
- [`strip_first_and_last_char()`](https://gdrplatform.github.io/gDRstyle/reference/strip_first_and_last_char.md)
  : String first and last characters of a string.
- [`validate_mae_with_schema()`](https://gdrplatform.github.io/gDRstyle/reference/validate_mae_with_schema.md)
  : Validate MAE against a schema.

## JSON validation

- [`validate_json()`](https://gdrplatform.github.io/gDRstyle/reference/validate_json.md)
  : Validate JSON against a schema.

## JSON const getters

- [`get_isobologram_columns()`](https://gdrplatform.github.io/gDRstyle/reference/get_isobologram_columns.md)
  : Get isobologram column names
- [`get_settings_from_json()`](https://gdrplatform.github.io/gDRstyle/reference/get_settings_from_json.md)
  : Get settings from JSON file In most common scenario the settings are
  stored in JSON file to avoid hardcoding

## Metadata management

- [`addClass()`](https://gdrplatform.github.io/gDRstyle/reference/addClass.md)
  : add arbitrary S3 class to an object
- [`modifyData()`](https://gdrplatform.github.io/gDRstyle/reference/modifyData.md)
  : modify assay with additional data

## Standardize MAE

- [`get_optional_coldata_fields()`](https://gdrplatform.github.io/gDRstyle/reference/get_optional_coldata_fields.md)
  : get optional colData fields
- [`get_optional_rowdata_fields()`](https://gdrplatform.github.io/gDRstyle/reference/get_optional_rowdata_fields.md)
  : get optional rowData fields
- [`refine_coldata()`](https://gdrplatform.github.io/gDRstyle/reference/refine_coldata.md)
  : refine colData
- [`refine_rowdata()`](https://gdrplatform.github.io/gDRstyle/reference/refine_rowdata.md)
  : refine rowData
- [`rename_DFrame()`](https://gdrplatform.github.io/gDRstyle/reference/rename_DFrame.md)
  : Rename DFrame
- [`rename_bumpy()`](https://gdrplatform.github.io/gDRstyle/reference/rename_bumpy.md)
  : Rename BumpyMatrix
- [`set_unique_cl_names()`](https://gdrplatform.github.io/gDRstyle/reference/set_unique_cl_names.md)
  : Set Unique Parental Identifiers
- [`set_unique_cl_names_dt()`](https://gdrplatform.github.io/gDRstyle/reference/set_unique_cl_names_dt.md)
  : Set unique primary cell line identifiers in the table
- [`set_unique_drug_names()`](https://gdrplatform.github.io/gDRstyle/reference/set_unique_drug_names.md)
  : Set Unique Drug Names
- [`set_unique_drug_names_dt()`](https://gdrplatform.github.io/gDRstyle/reference/set_unique_drug_names_dt.md)
  : Set unique primary drug identifiers in the table
- [`set_unique_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/set_unique_identifiers.md)
  : Set Unique Identifiers in MultiAssayExperiment
- [`set_unique_names_dt()`](https://gdrplatform.github.io/gDRstyle/reference/set_unique_names_dt.md)
  : Set unique primary identifiers in the data.frame-like objects
- [`standardize_mae()`](https://gdrplatform.github.io/gDRstyle/reference/standardize_mae.md)
  : Standardize MAE by switching from custom identifiers into
  gDR-default
- [`standardize_se()`](https://gdrplatform.github.io/gDRstyle/reference/standardize_se.md)
  : Standardize SE by switching from custom identifiers into gDR-default

## Duplicates

- [`get_assay_dt_duplicated_rows()`](https://gdrplatform.github.io/gDRstyle/reference/get_assay_dt_duplicated_rows.md)
  : Helper function to find duplicated rows in assay data
- [`get_assay_req_uniq_cols()`](https://gdrplatform.github.io/gDRstyle/reference/get_assay_req_uniq_cols.md)
  : get columns in the assay data required to have unique data
- [`get_duplicated_rows()`](https://gdrplatform.github.io/gDRstyle/reference/get_duplicated_rows.md)
  : Helper function to find duplicated rows
- [`has_assay_dt_duplicated_rows()`](https://gdrplatform.github.io/gDRstyle/reference/has_assay_dt_duplicated_rows.md)
  : check if assay data contains duplicated data
- [`has_dt_duplicated_rows()`](https://gdrplatform.github.io/gDRstyle/reference/has_dt_duplicated_rows.md)
  : check if data.table contains duplicated data
- [`throw_msg_if_duplicates()`](https://gdrplatform.github.io/gDRstyle/reference/throw_msg_if_duplicates.md)
  : throw message if assay data.table contains duplicated rows

## Utils

- [`MAEpply()`](https://gdrplatform.github.io/gDRstyle/reference/MAEpply.md)
  : Lapply through all the experiments in MultiAssayExperiment object
- [`apply_bumpy_function()`](https://gdrplatform.github.io/gDRstyle/reference/apply_bumpy_function.md)
  : Apply a function to every element of a bumpy matrix.
- [`assert_choices()`](https://gdrplatform.github.io/gDRstyle/reference/assert_choices.md)
  : assert choices
- [`average_biological_replicates_dt()`](https://gdrplatform.github.io/gDRstyle/reference/average_biological_replicates_dt.md)
  : Average biological replicates on the data table side.
- [`calc_sd()`](https://gdrplatform.github.io/gDRstyle/reference/calc_sd.md)
  : Calculate Standard Deviation or Return Zero
- [`cap_assay_infinities()`](https://gdrplatform.github.io/gDRstyle/reference/cap_assay_infinities.md)
  : Cap infinity values (Inf, -Inf) in the assay data
- [`.standardize_conc()`](https://gdrplatform.github.io/gDRstyle/reference/dot-standardize_conc.md)
  : Standardize concentration values.
- [`extend_normalization_type_name()`](https://gdrplatform.github.io/gDRstyle/reference/extend_normalization_type_name.md)
  : extend abbreviated normalization type
- [`geometric_mean()`](https://gdrplatform.github.io/gDRstyle/reference/geometric_mean.md)
  : Geometric mean
- [`get_env_var()`](https://gdrplatform.github.io/gDRstyle/reference/get_env_var.md)
  : safe wrapper of Sys.getenv()
- [`get_gDR_session_info()`](https://gdrplatform.github.io/gDRstyle/reference/get_gDR_session_info.md)
  : get gDR package and their version installed in the environment
- [`get_non_empty_assays()`](https://gdrplatform.github.io/gDRstyle/reference/get_non_empty_assays.md)
  : get_non_empty_assays
- [`get_synthetic_data()`](https://gdrplatform.github.io/gDRstyle/reference/get_synthetic_data.md)
  : Get synthetic data from gDRtestData package
- [`is_any_exp_empty()`](https://gdrplatform.github.io/gDRstyle/reference/is_any_exp_empty.md)
  : is_any_exp_empty
- [`is_exp_empty()`](https://gdrplatform.github.io/gDRstyle/reference/is_exp_empty.md)
  : is_exp_empty
- [`is_mae_empty()`](https://gdrplatform.github.io/gDRstyle/reference/is_mae_empty.md)
  : is_mae_empty
- [`loop()`](https://gdrplatform.github.io/gDRstyle/reference/loop.md) :
  Conditional lapply with optional batch processing.
- [`map_conc_to_standardized_conc()`](https://gdrplatform.github.io/gDRstyle/reference/map_conc_to_standardized_conc.md)
  : Create a mapping of concentrations to standardized concentrations.
- [`mcolData()`](https://gdrplatform.github.io/gDRstyle/reference/mcolData.md)
  : mcolData
- [`mrowData()`](https://gdrplatform.github.io/gDRstyle/reference/mrowData.md)
  : mrowData
- [`process_batch()`](https://gdrplatform.github.io/gDRstyle/reference/process_batch.md)
  : Process and save a batch of results.
- [`remove_drug_batch()`](https://gdrplatform.github.io/gDRstyle/reference/remove_drug_batch.md)
  : Remove batch substring from drug id
- [`shorten_normalization_type_name()`](https://gdrplatform.github.io/gDRstyle/reference/shorten_normalization_type_name.md)
  : shorten normalization type
- [`split_big_table_for_xlsx()`](https://gdrplatform.github.io/gDRstyle/reference/split_big_table_for_xlsx.md)
  : Split big table

## Test helpers

- [`gen_synthetic_data()`](https://gdrplatform.github.io/gDRstyle/reference/gen_synthetic_data.md)
  : gen_synthetic_data
- [`get_testdata()`](https://gdrplatform.github.io/gDRstyle/reference/get_testdata.md)
  : get_testdata
- [`get_testdata_codilution()`](https://gdrplatform.github.io/gDRstyle/reference/get_testdata_codilution.md)
  : get_testdata_codilution
- [`get_testdata_combo()`](https://gdrplatform.github.io/gDRstyle/reference/get_testdata_combo.md)
  : get_testdata_combo
