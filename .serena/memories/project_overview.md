---
name: gDRutils project overview
description: High-level overview of the gDRutils R package purpose, structure, and tech stack
type: project
---

# gDRutils

**Purpose**: R utility package for the gDR (drug response) platform. Provides helper functions for:
- Fitting dose-response curves
- Manipulating/converting data between long table and SummarizedExperiment/MAE structures
- Identifier get/set/validation
- Constants and defaults for the gDR platform

**Version**: 1.9.7 (Bioconductor package)
**License**: Artistic-2.0

## Tech Stack
- Language: R (>= 4.2)
- Key dependencies: BumpyMatrix, MultiAssayExperiment, SummarizedExperiment, data.table, drc, checkmate, S4Vectors, qs, jsonlite

## Code Structure
- `R/` — source files (27 files):
  - `utils.R` — large general utilities (~51KB)
  - `fit_curves.R` — dose-response curve fitting (~29KB)
  - `convert_mae_se_assay_to_dt.R` — MAE/SE to data.table conversion (~21KB)
  - `standardize_MAE.R` — MAE standardization (~17KB)
  - `merge_SE.R` — SummarizedExperiment merging (~13KB)
  - `identifiers.R`, `identifiers_list.R` — column identifier management
  - `assay_names.R`, `headers.R`, `headers_list.R` — assay/header naming
  - `combo.R` — combination drug data helpers
  - `concatentate_SEs.R`, `split_SE_components.R` — SE manipulation
  - `json_const.R`, `json_convert.R`, `json_validate.R` — JSON handling
  - `prettify.R` — display formatting
  - `se_metadata.R`, `manage_additional_metadata.R` — metadata management
  - `experiment_validators.R`, `validate_identifiers.R` — validation
  - `flatten.R`, `duplicates.R`, `global_cache.R` — misc utilities
- `tests/testthat/` — test files matching each R source file
- `man/` — roxygen2-generated documentation
- `vignettes/` — package vignettes
- `inst/` — installed files
