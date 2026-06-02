# split_SE_components

Divide the columns of an input data.table into treatment metadata,
condition metadata, experiment metadata, and assay data for further
analysis. This will most commonly be used to identify the different
components of a SummarizedExperiment object.

## Usage

``` r
split_SE_components(df_, nested_keys = NULL, combine_on = 1L)
```

## Arguments

- df\_:

  data.table with drug-response data

- nested_keys:

  character vector of keys to exclude from the row or column metadata,
  and to instead nest within an element of the matrix. See details.

- combine_on:

  integer value of `1` or `2`, indicating whether unrecognized columns
  should be combined on row or column respectively. Defaults to `1`.

## Value

named list containing different elements of a SummarizedExperiment; see
details.

## Details

Named list containing the following elements:

- "treatment_md": :

  treatment metadata

- "condition_md": :

  condition metadata

- "data_fields": :

  all data.table column names corresponding to fields nested within a
  BumpyMatrix cell

- "experiment_md": :

  metadata that is constant for all entries of the data.table

- "identifiers_md": :

  key identifier mappings

The `nested_keys` provides the user the opportunity to specify that they
would not like to use that metadata field as a differentiator of the
treatments, and instead, incorporate it into a nested `DataFrame` in the
BumpyMatrix matrix object.

In the event that if any of the `nested_keys` are constant throughout
the whole data.table, they will still be included in the DataFrame of
the BumpyMatrix and not in the experiment_metadata.

Columns within the `df_` will be identified through the following logic:
First, the known data fields and any specified `nested_keys` are
extracted. Following that, known cell and drug metadata fields are
detected, and any remaining columns that represent constant metadata
fields across all rows are extracted. Next, any cell line metadata will
be heuristically extracted. Finally, all remaining columns will be
combined on either the rows or columns as specified by `combine_on`.

## Examples

``` r
split_SE_components(data.table::data.table(clid = "CL1", Gnumber = "DrugA"))
#> $condition_md
#> DataFrame with 1 row and 1 column
#>            clid
#>     <character>
#> CL1         CL1
#> 
#> $treatment_md
#> DataFrame with 1 row and 1 column
#>           Gnumber
#>       <character>
#> DrugA       DrugA
#> 
#> $data_fields
#> character(0)
#> 
#> $experiment_md
#> DataFrame with 1 row and 0 columns
#> 
#> $identifiers_md
#> $identifiers_md$duration
#> [1] "Duration"
#> 
#> $identifiers_md$cellline
#> [1] "clid"
#> 
#> $identifiers_md$cellline_name
#> [1] "CellLineName"
#> 
#> $identifiers_md$cellline_tissue
#> [1] "Tissue"
#> 
#> $identifiers_md$cellline_ref_div_time
#> [1] "ReferenceDivisionTime"
#> 
#> $identifiers_md$cellline_parental_identifier
#> [1] "parental_identifier"
#> 
#> $identifiers_md$cellline_subtype
#> [1] "subtype"
#> 
#> $identifiers_md$drug
#> [1] "Gnumber"
#> 
#> $identifiers_md$drug_name
#> [1] "DrugName"
#> 
#> $identifiers_md$drug_moa
#> [1] "drug_moa"
#> 
#> $identifiers_md$untreated_tag
#> [1] "vehicle"   "untreated"
#> 
#> $identifiers_md$masked_tag
#> [1] "masked"
#> 
#> $identifiers_md$well_position
#> [1] "WellRow"    "WellColumn"
#> 
#> $identifiers_md$concentration
#> [1] "Concentration"
#> 
#> $identifiers_md$template
#> [1] "Template"  "Treatment"
#> 
#> $identifiers_md$barcode
#> [1] "Barcode" "Plate"  
#> 
#> $identifiers_md$drug2
#> [1] "Gnumber_2"
#> 
#> $identifiers_md$drug_name2
#> [1] "DrugName_2"
#> 
#> $identifiers_md$drug_moa2
#> [1] "drug_moa_2"
#> 
#> $identifiers_md$concentration2
#> [1] "Concentration_2"
#> 
#> $identifiers_md$drug3
#> [1] "Gnumber_3"
#> 
#> $identifiers_md$drug_name3
#> [1] "DrugName_3"
#> 
#> $identifiers_md$drug_moa3
#> [1] "drug_moa_3"
#> 
#> $identifiers_md$concentration3
#> [1] "Concentration_3"
#> 
#> $identifiers_md$data_source
#> [1] "data_source"
#> 
#> $identifiers_md$replicate
#> [1] "Replicate"
#> 
#> $identifiers_md$normalization_type
#> [1] "normalization_type"
#> 
#> 
```
