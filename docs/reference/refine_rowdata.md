# refine rowData

current improvements done on the rowData as a standardization step:

- set default value for optional rowData fields

## Usage

``` r
refine_rowdata(rd, se, default_v = "Undefined")
```

## Arguments

- rd:

  DataFrame with rowData

- se:

  a SummarizedExperiment object with drug-response data generate by gDR
  pipeline

- default_v:

  string with default value for optional columns in rowData

## Value

refined rowData

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
refine_rowdata(SummarizedExperiment::colData(mae[[1]]), mae[[1]])
#> DataFrame with 10 rows and 5 columns
#>                                        clid CellLineName      Tissue
#>                                 <character>  <character> <character>
#> CL00011_cellline_BA_tissue_x_26     CL00011  cellline_BA    tissue_x
#> CL00012_cellline_CA_tissue_x_30     CL00012  cellline_CA    tissue_x
#> CL00013_cellline_DA_tissue_x_34     CL00013  cellline_DA    tissue_x
#> CL00014_cellline_EA_tissue_x_38     CL00014  cellline_EA    tissue_x
#> CL00015_cellline_FA_tissue_x_42     CL00015  cellline_FA    tissue_x
#> CL00016_cellline_GB_tissue_y_46     CL00016  cellline_GB    tissue_y
#> CL00017_cellline_HB_tissue_y_50     CL00017  cellline_HB    tissue_y
#> CL00018_cellline_IB_tissue_y_54     CL00018  cellline_IB    tissue_y
#> CL00019_cellline_JB_tissue_z_58     CL00019  cellline_JB    tissue_z
#> CL00020_cellline_KB_tissue_z_62     CL00020  cellline_KB    tissue_z
#>                                 ReferenceDivisionTime    drug_moa
#>                                             <numeric> <character>
#> CL00011_cellline_BA_tissue_x_26                    26   Undefined
#> CL00012_cellline_CA_tissue_x_30                    30   Undefined
#> CL00013_cellline_DA_tissue_x_34                    34   Undefined
#> CL00014_cellline_EA_tissue_x_38                    38   Undefined
#> CL00015_cellline_FA_tissue_x_42                    42   Undefined
#> CL00016_cellline_GB_tissue_y_46                    46   Undefined
#> CL00017_cellline_HB_tissue_y_50                    50   Undefined
#> CL00018_cellline_IB_tissue_y_54                    54   Undefined
#> CL00019_cellline_JB_tissue_z_58                    58   Undefined
#> CL00020_cellline_KB_tissue_z_62                    62   Undefined
```
