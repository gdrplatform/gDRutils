# refine colData

current improvements done on the colData as a standardization step:

- set default value for optional colData fields

## Usage

``` r
refine_coldata(cd, se, default_v = "Undefined")
```

## Arguments

- cd:

  DataFrame with colData

- se:

  a SummarizedExperiment object with drug-response data generate by gDR
  pipeline

- default_v:

  string with default value for optional columns in colData

## Value

refined colData

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
refine_coldata(SummarizedExperiment::colData(mae[[1]]), mae[[1]])
#> DataFrame with 10 rows and 4 columns
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
#>                                 ReferenceDivisionTime
#>                                             <numeric>
#> CL00011_cellline_BA_tissue_x_26                    26
#> CL00012_cellline_CA_tissue_x_30                    30
#> CL00013_cellline_DA_tissue_x_34                    34
#> CL00014_cellline_EA_tissue_x_38                    38
#> CL00015_cellline_FA_tissue_x_42                    42
#> CL00016_cellline_GB_tissue_y_46                    46
#> CL00017_cellline_HB_tissue_y_50                    50
#> CL00018_cellline_IB_tissue_y_54                    54
#> CL00019_cellline_JB_tissue_z_58                    58
#> CL00020_cellline_KB_tissue_z_62                    62
```
