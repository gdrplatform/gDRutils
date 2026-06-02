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
#> $fit_parameters
#> NULL
#> 
#> $.internal
#> $.internal$date_processed
#> [1] "2025-07-29"
#> 
#> $.internal$session_info
#> R version 4.4.1 (2024-06-14)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sonoma 14.6.1
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Warsaw
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] gDRtestData_1.7.1           gDRcore_1.7.7              
#>  [3] gDRutils_1.7.14             gDRinternal_1.0.186        
#>  [5] logExtra_0.1.19             GeneDataScreenR_1.0.7      
#>  [7] testthat_3.2.3              BumpyMatrix_1.14.0         
#>  [9] MultiAssayExperiment_1.32.0 SummarizedExperiment_1.36.0
#> [11] Biobase_2.66.0              GenomicRanges_1.58.0       
#> [13] GenomeInfoDb_1.42.3         IRanges_2.40.1             
#> [15] S4Vectors_0.44.0            BiocGenerics_0.52.0        
#> [17] MatrixGenerics_1.18.1       matrixStats_1.5.0          
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0         
#>   [4] magrittr_2.0.3          TH.data_1.1-3           farver_2.1.2           
#>   [7] rmarkdown_2.29          fs_1.6.6                zlibbioc_1.52.0        
#>  [10] vctrs_0.6.5             memoise_2.0.1           prettydoc_0.4.1        
#>  [13] htmltools_0.5.8.1       S4Arrays_1.6.0          usethis_3.1.0          
#>  [16] plotrix_3.8-4           lambda.r_1.2.4          SparseArray_1.6.2      
#>  [19] Formula_1.2-5           htmlwidgets_1.6.4       desc_1.4.3             
#>  [22] plyr_1.8.9              sandwich_3.1-1          futile.options_1.0.1   
#>  [25] zoo_1.8-14              cachem_1.1.0            commonmark_1.9.5       
#>  [28] dcmad_0.0.6             mime_0.13               lifecycle_1.0.4        
#>  [31] pkgconfig_2.0.3         Matrix_1.7-3            R6_2.6.1               
#>  [34] fastmap_1.2.0           GenomeInfoDbData_1.2.13 shiny_1.10.0           
#>  [37] digest_0.6.37           rprojroot_2.0.4         pkgload_1.4.0          
#>  [40] httr_1.4.7              abind_1.4-8             compiler_4.4.1         
#>  [43] remotes_2.5.0           withr_3.0.2             backports_1.5.0        
#>  [46] BiocParallel_1.40.2     carData_3.0-5           qs_0.27.3              
#>  [49] pkgbuild_1.4.8          MASS_7.3-65             drc_3.0-1              
#>  [52] DelayedArray_0.32.0     sessioninfo_1.2.3       rjson_0.2.23           
#>  [55] gtools_3.9.5            tools_4.4.1             httpuv_1.6.16          
#>  [58] jsonvalidate_1.5.0      glue_1.8.0              promises_1.3.3         
#>  [61] grid_4.4.1              checkmate_2.3.2         generics_0.1.4         
#>  [64] gtable_0.3.6            data.table_1.17.6       RApiSerialize_0.1.4    
#>  [67] xml2_1.3.8              stringfish_0.16.0       car_3.1-3              
#>  [70] XVector_0.46.0          pillar_1.10.2           markdown_2.0           
#>  [73] stringr_1.5.1           later_1.4.2             splines_4.4.1          
#>  [76] dplyr_1.1.4             moments_0.14.1          lattice_0.22-7         
#>  [79] survival_3.8-3          tidyselect_1.2.1        miniUI_0.1.2           
#>  [82] knitr_1.50              futile.logger_1.4.3     xfun_0.52              
#>  [85] devtools_2.4.5          brio_1.1.5              stringi_1.8.7          
#>  [88] UCSC.utils_1.2.0        yaml_2.3.10             evaluate_1.0.4         
#>  [91] codetools_0.2-20        tibble_3.3.0            cli_3.6.5              
#>  [94] RcppParallel_5.1.10     xtable_1.8-4            roxygen2_7.3.2         
#>  [97] log4r_0.4.4             Rcpp_1.0.14             XML_3.99-0.18          
#> [100] parallel_4.4.1          ellipsis_0.3.2          ggplot2_3.5.2          
#> [103] profvis_0.4.0           urlchecker_1.0.1        mvtnorm_1.3-3          
#> [106] scales_1.4.0            purrr_1.0.4             crayon_1.5.3           
#> [109] rlang_1.1.6             multcomp_1.4-28         formatR_1.14           
#> 
#> 
```
