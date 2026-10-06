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
#> [1] "2026-09-17"
#> 
#> $.internal$session_info
#> R version 4.4.0 (2024-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 22.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Etc/UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] gDRcore_1.11.14    gDRtestData_1.11.7 testthat_3.2.1.1  
#> 
#> loaded via a namespace (and not attached):
#>  [1] farver_2.1.2                fastmap_1.2.0              
#>  [3] BumpyMatrix_1.12.0          TH.data_1.1-2              
#>  [5] promises_1.5.0              stringfish_0.19.2          
#>  [7] digest_0.6.39               mime_0.13                  
#>  [9] lifecycle_1.0.5             ellipsis_0.3.2             
#> [11] gDRutils_1.11.9             survival_3.5-8             
#> [13] magrittr_2.0.5              compiler_4.4.0             
#> [15] rlang_1.3.0                 drc_3.0-1                  
#> [17] tools_4.4.0                 plotrix_3.8-14             
#> [19] data.table_1.18.4           lambda.r_1.2.4             
#> [21] S4Arrays_1.4.1              htmlwidgets_1.6.4          
#> [23] pkgbuild_1.4.8              DelayedArray_0.30.1        
#> [25] RColorBrewer_1.1-3          pkgload_1.4.0              
#> [27] multcomp_1.4-31             abind_1.4-8                
#> [29] BiocParallel_1.38.0         miniUI_0.1.2               
#> [31] withr_3.0.3                 purrr_1.2.2                
#> [33] BiocGenerics_0.50.0         desc_1.4.3                 
#> [35] grid_4.4.0                  stats4_4.4.0               
#> [37] urlchecker_1.0.1            profvis_0.3.8              
#> [39] xtable_1.8-8                scales_1.4.0               
#> [41] gtools_3.9.5                MASS_7.3-60.2              
#> [43] MultiAssayExperiment_1.30.2 SummarizedExperiment_1.34.0
#> [45] cli_3.6.6                   mvtnorm_1.2-5              
#> [47] crayon_1.5.3                remotes_2.4.2.9000         
#> [49] RcppParallel_6.2.0          otel_0.2.0                 
#> [51] rstudioapi_0.18.0           httr_1.4.8                 
#> [53] sessioninfo_1.2.3           cachem_1.1.0               
#> [55] stringr_1.6.0               zlibbioc_1.50.0            
#> [57] splines_4.4.0               parallel_4.4.0             
#> [59] XVector_0.44.0              formatR_1.14               
#> [61] matrixStats_1.5.0           vctrs_0.7.3                
#> [63] devtools_2.4.5              Matrix_1.7-0               
#> [65] sandwich_3.1-3              jsonlite_2.0.0             
#> [67] carData_3.0-6               car_3.1-5                  
#> [69] IRanges_2.38.0              S4Vectors_0.42.0           
#> [71] Formula_1.2-6               glue_1.8.1                 
#> [73] codetools_0.2-20            stringi_1.8.9              
#> [75] futile.logger_1.4.9         later_1.4.8                
#> [77] GenomeInfoDb_1.40.1         GenomicRanges_1.56.1       
#> [79] UCSC.utils_1.0.0            htmltools_0.5.9            
#> [81] brio_1.1.5                  GenomeInfoDbData_1.2.12    
#> [83] R6_2.6.1                    rprojroot_2.1.1            
#> [85] Biobase_2.64.0              shiny_1.13.0               
#> [87] lattice_0.22-6              futile.options_1.0.1       
#> [89] backports_1.5.1             memoise_2.0.1              
#> [91] httpuv_1.6.17               Rcpp_1.1.2                 
#> [93] SparseArray_1.4.8           checkmate_2.3.4            
#> [95] qs2_0.2.2                   fs_2.1.0                   
#> [97] MatrixGenerics_1.16.0       zoo_1.9-0                  
#> [99] usethis_2.2.3              
#> 
#> 
```
