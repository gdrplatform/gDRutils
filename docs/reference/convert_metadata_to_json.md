# Convert experiment metadata to JSON format for elasticsearch indexing.

Convert experiment metadata to JSON format for elasticsearch indexing.

## Usage

``` r
convert_metadata_to_json(se)
```

## Arguments

- se:

  SummarizedExperiment object.

## Value

JSON string capturing experiment metadata.

## Examples

``` r
md <- list(title = "my awesome experiment",
  description = "description of experiment",
  sources = list(list(name = "GeneData_Screener", id = "QCS-12345")))
se <- SummarizedExperiment::SummarizedExperiment(metadata = md)
convert_metadata_to_json(se)
#>  
```
