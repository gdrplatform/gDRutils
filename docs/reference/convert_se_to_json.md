# Convert a SummarizedExperiment object to a JSON document.

Convert a SummarizedExperiment object to a JSON document.

## Usage

``` r
convert_se_to_json(se)
```

## Arguments

- se:

  SummarizedExperiment object.

## Value

String representation of a JSON document.

## Examples

``` r
md <- list(title = "my awesome experiment",
  description = "description of experiment",
  source = list(name = "GeneData_Screener", id = "QCS-12345"))
rdata <- data.table::data.table(
 mydrug = letters,
  mydrugname = letters,
  mydrugmoa = letters,
  Duration = 1)
cdata <- data.table::data.table(mycellline = letters, mycelllinename = letters,
 mycelllinetissue = letters, cellline_ref_div_time = letters)
identifiers <- list(cellline = "mycellline",
                    cellline_name = "mycelllinename",
                    cellline_tissue = "mycelllinetissue",
                    cellline_ref_div_time = "cellline_ref_div_time",
                    drug = "mydrug",
                    drug_name = "mydrugname",
                    drug_moa = "mydrugmoa",
                    duration = "Duration")
se <- SummarizedExperiment::SummarizedExperiment(rowData = rdata,
                                                 colData = cdata)
se <- set_SE_experiment_metadata(se, md)
se <- set_SE_identifiers(se, identifiers)
convert_se_to_json(se)
#> {"title":"my awesome experiment","description":"description of experiment","source":{"name":"GeneData_Screener","id":"QCS-12345"},"drug":["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"],"drug_name":["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"],"drug_moa":["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"],"duration":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"misc_rowdata":{},"cellline":["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"],"cellline_name":["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"],"cellline_tissue":["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"],"cellline_ref_div_time":["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"],"misc_coldata":{}} 
```
