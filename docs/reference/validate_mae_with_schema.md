# Validate MAE against a schema.

Validate MAE object against a schema. Currently only SEs are validated
TODO: add mae.json schema and validate full MAE object

## Usage

``` r
validate_mae_with_schema(
  mae,
  schema_package = Sys.getenv("SCHEMA_PACKAGE", "gDRutils"),
  schema_dir_path = Sys.getenv("SCHEMA_DIR_PATH", "schemas"),
  schema = c(se = "se.json", mae = "mae.json")
)
```

## Arguments

- mae:

  MultiAssayExperiment object

- schema_package:

  string name of the package with JSON schema files

- schema_dir_path:

  path to the dir with JSON schema files

- schema:

  named charvec with filenames of schemas to validate against.

## Value

Boolean of whether or not mae is valid

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
validate_mae_with_schema(mae)
#> $mae
#> [1] TRUE
#> attr(,"exit_code")
#> [1] 0
#> 
#> $`experiment:single-agent`
#> [1] FALSE
#> attr(,"error")
#> [1] "JSON validation of data model failed"
#> attr(,"exit_code")
#> [1] 2
#> attr(,"derror")
#>   instancePath schemaPath  keyword missingProperty
#> 1              #/required required           title
#> 2              #/required required     description
#> 3              #/required required experimentalist
#> 4              #/required required         sources
#>                                         message       schema
#> 1           must have required property 'title' title, d....
#> 2     must have required property 'description' title, d....
#> 3 must have required property 'experimentalist' title, d....
#> 4         must have required property 'sources' title, d....
#>                     parentSchema.$schema   parentSchema.title
#> 1 http://json-schema.org/draft-07/schema SummarizedExperiment
#> 2 http://json-schema.org/draft-07/schema SummarizedExperiment
#> 3 http://json-schema.org/draft-07/schema SummarizedExperiment
#> 4 http://json-schema.org/draft-07/schema SummarizedExperiment
#>   parentSchema.description parentSchema.type parentSchema.properties.type
#> 1                                     object                       string
#> 2                                     object                       string
#> 3                                     object                       string
#> 4                                     object                       string
#>   parentSchema.properties.description.type
#> 1                                   string
#> 2                                   string
#> 3                                   string
#> 4                                   string
#>        parentSchema.properties.description.description
#> 1 A free text description of the experiment performed.
#> 2 A free text description of the experiment performed.
#> 3 A free text description of the experiment performed.
#> 4 A free text description of the experiment performed.
#>   parentSchema.properties.experimentalist.type
#> 1                                       string
#> 2                                       string
#> 3                                       string
#> 4                                       string
#>   parentSchema.properties.experimentalist.description
#> 1                 Unix id or name of experimentalist.
#> 2                 Unix id or name of experimentalist.
#> 3                 Unix id or name of experimentalist.
#> 4                 Unix id or name of experimentalist.
#>   parentSchema.properties.experimentalist.examples
#> 1                                            eqlin
#> 2                                            eqlin
#> 3                                            eqlin
#> 4                                            eqlin
#>                   parentSchema.properties.duration.description
#> 1 Number of days cells were treated with perturbation in days.
#> 2 Number of days cells were treated with perturbation in days.
#> 3 Number of days cells were treated with perturbation in days.
#> 4 Number of days cells were treated with perturbation in days.
#>   parentSchema.properties.duration.type parentSchema.properties.duration.type
#> 1                                 array                               integer
#> 2                                 array                               integer
#> 3                                 array                               integer
#> 4                                 array                               integer
#>   parentSchema.properties.cellline.description
#> 1            Identifier of the given cell line
#> 2            Identifier of the given cell line
#> 3            Identifier of the given cell line
#> 4            Identifier of the given cell line
#>   parentSchema.properties.cellline.type parentSchema.properties.cellline.type
#> 1                                 array                                string
#> 2                                 array                                string
#> 3                                 array                                string
#> 4                                 array                                string
#>   parentSchema.properties.cellline_name.description
#> 1                       Name of the given cell line
#> 2                       Name of the given cell line
#> 3                       Name of the given cell line
#> 4                       Name of the given cell line
#>   parentSchema.properties.cellline_name.type
#> 1                                      array
#> 2                                      array
#> 3                                      array
#> 4                                      array
#>   parentSchema.properties.cellline_name.type
#> 1                                     string
#> 2                                     string
#> 3                                     string
#> 4                                     string
#>   parentSchema.properties.cellline_tissue.description
#> 1                Tissue origin of the given cell line
#> 2                Tissue origin of the given cell line
#> 3                Tissue origin of the given cell line
#> 4                Tissue origin of the given cell line
#>   parentSchema.properties.cellline_tissue.type
#> 1                                        array
#> 2                                        array
#> 3                                        array
#> 4                                        array
#>   parentSchema.properties.cellline_tissue.type
#> 1                                       string
#> 2                                       string
#> 3                                       string
#> 4                                       string
#>   parentSchema.properties.cellline_ref_div_time.description
#> 1            Reference division time of the given cell line
#> 2            Reference division time of the given cell line
#> 3            Reference division time of the given cell line
#> 4            Reference division time of the given cell line
#>   parentSchema.properties.cellline_ref_div_time.type
#> 1                                              array
#> 2                                              array
#> 3                                              array
#> 4                                              array
#>   parentSchema.properties.cellline_ref_div_time.items.type
#> 1                                             integer,....
#> 2                                             integer,....
#> 3                                             integer,....
#> 4                                             integer,....
#>   parentSchema.properties.cellline_ref_div_time.items.pattern
#> 1                                                   [0-9]+|NA
#> 2                                                   [0-9]+|NA
#> 3                                                   [0-9]+|NA
#> 4                                                   [0-9]+|NA
#>   parentSchema.properties.drug.description parentSchema.properties.drug.type
#> 1             Identifier of the given drug                             array
#> 2             Identifier of the given drug                             array
#> 3             Identifier of the given drug                             array
#> 4             Identifier of the given drug                             array
#>   parentSchema.properties.drug.type
#> 1                            string
#> 2                            string
#> 3                            string
#> 4                            string
#>   parentSchema.properties.drug_name.description
#> 1                        Name of the given drug
#> 2                        Name of the given drug
#> 3                        Name of the given drug
#> 4                        Name of the given drug
#>   parentSchema.properties.drug_name.type parentSchema.properties.drug_name.type
#> 1                                  array                                 string
#> 2                                  array                                 string
#> 3                                  array                                 string
#> 4                                  array                                 string
#>   parentSchema.properties.drug_moa.description
#> 1        Mechanism of action of the given drug
#> 2        Mechanism of action of the given drug
#> 3        Mechanism of action of the given drug
#> 4        Mechanism of action of the given drug
#>   parentSchema.properties.drug_moa.type parentSchema.properties.drug_moa.type
#> 1                                 array                                string
#> 2                                 array                                string
#> 3                                 array                                string
#> 4                                 array                                string
#>   parentSchema.properties.sources.type
#> 1                                array
#> 2                                array
#> 3                                array
#> 4                                array
#>   parentSchema.properties.sources.description
#> 1                                 Data source
#> 2                                 Data source
#> 3                                 Data source
#> 4                                 Data source
#>   parentSchema.properties.sources.items.type
#> 1                                     object
#> 2                                     object
#> 3                                     object
#> 4                                     object
#>   parentSchema.properties.sources.items.required
#> 1                                       name, id
#> 2                                       name, id
#> 3                                       name, id
#> 4                                       name, id
#>   parentSchema.properties.sources.items.additionalProperties
#> 1                                                      FALSE
#> 2                                                      FALSE
#> 3                                                      FALSE
#> 4                                                      FALSE
#>   parentSchema.properties.sources.items.properties.name.type
#> 1                                                     string
#> 2                                                     string
#> 3                                                     string
#> 4                                                     string
#>   parentSchema.properties.sources.items.properties.name.examples
#> 1                                                   GeneData....
#> 2                                                   GeneData....
#> 3                                                   GeneData....
#> 4                                                   GeneData....
#>   parentSchema.properties.sources.items.properties.id.type
#> 1                                                   string
#> 2                                                   string
#> 3                                                   string
#> 4                                                   string
#>   parentSchema.properties.sources.items.properties.id.description
#> 1                                   An identifier for the source.
#> 2                                   An identifier for the source.
#> 3                                   An identifier for the source.
#> 4                                   An identifier for the source.
#>   parentSchema.properties.sources.items.properties.id.examples
#> 1                                                    IDS-12345
#> 2                                                    IDS-12345
#> 3                                                    IDS-12345
#> 4                                                    IDS-12345
#>   parentSchema.required    data.drug data.drug_name data.drug_moa data.duration
#> 1          title, d.... G00002, ....   drug_002....  moa_A, m....  72, 72, ....
#> 2          title, d.... G00002, ....   drug_002....  moa_A, m....  72, 72, ....
#> 3          title, d.... G00002, ....   drug_002....  moa_A, m....  72, 72, ....
#> 4          title, d.... G00002, ....   drug_002....  moa_A, m....  72, 72, ....
#>   data.cellline data.cellline_name data.cellline_tissue
#> 1  CL00011,....       cellline....         tissue_x....
#> 2  CL00011,....       cellline....         tissue_x....
#> 3  CL00011,....       cellline....         tissue_x....
#> 4  CL00011,....       cellline....         tissue_x....
#>   data.cellline_ref_div_time dataPath
#> 1               26, 30, ....         
#> 2               26, 30, ....         
#> 3               26, 30, ....         
#> 4               26, 30, ....         
#> 
```
