---
name: gDRutils suggested commands
description: Key commands for developing, testing, linting, and building the gDRutils R package
type: project
---

# Suggested Commands for gDRutils Development

## Testing
```r
# Run all tests
testthat::test_package("gDRutils")
# or from shell:
Rscript -e 'testthat::test_package("gDRutils")'

# Run a specific test file
Rscript -e 'testthat::test_file("tests/testthat/test-fit_curves.R")'
```

## Linting
```r
lintr::lint_package()
# or from shell:
Rscript -e 'lintr::lint_package()'
```

## Documentation (Roxygen2)
```r
roxygen2::roxygenise()
# or:
devtools::document()
```

## Build & Check
```bash
R CMD build .
R CMD check gDRutils_*.tar.gz
# or:
Rscript -e 'rcmdcheck::rcmdcheck()'
```

## Install locally
```r
devtools::install()
```

## Load for development
```r
devtools::load_all()
```
