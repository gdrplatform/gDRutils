---
name: gDRutils code style and conventions
description: R coding style, naming conventions, and documentation patterns used in gDRutils
type: project
---

# Code Style & Conventions

## Naming
- Public functions: `snake_case` (e.g., `fit_curves`, `average_biological_replicates_dt`)
- Private/internal functions: prefixed with `.` (e.g., `.applyLogisticFit`, `.checkNonNaAvgNorm`)
- S4 classes and methods follow Bioconductor conventions

## Documentation
- Roxygen2 with markdown enabled (`Roxygen: list(markdown = TRUE)`)
- All exported functions documented with `@param`, `@return`, `@examples`
- RoxygenNote: 7.3.3

## Style
- lintr used for linting
- ByteCompile: TRUE
- Data manipulation via `data.table` (not dplyr)
- Argument validation via `checkmate`
- Bioconductor-style package (biocViews: Software, Infrastructure)

## Testing
- testthat framework
- One test file per source file (`test-<filename>.R`)
- Setup file at `tests/testthat/setup.R`
