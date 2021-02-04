packagePathForLintr <- commandArgs(trailingOnly = TRUE)[1]
packagePathForLintr <- ifelse(is.na(packagePathForLintr), getwd(), packagePathForLintr)
message(">>>>>>>> packagePathForLintr: ", packagePathForLintr)

# look for forbidden lines
preventPattern <- function(patternToPrevent, packagePathForLintr) {
  message(sprintf("Checking pattern \"%s\" for package \"%s\"...", patternToPrevent, packagePathForLintr))
  
  out <- suppressWarnings({
    system(sprintf("cd %s; grep -R --include='*.R' '%s'", 
                   packagePathForLintr, patternToPrevent), 
           intern = TRUE)
  })
  
  tryCatch({
    stopifnot(length(out) == 0)
  }, error = function(err) {
    message(sprintf("Found lines: %s", out))
    stop(sprintf("'%s' code has been found in code!", patternToPrevent))
  })
}
preventPattern("browser()", packagePathForLintr)

tryCatch({
  # load pkg first to avoid non decalred variables/function
  suppressMessages(devtools::load_all(packagePathForLintr));
}, error = function(e) {
  message("It was unable to load package from directory: ", packagePathForLintr,
          "\nERROR: ", e)
  quit(status = 1)
})


library(lintr)

files_to_lint <- c(
  dir(paste0(packagePathForLintr, "/R"), full.names = TRUE),
  dir(paste0(packagePathForLintr, "/tests"), full.names = TRUE, recursive = TRUE),
  dir(paste0(packagePathForLintr, "/inst/shiny"), full.names = TRUE, recursive = FALSE, pattern = "*.R")
)

# ignore extensions
ignore.ext <- c("md", "html")
if (!is.null(ignore.ext) && any(vapply(ignore.ext, function(x) nchar(x) > 0, FUN.VALUE = logical(1)))) {
  pattern <- paste0(sprintf(".*\\.%s$", ignore.ext), collapse = "|")
  files_to_lint <- files_to_lint[!grepl(pattern, files_to_lint, ignore.case = TRUE)]
}

linters_config <- lintr::with_defaults(
  #object_usage_linter = NULL,
  #absolute_paths_linter = NULL,
  #assignment_linter = NULL,
  #closed_curly_linter = NULL,
  #commas_linter = NULL,
  cyclocomp_linter = NULL,
  #infix_spaces_linter = NULL,
  line_length_linter = lintr::line_length_linter(120),
  #no_tab_linter = NULL,
  camel_case_linter = NULL,
  snake_case_linter = NULL,
  object_name_linter = NULL,
  multiple_dots_linter = NULL,
  #multiple_dots_linter = NULL,
  #object_length_linter = NULL,
  #open_curly_linter = NULL,
  seq_linter = NULL,
  #single_quotes_linter = NULL,
  #spaces_inside_linter = NULL,
  #spaces_left_parentheses_linter = NULL,
  trailing_blank_lines_linter = NULL,
  trailing_whitespace_linter = NULL,
  #assignment_linter = NULL,
  object_usage_linter = NULL,
  object_length_linter = NULL
)

lint_file <- function(filepath) {
  print(paste("Linting file:", filepath))
  result <- lintr::lint(filepath, linters = linters_config)
  if (length(result) > 0) {
    print(result) # Show linter messages
    stop(paste0("Linter fails on file:", filepath))
  }
}

lint_all_files <- function(files_to_lint) {  
  for (filepath in files_to_lint) {
    lint_file(filepath)
  }
  print("All files OK!")
}

tryCatch({
  lint_all_files(files_to_lint)
}, error = function(e) {
  message(e)
  quit(status = 1)
})