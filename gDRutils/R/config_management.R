REPLACEMENTS <- new.env(parent = emptyenv())

#' Set replacements
#'
#' @param ... replacements as name = value parameters
#'
#' @export
set_replacements <- function(...) {
  replacements <- list(...)
  checkmate::assert_list(replacements, any.missing = FALSE)
  checkmate::assert_character(names(replacements), min.chars = 1)

  for (n in names(replacements)) REPLACEMENTS[[n]] <- replacements[[n]]
}

#' Get replacements
#'
#' @return named list of replacements
#'
#' @export
get_replacements <- function() {
  r <- lapply(names(REPLACEMENTS), function(n) REPLACEMENTS[[n]])
  names(r) <- names(REPLACEMENTS)
  r
}

#' Replace {x} with values defined in replacements
#'
#' @param s string that may contain replacements
#' @param replacements named list containing replacements e.g list(from = 'to something new')
#' @return replaced string
#'
#' @export
str_replace <- function(s, replacements) {
  checkmate::assert_string(s)
  checkmate::assert_list(replacements, any.missing = FALSE)

  if (length(replacements)) {
    checkmate::assert_character(names(replacements), min.chars = 1)

    for (name in names(replacements)) s <- gsub(paste0("{", name, "}"), replacements[[name]], s, fixed = TRUE)
  }
  s
}

#' Load configuration
#'
#' @param ... optional paths to configuration files and/or paths to directories that may contain these files
#' @param use_default load default.yaml config from this package
#' @param replacements named list containing replacements e.g list(from = 'to something new')
#' @param evaluate list with flags responsible for evaluations R, PKG, ENV strings
#' @param dump flag responsible for displaying the contents of the final configuration
#' @return Config "object" (R list)
#'
#' @examples
#' \dontrun{
#' --- example config files ---
#'
#' gDRutils/config/default.yaml
#'   default: test
#'
#' a.yaml
#'   x1: a
#'   x2:
#'     x21: a
#'     x22: a
#'   x3: ["a1" , "a2"]
#'   r1: "{from_str_1}"
#'   r2: "{from_str_2}"
#'
#' b.yaml
#'   x1: b
#'   x2:
#'     x21: b
#'     x23: b
#'   e1: "`PKG:gDRviz`/subdir"
#'   e2: "`ENV:X`"
#'   e2: `R:paste0("evaluated")`
#'
#' -- R code --
#'
#' load_configuration(
#'  "a.yaml",
#'  repacements = list(from_str_1 = "to str 1", from_str_3 = "doesn't exist")
#' )
#'   default: test
#'   x1: a
#'   x2:
#'     x21: a
#'     x22: a
#'   x3: ["a1" , "a2"]
#'   r1: "to_str_1"
#'   r2: "{from_str_2}"
#'
#' load_configuration(
#'   "a.yaml", "b.yaml", use_default = FALSE,
#'   evaluate = list(PKG = TRUE, ENV = FALSE, R = TRUE)
#' )
#'   x1: b
#'   x2:
#'     x21: b
#'     x22: a
#'     x23: b
#'   x3: ["a1" , "a2"]
#'   r1: "{from_str_1}"
#'   r2: "{from_str_2}"
#'   e1: "/usr/local/lib/R/library/gDRviz/subdir"
#'   e2: "`ENV:X`"
#'   e2: evaluated
#'
#' Sys.setenv(X = "X value")
#' set_replacements(from_str_2 = "to two")
#' load_configuration(
#'   "a.yaml", "b.yaml",
#'   evaluate = list(PKG = FALSE, ENV = TRUE, R = FALSE)
#' )
#'   default: test
#'   x1: b
#'   x2:
#'     x21: b
#'     x22: a
#'     x23: b
#'   x3: ["a1" , "a2"]
#'   r1: "{from_str_1}"
#'   r2: "to_two"
#'   e1: "`PKG:gDRviz`/subdir"
#'   e2: "X value"
#'   e2: `R:paste0("evaluated")`
#' }
#'
#' @export
load_configuration <- function(
  ...,
  use_default = TRUE,
  replacements = get_replacements(),
  evaluate = list(R = FALSE, PKG = TRUE, ENV = TRUE),
  dump = FALSE
) {
  checkmate::assert_flag(use_default)
  checkmate::assert_flag(dump)

  config <- list()

  config_sources <- unlist(list(...))
  files <- c(
    if (use_default) system.file("configs", "default.yaml", package = "gDRutils"),
    unlist(lapply(config_sources, function(c) {
      if (dir.exists(c)) {
        get_path_configs(c)
      } else if (file.exists(c)) {
        c
      } else {
        c_by_name <- system.file("configs", paste0(c, ".yaml"), package = "gDRutils")
        if (file.exists(c_by_name)) c_by_name
        else c
      }
    }))
  )

  cat("loading configs:", paste0(files, collapse = ", "))

  for (i in seq_len(length(files))) {
    file <- read_config(files[[i]])
    if (is.null(file)) next

    # string replacements {..}
    file <- str_replace(file, replacements)

    # string evaluation `..`
    file <- str_eval(file, evaluate, paste0("file '", files[i], "'"))

    # merge
    tryCatch(
      config <- suppressWarnings(merge_configs(config, yaml::yaml.load(file))),
      error = function(e) {
        stop("error loading file", paste0("'", files[i], "',"), "message:", e$message)
      }
    )
  }

  if (dump) cat(paste0("configuration:\n", dump_config(config)))

  config
}

#' Get config files from directory
#'
#' @param dir path to a directory which may contain configuration files
#' @return NULL or character vector with config files
#'
#' @export
get_path_configs <- function(dir) {
  checkmate::assert_string(dir)

  sort(list.files(dir, pattern = "\\.ya?ml$", full.names = TRUE))
}

#' Get config files from environment variable(s)
#'
#' @param prefix beginning of the name of the environment variables
#' @return NULL or character vector of config files
#'
#' @export
get_env_configs <- function(prefix = "GDR_CONFIG") {
  checkmate::assert_string(prefix)

  env <- Sys.getenv()
  as.character(env[sort(names(env)[grepl(paste0("^", prefix), names(env), perl = TRUE)])])
}

#' Get config element by name
#'
#' @param config config object (R list)
#' @param path string or character vector containing the path to the config item
#' @param default default value for missing item
#' @param stop_on_missing call "stop" if required item is missing
#' @return config item
#'
#' @export
get_config_item <- function(config, path, default = NULL, stop_on_missing = FALSE) {
  checkmate::assert_list(config)
  checkmate::assert_character(path, min.chars = 1, min.len = 1)
  checkmate::assert_flag(stop_on_missing)

  item <- config
  for (name in path) {
    if (!any(name %in% names(item))) {
      stopifnot("missing required config item" = !stop_on_missing)
      return(default)
    }
    item <- item[[name]]
  }
  item
}

#' A function that creates a configuration dump
#'
#' @param config config (R list)
#' @param tab tab (prefix) string
#' @param eol end of line string
#' @param level current config level (for internal use)
#' @return string "tree" representation of the config
#'
#' @export
dump_config <- function(config, tab = "  ", eol = "\n", level = 0) {
  checkmate::assert_list(config)
  checkmate::assert_string(tab)
  checkmate::assert_string(eol)
  checkmate::assert_int(level)

  config <- config[order(names(config))]
  prefix <- paste0(rep(tab, level), collapse = "")
  out <- c()
  for (name in names(config)) {
    value <- config[[name]]
    if (is.list(value)) {
      out <- c(
        out,
        paste0(prefix, name, ":"),
        dump_config(
          config = value,
          tab = tab,
          eol = eol,
          level = level + 1
        )
      )
    } else {
      out <- c(
        out,
        paste0(prefix, name, ": ", paste0(value, collapse = ", "))
      )
    }
  }
  paste0(out, collapse = eol)
}


#' Function that saves the configuration to the yaml file
#'
#' @param config configuration "object"
#' @param path path to the output yaml file
#' @return TRUE on success, FALSE on error
#'
#' @export
write_config <- function(config, path) {
  checkmate::assert_list(config)
  checkmate::assert_string(path)

  tryCatch(
    expr = {
      if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
      suppressWarnings(yaml::write_yaml(config, path))
      cat(sprintf("configuration saved to file '%s'", path))
      TRUE
    },
    error = function(e) {
      stop(sprintf(
        "unable to save configuration to file '%s', original error message: '%s'", path, e$message
      ))
      FALSE
    }
  )
}

#' Function that read the configuration from the yaml file
#'
#' @param path path to the yaml file
#' @param as_yaml encode yaml content flag
#' @return NULL if error occurred, string if as_yaml = FALSE, otherwise config object
#'
#' @export
read_config <- function(path, as_yaml = FALSE) {
  checkmate::assert_string(path)
  checkmate::assert_flag(as_yaml)

  if (as_yaml) {
    return(tryCatch(
      expr = suppressWarnings(yaml::read_yaml(file)),
      error = function(e) {
        stop(sprintf("cannot read yaml file '%s', error message: '%s'", path, e$message))
        NULL
      }
    ))
  }

  tryCatch(
    expr = suppressWarnings(readChar(path, file.info(path)$size)),
    error = function(e) {
      stop(sprintf("cannot open file '%s', error message: '%s'", path, e$message))
      NULL
    }
  )
}

#' Helper that creates path to user config file by user config dir and username
#'
#' @param username user name
#' @param users_settings_dir path to directory where user setting file is (or will be) stored
#' @return path to user config file
#'
#' @export
user_config_file <- function(username, users_settings_dir) {
  checkmate::assert_string(username, min.chars = 1)
  checkmate::assert_string(users_settings_dir)

  file.path(users_settings_dir, paste0(username, ".yaml"))
}


#' Evaluate `PKG/ENV/R : ...` expressions in string
#'
#' @param s string that may contain `PKG/ENV/R : ...` expressions
#' @param evaluate list with flags responsible for evaluations R, PKG, ENV strings
#' @param s_name friendly name of s string that will be used in log messages
#' @return evaluated string
#'
#' @export
str_eval <- function(s, evaluate, s_name = paste0("'", s, "'")) {
  checkmate::assert_string(s)
  checkmate::assert_list(evaluate, types = "logical")
  checkmate::assert_names(names(evaluate), permutation.of = c("R", "PKG", "ENV"))
  checkmate::assert_string(s_name)

  if (!nzchar(s) || !length(evaluate)) return(s)

  offset <- 1
  while ((reg_res <- regexpr("`(PKG|R|ENV):[^`]+`", substr(s, offset, nchar(s)))) > 0) {
    expr_str <- substr(
      s,
      offset + reg_res,
      offset + reg_res + attr(reg_res, "match.length") - 3
    )
    expr_parts <- unlist(strsplit(expr_str, ":", fixed = TRUE))
    expr_type <- expr_parts[1]
    expr_val <- paste(expr_parts[-1], collapse = ":")

    if (evaluate[[expr_type]]) {
      res <- tryCatch(
        expr = {
          switch(expr_type,
            PKG = system.file(package = expr_val, mustWork = TRUE),
            ENV = Sys.getenv(expr_val),
            R = {
              res <- as.character(eval(parse(text = expr_val)))
              if (!length(res)) res <- ""
              res
            }
          )
        },
        error = function(e) {
          stop(sprintf(
            "error during evaluation `%s` in %s",
            expr_str,
            s_name
          ))
        }
      )
    } else {
      cat(sprintf(
        "'%s' evaluation disabled, `%s` in %s will not be proceed",
        expr_type,
        expr_str,
        s_name
      ))
      res <- NULL
    }

    if (is.null(res)) {
      offset <- offset + reg_res - 1 + attr(reg_res, "match.length")
    } else {
      s <- gsub(paste0("`", expr_str, "`"), res, s, fixed = TRUE)
      offset <- offset + reg_res - 1 + nchar(res)
    }
  }

  s
}


# FUNCTION COPIED FROM config R PACAKGE ON 17.06.2021
# Commit id: 2a5bb6ae4168be11abd296c875008aea8e8d52d7
# https://github.com/rstudio/config/blob/main/R/merge.R
# Remove it when https://github.com/rstudio/config/issues/37 is fixed
#
#' Merge two "config" lists
#'
#' @param base_list base list
#' @param overlay_list overlay list
#' @param recursive merge recursive flag
#' @return merged config list
#'
#' @export
merge_lists <- function(base_list, overlay_list, recursive = TRUE) {
  if (length(base_list) == 0) {
    overlay_list
  } else if (length(overlay_list) == 0) {
    base_list
  } else {
    merged_list <- base_list
    for (name in names(overlay_list)) {
      base <- base_list[[name]]
      overlay <- overlay_list[[name]]
      if (is.list(base) && is.list(overlay) && recursive) {
        merged_list[[name]] <- merge_lists(base, overlay)
      } else {
        merged_list[[name]] <- NULL
        merged_list <- append(merged_list,
                              overlay_list[which(names(overlay_list) %in% name)])
      }
    }
    merged_list
  }
}

#' Merge two "config" lists
#'
#' @param base_config base config list
#' @param overlay_config overlay config list
#' @param recursive merge recursive flag
#' @return merged config lists
#'
#' @export
merge_configs <- function(base_config, overlay_config, recursive = TRUE) {
  # internal helper unctions
  remove_exclamation_mark <- function(str) {
    if (substring(str, nchar(str)) == "!") substr(str, 1, nchar(str) - 1)
    else str
  }
  remove_exclamation_marks <- function(l) {
    l_names <- c()
    for (name in names(l)) {
      l_names <- c(l_names, remove_exclamation_mark(name))
      if (is.list(l[[name]])) l[[name]] <- remove_exclamation_marks(l[[name]])
    }
    names(l) <- l_names
    l
  }

  # nothing to merge
  if (length(base_config) == 0) return(remove_exclamation_marks(overlay_config))
  if (length(overlay_config) == 0) return(remove_exclamation_marks(base_config))

  # merge configs
  merged_config <- base_config <- remove_exclamation_marks(base_config)
  for (overlay_name in names(overlay_config)) {
    name <- remove_exclamation_mark(overlay_name)

    base <- base_config[[name]]
    overlay <- overlay_config[[overlay_name]]

    merged_config[[name]] <-
      if (is.list(base) && is.list(overlay) && recursive && name == overlay_name) {
        merge_configs(base, overlay)
      } else {
        remove_exclamation_marks(overlay)
      }
  }
  merged_config
}
