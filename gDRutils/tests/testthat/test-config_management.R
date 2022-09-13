testthat::context("Test config management")

test_that("set_replacements & get_replacements", {
  set_replacements(from1 = "to1", from2 = "to2")
  expect_equal(list(from1 = "to1", from2 = "to2"), get_replacements())
})

test_that("str_replace", {
  expect_equal(
    str_replace("a{b}c{d}e{f}", list(b = 1, f = 2)),
    "a1c{d}e2"
  )
})

test_that("str_eval", {
  s <- "`PKG:X`, `PKG:base`, `ENV:X`, `R:paste0('evaluated')`"

  expect_equal(
    str_eval(s, list(PKG = TRUE, ENV = FALSE, R = TRUE)),
    "`PKG:X`, /usr/local/lib/R/library/base, `ENV:X`, evaluated"
  )

  expect_equal(
    withr::with_envvar(
      c(X = "Y"),
      str_eval(s, list(PKG = FALSE, ENV = TRUE, R = FALSE))
    ),
    "`PKG:X`, `PKG:base`, Y, `R:paste0('evaluated')`"
  )
})

test_that("load_configuration", {
  set_replacements(r1 = "one")
  CONFIG <- load_configuration("testthat/a", "testthat/b", use_default = FALSE)
  expect_equal(
    CONFIG,
    list(
      x1 = "b",
      x2 = list(x21 = "a", x22 = "b", x23 = "b"),
      o = list(o2 = 2, o3 = 2),
      r1 = "one",
      r2 = "{r2}",
      x3 = "b",
      e = "/usr/local/lib/R/library/base"
    )
  )
})


test_that("load_configuration() prints a helpful error message", {
  a <- capture_messages(
    load_configuration("testthat/broken", use_default = FALSE)
  )

  expect_true(
    any(grepl(
      "error loading.*broken.yaml.*at line 2",
      a
    ))
  )
})


test_that("get_path_configs", {
  expect_equal(
    basename(get_path_configs(system.file("configs", "testthat", package = "gDRutils"))),
    c("a.yaml", "b.yaml", "broken.yaml")
  )
})

test_that("get_env_configs", {
  expect_equal(
    withr::with_envvar(
      c(GDR_CONFIG_TT_1 = "config/zzz", GDR_CONFIG_TT_2 = "config/aaa"),
      get_env_configs("GDR_CONFIG_TT")
    ),
    c("config/zzz", "config/aaa")
  )
})

test_that("get_config_item", {
  CONFIG <- list(
    x1 = "one",
    x2 = list(x21 = "a", x22 = "b", x23 = "c")
  )
  expect_equal(get_config_item(CONFIG, "x1"), "one")
  expect_equal(get_config_item(CONFIG, "undefined"), NULL)
  expect_equal(get_config_item(CONFIG, c("x2", "x22")), "b")
  expect_error(get_config_item(CONFIG, "undefined", stop_on_missing = TRUE))
})
