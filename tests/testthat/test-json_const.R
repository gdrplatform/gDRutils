
test_that("get_settings_from_json works as expected", {
  expect_error(
    get_settings_from_json(json_path = "/no/such/path"),
    "Assertion on 'json_path' failed"
  )
  json_path <-
    system.file(package = "gDRutils", "test_settings.json")
  expect_error(
    get_settings_from_json(s = "no_such_entry", json_path = json_path),
    "Assertion on 's' failed"
  )
  
  s <- get_settings_from_json(json_path = json_path)
  expect_true(is.list(s))
  exp_names <-
    c("DICT",
      "DICT_WITH_LOGICAL",
      "DICT_WITH_LISTS",
      "LIST",
      "STRING")
  expect_identical(sort(names(s)), sort(exp_names))
  expect_true(is.logical(s$DICT_WITH_LOGICAL$axisTitleText))
  expect_identical(s$DICT_WITH_LISTS$`GR AOC within set range`, c(0L, 1L))
  
  s2 <- get_settings_from_json(s = "DICT", json_path = json_path)
  expect_identical(s$DICT, s2)
  
  expect_no_error(get_settings_from_json())
})



test_that("get_isobologram_columns as expected", {
  expect_equal(
    get_isobologram_columns(),
    c("Iso_Level", "Pos_x", "Pos_x_Ref", "Pos_y", "Pos_y_Ref", "Log10_Ratio_Conc", "Log2_CI"))
  
  expect_equal(
    get_isobologram_columns(prettify = FALSE),
    c("iso_level", "pos_x", "pos_x_ref", "pos_y", "pos_y_ref", "log10_ratio_conc", "log2_CI"))
  
  expect_equal(get_isobologram_columns("iso_level", prettify = TRUE), "Iso_Level")
  expect_equal(get_isobologram_columns("iso_level", prettify = FALSE), "iso_level")
  
  expect_error(get_isobologram_columns(1))
  expect_error(get_isobologram_columns(prettify = 2))
})
