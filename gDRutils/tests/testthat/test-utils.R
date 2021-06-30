test_that(".clean_key_inputs works as expected", {
  keys <- LETTERS[1:5]
  cols <- LETTERS[1:3]
  expect_warning(out <- gDRutils:::.clean_key_inputs(keys, cols))
  expect_equal(out, cols)
})

