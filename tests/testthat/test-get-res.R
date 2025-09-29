test_that("get_tidy_result works", {
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_s3_class(get_tidy_result(result, "main_test"), "tbl_df")
  expect_type(get_tidy_result(result), "list")
})

test_that("get_raw_result works", {
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_type(get_raw_result(result), "list")
})

test_that("get_tidy_result raises error if which is not in tidy_result", {
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_error(get_tidy_result(result, "nonexistent"), "must be one of")
})