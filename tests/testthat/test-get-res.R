test_that("get_tidy_result works", {
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_s3_class(get_tidy_result(result, "main_test"), "tbl_df")
})

test_that("get_raw_result works", {
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_type(get_raw_result(result), "list")
})

test_that("get_tidy_result raises error if which is not in tidy_result", {
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_error(get_tidy_result(result, "nonexistent"), "must be one of")
})

test_that("get_tidy_result raises error if tidy_result is a list but no `which` is provided", {
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_error(get_tidy_result(result), "must be provided for")
})