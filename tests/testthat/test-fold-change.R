test_that("gly_fold_change works with basic 2-group comparison", {
  # Use test_gp_exp and filter to 2 groups
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing

  # Run fold change calculation with add_info = FALSE for basic test
  result <- suppressMessages(gly_fold_change(exp_2group, add_info = FALSE))

  # Test core functionality
  expect_s3_class(result, "glystats_fc_res")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 2)
  expect_setequal(colnames(result), c("variable", "log2fc"))
  expect_type(result$variable, "character")
  expect_type(result$log2fc, "double")

  # Check that all variables are included
  expect_setequal(result$variable, glyexp::get_var_info(exp_2group)$variable)

  # Test with add_info = TRUE (default)
  result_with_info <- suppressMessages(gly_fold_change(exp_2group))
  expect_s3_class(result_with_info, "glystats_fc_res")
  expect_equal(nrow(result_with_info), 10)
  expect_true(ncol(result_with_info) > 2)  # Should have more columns with var_info
  expect_true("variable" %in% colnames(result_with_info))
  expect_true("log2fc" %in% colnames(result_with_info))
})

test_that("gly_fold_change works with custom group column", {
  # Create test data with custom group column
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::mutate_obs(treatment = group)  # Add a custom group column
  
  # Test with custom group column and add_info = FALSE
  result <- suppressMessages(gly_fold_change(exp_2group, group_col = "treatment", add_info = FALSE))

  expect_s3_class(result, "glystats_fc_res")
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 2)
  expect_setequal(colnames(result), c("variable", "log2fc"))
})

test_that("gly_fold_change handles factor groups correctly", {
  # Create test data with factor groups
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::mutate_obs(group = factor(group, levels = c("H", "C")))  # H as reference
  
  result <- suppressMessages(gly_fold_change(exp_2group))
  
  expect_s3_class(result, "glystats_fc_res")
  expect_equal(nrow(result), 5)
  expect_type(result$log2fc, "double")
  
  # Check that fold change is calculated correctly (C vs H, H is reference)
  # All values should be finite (no NA/Inf)
  expect_true(all(is.finite(result$log2fc)))
})

test_that("gly_fold_change error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)
  
  # Test with non-glyexp_experiment object
  expect_error(gly_fold_change("not_an_experiment"))
  
  # Test with non-existent group column
  expect_error(gly_fold_change(exp_small, group_col = "nonexistent"), 
               "not found in sample information")
  
  # Test with non-string group_col
  expect_error(gly_fold_change(exp_small, group_col = 123))
})

test_that("gly_fold_change group validation", {
  # Test with more than 2 groups
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)
  
  expect_error(suppressMessages(gly_fold_change(exp_3group)), "exactly 2 levels")
  
  # Test with only 1 group
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)
  
  expect_error(suppressMessages(gly_fold_change(exp_1group)), "exactly 2 levels")
})

test_that("gly_fold_change outputs informative messages", {
  # Test that function outputs group information
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 3)
  
  # Capture messages
  messages <- capture.output(
    result <- gly_fold_change(exp_2group),
    type = "message"
  )
  
  # Should contain group information
  expect_true(any(stringr::str_detect(messages, "Group 1")))
  expect_true(any(stringr::str_detect(messages, "Group 2")))
})

test_that("gly_fold_change handles edge cases", {
  # Test with minimal data (just 1 variable)
  exp_minimal <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 1)
  
  result <- suppressMessages(gly_fold_change(exp_minimal, add_info = FALSE))

  expect_s3_class(result, "glystats_fc_res")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  
  # Test with many variables to ensure performance
  exp_large <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H"))
  
  expect_no_error(suppressMessages(gly_fold_change(exp_large)))
}) 