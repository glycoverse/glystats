test_that("gly_ttest works with t-test method", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10) # Use smaller subset for faster testing

  # Run DEA with t-test
  result <- suppressMessages(gly_ttest(exp_2group))

  # Test core functionality
  expect_s3_class(result, c("glystats_ttest_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))
  expect_equal(nrow(result$tidy_result), 10)
  expect_true("p_adj" %in% colnames(result$tidy_result)) # p_adj should exist
  expect_true("log2fc" %in% colnames(result$tidy_result)) # log2fc should exist
  expect_type(result$tidy_result$log2fc, "double") # log2fc should be numeric

  # Test raw_result
  expect_type(result$raw_result, "list")
  expect_equal(length(result$raw_result), 10)
  expect_true(all(purrr::map_lgl(result$raw_result, ~ inherits(.x, "htest"))))
})

test_that("gly_ttest assigns NA for failed variables", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)
  exp_2group$expr_mat[1:3, ] <- NA # This will lead to stats::t.test() failing
  na_vars <- exp_2group$var_info$variable[1:3]

  # Run DEA with t-test
  expect_warning(result <- suppressMessages(gly_ttest(exp_2group)))

  # Test results
  main_test_raw <- result$raw_result
  expect_true(all(is.na(main_test_raw[na_vars])))
  main_test_tidy <- result$tidy_result
  p_values <- main_test_tidy |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(p_val)
  expect_true(all(is.na(p_values)))
})

test_that("gly_wilcox works with wilcoxon method", {
  # Use test_gp_exp and filter to 2 groups for wilcoxon test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 10) # Use smaller subset for faster testing

  # Run DEA with wilcoxon test
  result <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group)))

  # Test core functionality
  expect_s3_class(result, c("glystats_wilcox_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))
  expect_equal(nrow(result$tidy_result), 10)
  expect_true("log2fc" %in% colnames(result$tidy_result)) # Wilcoxon should now have log2fc
  expect_type(result$tidy_result$log2fc, "double") # log2fc should be numeric

  # Test raw_result
  expect_type(result$raw_result, "list")
  expect_equal(length(result$raw_result), 10)
  expect_true(all(purrr::map_lgl(result$raw_result, ~ inherits(.x, "htest"))))
})

test_that("gly_ttest_ works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B"), each = 5))

  # Test function execution
  suppressMessages({
    result <- gly_ttest_(expr_mat, groups)
  })

  # Verify results
  expect_s3_class(result, "glystats_ttest_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_true(nrow(result$tidy_result) > 0)
  expect_true("log2fc" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)
})

test_that("gly_wilcox_ works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B"), each = 5))

  # Test function execution
  suppressMessages({
    result <- gly_wilcox_(expr_mat, groups)
  })

  # Verify results
  expect_s3_class(result, "glystats_wilcox_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_true(nrow(result$tidy_result) > 0)
  expect_true("log2fc" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)
})

test_that("gly_ttest and gly_wilcox basic functionality works", {
  # Test both methods work with test_gp_exp
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5) # Use very small subset

  # 2-group methods
  exp_2group <- exp_small |> glyexp::filter_obs(group %in% c("C", "H"))
  expect_no_error(suppressMessages(gly_ttest(exp_2group)))
  expect_no_error(suppressMessages(suppressWarnings(gly_wilcox(exp_2group))))
})

test_that("gly_ttest and gly_wilcox error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)

  # Test various error conditions - group column not found
  expect_error(
    suppressMessages(gly_ttest(exp_small, group_col = "nonexistent")),
    "not found in sample information"
  )
  expect_error(
    suppressMessages(gly_wilcox(exp_small, group_col = "nonexistent")),
    "not found in sample information"
  )
})

test_that("gly_ttest and gly_wilcox group validation", {
  # Test with 3 groups (using C, H, M from test_gp_exp)
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with 1 group
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)

  # Test 3 groups with 2-group methods
  expect_error(suppressMessages(gly_ttest(exp_3group)), "exactly 2 levels")
  expect_error(suppressMessages(gly_wilcox(exp_3group)), "exactly 2 levels")

  # Test 1 group with 2-group methods
  expect_error(suppressMessages(gly_ttest(exp_1group)), "exactly 2 levels")
  expect_error(suppressMessages(gly_wilcox(exp_1group)), "exactly 2 levels")
})

test_that("gly_ttest ref_group parameter works", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with default (no ref_group)
  result_default <- suppressMessages(gly_ttest(exp_2group))

  # Test with ref_group = "H" (should reverse the comparison)
  result_ref_h <- suppressMessages(gly_ttest(exp_2group, ref_group = "H"))

  # Test with ref_group = "C" (should be same as default since C is first alphabetically)
  result_ref_c <- suppressMessages(gly_ttest(exp_2group, ref_group = "C"))

  # Check that log2fc values are negated when reference group changes
  expect_equal(
    result_default$tidy_result$log2fc,
    -result_ref_h$tidy_result$log2fc,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$log2fc,
    result_ref_c$tidy_result$log2fc,
    tolerance = 1e-10
  )

  # Test invalid ref_group
  expect_error(
    suppressMessages(gly_ttest(exp_2group, ref_group = "invalid")),
    "Must be element of set"
  )
})

test_that("gly_wilcox ref_group parameter works", {
  # Use test_gp_exp and filter to 2 groups for wilcoxon test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with default (no ref_group)
  result_default <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group)))

  # Test with ref_group = "Y" (should reverse the comparison)
  result_ref_y <- suppressMessages(suppressWarnings(gly_wilcox(
    exp_2group,
    ref_group = "Y"
  )))

  # Test with ref_group = "M" (should be same as default since M is first alphabetically)
  result_ref_m <- suppressMessages(suppressWarnings(gly_wilcox(
    exp_2group,
    ref_group = "M"
  )))

  # Check that log2fc values are negated when reference group changes
  expect_equal(
    result_default$tidy_result$log2fc,
    -result_ref_y$tidy_result$log2fc,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$log2fc,
    result_ref_m$tidy_result$log2fc,
    tolerance = 1e-10
  )

  # Test invalid ref_group
  expect_error(
    suppressMessages(gly_wilcox(exp_2group, ref_group = "invalid")),
    "Must be element of set"
  )
})
