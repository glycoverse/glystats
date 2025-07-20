test_that("gly_ttest works with t-test method", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing

  # Run DEA with t-test
  result <- suppressMessages(gly_ttest(exp_2group))

  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_ttest", "glystats_dea_res", "glystats_res"))
  expect_equal(nrow(result), 10)
  expect_true("p_adj" %in% colnames(result))  # p_adj should exist
  expect_true("log2fc" %in% colnames(result))  # log2fc should exist
  expect_type(result$log2fc, "double")  # log2fc should be numeric
})

test_that("gly_wilcox works with wilcoxon method", {
  # Use test_gp_exp and filter to 2 groups for wilcoxon test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing

  # Run DEA with wilcoxon test
  result <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group)))

  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_wilcoxon", "glystats_dea_res", "glystats_res"))
  expect_equal(nrow(result), 10)
  expect_true("log2fc" %in% colnames(result))  # Wilcoxon should now have log2fc
  expect_type(result$log2fc, "double")  # log2fc should be numeric
})

test_that("gly_ttest and gly_wilcox basic functionality works", {
  # Test both methods work with test_gp_exp
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)  # Use very small subset
  
  # 2-group methods
  exp_2group <- exp_small |> glyexp::filter_obs(group %in% c("C", "H"))
  expect_no_error(suppressMessages(gly_ttest(exp_2group)))
  expect_no_error(suppressMessages(suppressWarnings(gly_wilcox(exp_2group))))
})

test_that("gly_ttest and gly_wilcox error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)
  
  # Test various error conditions - group column not found
  expect_error(suppressMessages(gly_ttest(exp_small, group_col = "nonexistent")),
               "not found in sample information")
  expect_error(suppressMessages(gly_wilcox(exp_small, group_col = "nonexistent")),
               "not found in sample information")
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
  expect_equal(result_default$log2fc, -result_ref_h$log2fc, tolerance = 1e-10)
  expect_equal(result_default$log2fc, result_ref_c$log2fc, tolerance = 1e-10)

  # Test invalid ref_group
  expect_error(suppressMessages(gly_ttest(exp_2group, ref_group = "invalid")),
               "Must be element of set")
})

test_that("gly_wilcox ref_group parameter works", {
  # Use test_gp_exp and filter to 2 groups for wilcoxon test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with default (no ref_group)
  result_default <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group)))

  # Test with ref_group = "Y" (should reverse the comparison)
  result_ref_y <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group, ref_group = "Y")))

  # Test with ref_group = "M" (should be same as default since M is first alphabetically)
  result_ref_m <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group, ref_group = "M")))

  # Check that log2fc values are negated when reference group changes
  expect_equal(result_default$log2fc, -result_ref_y$log2fc, tolerance = 1e-10)
  expect_equal(result_default$log2fc, result_ref_m$log2fc, tolerance = 1e-10)

  # Test invalid ref_group
  expect_error(suppressMessages(gly_wilcox(exp_2group, ref_group = "invalid")),
               "Must be element of set")
})
