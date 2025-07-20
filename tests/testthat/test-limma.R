test_that("gly_limma works with limma method", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Use test_gp_exp and filter to 2 groups for limma
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing

  # Run DEA with limma
  result <- suppressMessages(gly_limma(exp_2group))

  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_limma", "glystats_dea_res", "glystats_res"))
  expect_equal(nrow(result), 10)
  expect_true("p_adj" %in% colnames(result))  # p_adj should exist
  expect_true("log2fc" %in% colnames(result))  # log2fc should exist
  expect_type(result$log2fc, "double")  # log2fc should be numeric
  # coefficient列已被重命名为log2fc，所以不需要单独的coefficient列
  expect_true("t" %in% colnames(result))  # t-statistic
  expect_true("b" %in% colnames(result))  # log-odds (B-statistic)
})

test_that("gly_limma basic functionality works", {
  # Only test limma if available
  skip_if_not_installed("limma")
  
  # Test limma works with test_gp_exp
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)  # Use very small subset
  
  # 2-group method
  exp_2group <- exp_small |> glyexp::filter_obs(group %in% c("C", "H"))
  expect_no_error(suppressMessages(gly_limma(exp_2group)))
})

test_that("gly_limma error handling", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")
  
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)
  
  # Test various error conditions - group column not found
  expect_error(suppressMessages(gly_limma(exp_small, group_col = "nonexistent")),
               "not found in sample information")
})

test_that("gly_limma group validation", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")
  
  # Test with 3 groups (using C, H, M from test_gp_exp)
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test with 1 group 
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)
  
  # Test 3 groups with 2-group method
  expect_error(suppressMessages(gly_limma(exp_3group)), "exactly 2 levels")
  
  # Test 1 group with 2-group method
  expect_error(suppressMessages(gly_limma(exp_1group)), "exactly 2 levels")
})

test_that("gly_limma ref_group parameter works", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Use test_gp_exp and filter to 2 groups for limma
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with default (no ref_group)
  result_default <- suppressMessages(gly_limma(exp_2group))

  # Test with ref_group = "H" (should reverse the comparison)
  result_ref_h <- suppressMessages(gly_limma(exp_2group, ref_group = "H"))

  # Test with ref_group = "C" (should be same as default since C is first alphabetically)
  result_ref_c <- suppressMessages(gly_limma(exp_2group, ref_group = "C"))

  # Check that log2fc values are negated when reference group changes
  expect_equal(result_default$log2fc, -result_ref_h$log2fc, tolerance = 1e-10)
  expect_equal(result_default$log2fc, result_ref_c$log2fc, tolerance = 1e-10)

  # Test invalid ref_group
  expect_error(suppressMessages(gly_limma(exp_2group, ref_group = "invalid")),
               "Must be element of set")
})

test_that("gly_limma works with real data", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")
  
  # This test uses a subset of test_gp_exp to ensure integration works
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 20)  # Use reasonable subset for testing
    
  result <- suppressMessages(gly_limma(exp_2group))
  expect_s3_class(result, c("glystats_dea_res_limma", "glystats_dea_res", "glystats_res"))
  expect_true(tibble::is_tibble(result))
  expect_true("log2fc" %in% colnames(result))
  expect_true("p_adj" %in% colnames(result))
})
