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

  # Test with 1 group
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)

  # Test 1 group should still fail (need at least 2 groups)
  expect_error(suppressMessages(gly_limma(exp_1group)), "at least 2 levels")
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

test_that("gly_limma works with multi-group data", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Test with 3 groups (using C, H, M from test_gp_exp)
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing

  result <- suppressMessages(gly_limma(exp_3group))

  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_limma", "glystats_dea_res", "glystats_res"))
  expect_true(tibble::is_tibble(result))
  expect_true("contrast" %in% colnames(result))  # Should have contrast column
  expect_true("log2fc" %in% colnames(result))
  expect_true("p_adj" %in% colnames(result))

  # Should have 3 pairwise comparisons: C_vs_H, C_vs_M, H_vs_M
  expect_equal(length(unique(result$contrast)), 3)
  expect_true(all(c("C_vs_H", "C_vs_M", "H_vs_M") %in% result$contrast))

  # Each contrast should have the same number of variables
  expect_equal(nrow(result), 10 * 3)  # 10 variables * 3 contrasts
})

test_that("gly_limma multi-group generates correct contrasts", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Test with 4 groups to verify all pairwise comparisons (using C, H, M, Y)
  exp_4group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M", "Y")) |>
    glyexp::slice_sample_var(n = 5)  # Use very small subset for speed

  result <- suppressMessages(gly_limma(exp_4group))

  # Should have 6 pairwise comparisons for 4 groups: C(4,2) = 6
  expect_equal(length(unique(result$contrast)), 6)

  # Check all expected contrasts are present
  expected_contrasts <- c("C_vs_H", "C_vs_M", "C_vs_Y", "H_vs_M", "H_vs_Y", "M_vs_Y")
  expect_setequal(unique(result$contrast), expected_contrasts)
})

test_that("gly_limma custom contrasts work with hyphen format", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Test with 4 groups using custom contrasts (hyphen format)
  exp_4group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M", "Y")) |>
    glyexp::slice_sample_var(n = 5)

  # Test custom contrasts with hyphen format
  custom_contrasts <- c("H-C", "H-M", "H-Y")
  result <- suppressMessages(gly_limma(exp_4group, contrasts = custom_contrasts))

  # Should have exactly 3 contrasts as specified
  expect_equal(length(unique(result$contrast)), 3)
  expect_setequal(unique(result$contrast), c("H_vs_C", "H_vs_M", "H_vs_Y"))
  expect_equal(nrow(result), 5 * 3)  # 5 variables * 3 contrasts
})

test_that("gly_limma custom contrasts work with _vs_ format", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Test with 4 groups using custom contrasts (_vs_ format)
  exp_4group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M", "Y")) |>
    glyexp::slice_sample_var(n = 5)

  # Test custom contrasts with _vs_ format
  custom_contrasts <- c("H_vs_C", "M_vs_C", "Y_vs_C")
  result <- suppressMessages(gly_limma(exp_4group, contrasts = custom_contrasts))

  # Should have exactly 3 contrasts as specified
  expect_equal(length(unique(result$contrast)), 3)
  expect_setequal(unique(result$contrast), c("H_vs_C", "M_vs_C", "Y_vs_C"))
  expect_equal(nrow(result), 5 * 3)  # 5 variables * 3 contrasts
})

test_that("gly_limma contrasts error handling works", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Test with 3 groups
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)

  # Test invalid group names
  expect_error(
    suppressMessages(gly_limma(exp_3group, contrasts = c("H-C", "X-Y"))),
    "Group.*not found in data"
  )

  # Test invalid contrast format
  expect_error(
    suppressMessages(gly_limma(exp_3group, contrasts = c("H_C"))),
    "Invalid contrast format"
  )
})

test_that("gly_limma handles group names with hyphens correctly", {
  # Skip test if limma is not available
  skip_if_not_installed("limma")

  # Create test data with group names containing hyphens
  exp_hyphen <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)

  # Modify group names to include hyphens
  sample_info_modified <- glyexp::get_sample_info(exp_hyphen)
  sample_info_modified$group <- factor(
    ifelse(sample_info_modified$group == "C", "Control-1",
           ifelse(sample_info_modified$group == "H", "High-dose", "Medium-dose")),
    levels = c("Control-1", "High-dose", "Medium-dose")
  )

  # Create new experiment with modified sample info
  exp_hyphen_modified <- glyexp::experiment(
    expr_mat = glyexp::get_expr_mat(exp_hyphen),
    sample_info = sample_info_modified,
    var_info = glyexp::get_var_info(exp_hyphen),
    exp_type = "glycoproteomics",
    glycan_type = "N"
  )

  # Test that hyphen format fails with helpful error message
  expect_error(
    suppressMessages(gly_limma(exp_hyphen_modified, contrasts = c("High-dose-Control-1"))),
    "Group names contain hyphens.*Use the format.*_vs_"
  )

  # Test that _vs_ format works
  result <- suppressMessages(gly_limma(exp_hyphen_modified, contrasts = c("High-dose_vs_Control-1")))
  expect_equal(length(unique(result$contrast)), 1)
  expect_equal(unique(result$contrast), "High-dose_vs_Control-1")
})
