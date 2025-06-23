test_that("gly_plsda works with basic functionality", {
  # Skip if mixOmics is not available
  skip_if_not_installed("mixOmics")
  
  # Note: this integration test only makes sure the function runs,
  # it doesn't promise the result is correct.
  
  # Use test_gp_exp with 2 groups for PLS-DA (filter to manageable size)
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing
  
  result <- gly_plsda(exp_2group)
  
  # Test basic structure
  expect_s3_class(result, c("glystats_plsda_res", "glystats_res"))
  expect_type(result, "list")
  expect_named(result, c("samples", "variables", "components"))
  
  # Check samples tibble
  expect_s3_class(result$samples, "tbl_df")
  expect_true("sample" %in% names(result$samples))
  expect_true("comp1" %in% names(result$samples))
  expect_true("comp2" %in% names(result$samples))
  
  # Check variables tibble
  expect_s3_class(result$variables, "tbl_df")
  expect_true("variable" %in% names(result$variables))
  expect_true("comp1" %in% names(result$variables))
  expect_true("comp2" %in% names(result$variables))
  
  # Check components tibble
  expect_s3_class(result$components, "tbl_df")
  expect_true("component" %in% names(result$components))
  expect_true("explained_variance" %in% names(result$components))
  expect_true("cumulative_variance" %in% names(result$components))
})

test_that("gly_plsda works with different number of components", {
  skip_if_not_installed("mixOmics")
  
  # Use test_gp_exp with small subset
  exp_small <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test with 3 components
  result <- gly_plsda(exp_small, ncomp = 3)
  
  expect_true("comp3" %in% names(result$samples))
  expect_true("comp3" %in% names(result$variables))
  expect_equal(nrow(result$components), 3)
})

test_that("gly_plsda works with custom group column", {
  skip_if_not_installed("mixOmics")
  
  # Create experiment with different group column name
  exp_custom <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::rename_obs(treatment = group)
  
  # Test with custom group column
  result <- gly_plsda(exp_custom, group_col = "treatment")
  
  expect_s3_class(result, c("glystats_plsda_res", "glystats_res"))
  expect_true("sample" %in% names(result$samples))
})

test_that("gly_plsda handles errors correctly", {
  skip_if_not_installed("mixOmics")
  
  exp_small <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::select_obs(-group)  # Remove group column
  
  # Should error when default "group" column doesn't exist
  expect_error(
    gly_plsda(exp_small),
    "Group variable group not found"
  )
  
  # Should error when specified group variable doesn't exist
  expect_error(
    gly_plsda(test_gp_exp, group_col = "nonexistent"),
    "Group variable nonexistent not found"
  )
}) 