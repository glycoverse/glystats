test_that("gly_oplsda works with basic functionality", {
  # Skip if ropls is not available
  skip_if_not_installed("ropls")
  
  # Note: this integration test only makes sure the function runs,
  # it doesn't promise the result is correct.
  
  # Use test_gp_exp with 2 groups for OPLS-DA (filter to manageable size)
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing
  
  result <- gly_oplsda(exp_2group)
  
  # Test basic structure
  expect_s3_class(result, "gly_oplsda")
  expect_type(result, "list")
  expect_named(result, c("samples", "variables", "components"))
  
  # Check samples tibble - basic required columns
  expect_s3_class(result$samples, "tbl_df")
  expect_true("sample" %in% names(result$samples))
  expect_true("group" %in% names(result$samples))
  
  # Check variables tibble - basic structure
  expect_s3_class(result$variables, "tbl_df")
  
  # Check components tibble - basic structure
  expect_s3_class(result$components, "tbl_df")
  
  # If components exist, they should have the right structure
  if (nrow(result$components) > 0) {
    expect_true("component" %in% names(result$components))
    expect_true("type" %in% names(result$components))
    expect_true("explained_variance" %in% names(result$components))
    expect_true("cumulative_variance" %in% names(result$components))
  }
})

test_that("gly_oplsda works with different number of components", {
  skip_if_not_installed("ropls")
  
  # Use test_gp_exp with small subset
  exp_small <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test with 2 predictive components and 1 orthogonal component
  # Note: ropls may automatically adjust the number of components based on significance
  # Suppress the expected warning about component number adjustment
  result <- suppressWarnings(gly_oplsda(exp_small, predI = 2, orthoI = 1))
  
  # Basic structure checks
  expect_s3_class(result, "gly_oplsda")
  expect_type(result, "list")
  expect_named(result, c("samples", "variables", "components"))
  
  # Check that we have the required columns even if no components are significant
  expect_true("sample" %in% names(result$samples))
  expect_true("group" %in% names(result$samples))
  
  # Check variables tibble structure
  expect_s3_class(result$variables, "tbl_df")
  expect_s3_class(result$components, "tbl_df")
})

test_that("gly_oplsda works with custom group column", {
  skip_if_not_installed("ropls")
  
  # Create experiment with different group column name
  exp_custom <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::rename_obs(treatment = group)
  
  # Test with custom group column
  result <- gly_oplsda(exp_custom, group_col = "treatment")
  
  expect_true("treatment" %in% names(result$samples))
})

test_that("gly_oplsda works with different scaling options", {
  skip_if_not_installed("ropls")
  
  # Use test_gp_exp with small subset
  exp_small <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test with no centering or scaling
  result_none <- gly_oplsda(exp_small, center = FALSE, scale = FALSE)
  expect_s3_class(result_none, "gly_oplsda")
  
  # Test with only centering
  result_center <- gly_oplsda(exp_small, center = TRUE, scale = FALSE)
  expect_s3_class(result_center, "gly_oplsda")
  
  # Test with only scaling
  result_scale <- gly_oplsda(exp_small, center = FALSE, scale = TRUE)
  expect_s3_class(result_scale, "gly_oplsda")
})

test_that("gly_oplsda handles errors correctly", {
  skip_if_not_installed("ropls")
  
  exp_small <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::select_obs(-group)  # Remove group column
  
  # Should error when default "group" column doesn't exist
  expect_error(
    gly_oplsda(exp_small),
    "Group variable group not found"
  )
  
  # Should error when specified group variable doesn't exist
  expect_error(
    gly_oplsda(test_gp_exp, group_col = "nonexistent"),
    "Group variable nonexistent not found"
  )
}) 