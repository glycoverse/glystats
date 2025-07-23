test_that("gly_oplsda works with valid Topliss ratio", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test with a dataset that satisfies Topliss ratio (30 samples, 5 variables, ratio = 6)
  suppressMessages(suppressWarnings({
    capture.output({
      oplsda_res <- gly_oplsda(exp_topliss_valid())
    }, type = "output")
  }))

  expect_s3_class(oplsda_res, c("glystats_oplsda_res", "glystats_res"))
  expect_type(oplsda_res, "list")
  expect_setequal(names(oplsda_res), c("tidy_result", "raw_result"))

  # Check tidy_result structure
  expect_s3_class(oplsda_res$tidy_result, c("glystats_oplsda_res", "glystats_res"))
  expect_type(oplsda_res$tidy_result, "list")
  expect_setequal(names(oplsda_res$tidy_result), c("samples", "variables", "variance", "vip"))

  # Check samples tibble
  expect_s3_class(oplsda_res$tidy_result$samples, "tbl_df")
  expect_true("group" %in% colnames(oplsda_res$tidy_result$samples))
  expect_true("p1" %in% colnames(oplsda_res$tidy_result$samples))  # predictive component 1

  # Check variables tibble
  expect_s3_class(oplsda_res$tidy_result$variables, "tbl_df")
  expect_true("variable" %in% colnames(oplsda_res$tidy_result$variables))
  expect_true("p1" %in% colnames(oplsda_res$tidy_result$variables))

  # Check variance tibble
  expect_s3_class(oplsda_res$tidy_result$variance, "tbl_df")
  expect_true("component" %in% colnames(oplsda_res$tidy_result$variance))
  expect_true("prop_var_explained" %in% colnames(oplsda_res$tidy_result$variance))
  expect_true("cumulative_prop_var" %in% colnames(oplsda_res$tidy_result$variance))

  # Check VIP tibble
  expect_s3_class(oplsda_res$tidy_result$vip, "tbl_df")
  expect_true("variable" %in% colnames(oplsda_res$tidy_result$vip))
  expect_true("VIP" %in% colnames(oplsda_res$tidy_result$vip))
  expect_true(all(oplsda_res$tidy_result$vip$VIP >= 0))  # VIP scores should be non-negative

  # Check raw_result
  expect_s4_class(oplsda_res$raw_result, "opls")
})

test_that("gly_oplsda validates Topliss ratio", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test that function correctly rejects datasets with insufficient n/p ratio
  # The test dataset has 6 samples and 500 variables (ratio = 0.012 << 5)
  expect_error(
    suppressMessages(gly_oplsda(exp_2groups())),
    "Insufficient sample-to-variable ratio"
  )

  # Test that the error message contains helpful information
  expect_error(
    suppressMessages(gly_oplsda(exp_2groups())),
    "Topliss ratio principle"
  )

  # Test that the error message suggests solutions
  expect_error(
    suppressMessages(gly_oplsda(exp_2groups())),
    "Collecting more samples"
  )
})

test_that("gly_oplsda works with orthogonal components", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test with valid dataset and orthogonal components
  suppressMessages(suppressWarnings({
    capture.output({
      oplsda_res <- gly_oplsda(exp_topliss_valid(), ortho_i = 1)
    }, type = "output")
  }))

  expect_s3_class(oplsda_res, c("glystats_oplsda_res", "glystats_res"))

  # Check that orthogonal components are present if model was built successfully
  if (ncol(oplsda_res$tidy_result$variables) > 0) {  # Model was built
    # May have orthogonal components in samples
    expect_true("p1" %in% colnames(oplsda_res$tidy_result$samples))
  }
})

test_that("gly_oplsda raw_result is accessible", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test raw_result access with valid dataset
  suppressMessages(suppressWarnings({
    capture.output({
      oplsda_res <- gly_oplsda(exp_topliss_valid())
    }, type = "output")
  }))

  expect_s4_class(oplsda_res$raw_result, "opls")
})

test_that("gly_oplsda validates inputs", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test invalid pred_i
  expect_error(
    suppressMessages(gly_oplsda(exp_2groups(), pred_i = 0)),
    "Assertion on 'pred_i' failed"
  )

  # Test invalid ortho_i
  expect_error(
    suppressMessages(gly_oplsda(exp_2groups(), ortho_i = -1)),
    "Assertion on 'ortho_i' failed"
  )

  # Test multi-group data (should fail at group validation before Topliss ratio)
  expect_error(
    suppressMessages(gly_oplsda(test_gp_exp)),
    "group must have exactly 2 levels for"
  )

  # Test invalid group column (should fail before Topliss ratio check)
  expect_error(
    suppressMessages(gly_oplsda(exp_2groups(), group_col = "nonexistent")),
    "not found in sample information"
  )
})

test_that("gly_oplsda handles different scaling options", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test with scaling
  suppressMessages(suppressWarnings({
    capture.output({
      oplsda_scaled <- gly_oplsda(exp_topliss_valid(), scale = TRUE)
    }, type = "output")
  }))
  expect_s3_class(oplsda_scaled, c("glystats_oplsda_res", "glystats_res"))

  # Test without scaling
  suppressMessages(suppressWarnings({
    capture.output({
      oplsda_unscaled <- gly_oplsda(exp_topliss_valid(), scale = FALSE)
    }, type = "output")
  }))
  expect_s3_class(oplsda_unscaled, c("glystats_oplsda_res", "glystats_res"))

  # Results should be different when scaling is different
  expect_false(identical(oplsda_scaled$tidy_result$samples, oplsda_unscaled$tidy_result$samples))
})

test_that("gly_oplsda validates Topliss ratio with invalid data", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test that function rejects datasets with insufficient n/p ratio
  expect_error(
    suppressMessages(gly_oplsda(exp_2groups())),
    "Insufficient sample-to-variable ratio"
  )
})

test_that("gly_oplsda_ works correctly", {
  skip_if_not_installed("ropls")

  exp <- exp_topliss_valid()
  expr_mat <- glyexp::get_expr_mat(exp)
  groups <- factor(glyexp::get_sample_info(exp)$group)

  # Test function execution
  suppressMessages(suppressWarnings({
    capture.output({
      result <- gly_oplsda_(expr_mat, groups)
    }, type = "output")
  }))

  # Verify results structure
  expect_s3_class(result, c("glystats_oplsda_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_s3_class(result$tidy_result, c("glystats_oplsda_res", "glystats_res"))
  expect_s4_class(result$raw_result, "opls")
})
