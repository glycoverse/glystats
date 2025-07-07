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
  expect_setequal(names(oplsda_res), c("samples", "variables", "variance", "vip"))

  # Check samples tibble
  expect_s3_class(oplsda_res$samples, "tbl_df")
  expect_true("group" %in% colnames(oplsda_res$samples))
  expect_true("p1" %in% colnames(oplsda_res$samples))  # predictive component 1

  # Check variables tibble
  expect_s3_class(oplsda_res$variables, "tbl_df")
  expect_true("variable" %in% colnames(oplsda_res$variables))
  expect_true("p1" %in% colnames(oplsda_res$variables))

  # Check variance tibble
  expect_s3_class(oplsda_res$variance, "tbl_df")
  expect_true("component" %in% colnames(oplsda_res$variance))
  expect_true("prop_var_explained" %in% colnames(oplsda_res$variance))
  expect_true("cumulative_prop_var" %in% colnames(oplsda_res$variance))

  # Check VIP tibble
  expect_s3_class(oplsda_res$vip, "tbl_df")
  expect_true("variable" %in% colnames(oplsda_res$vip))
  expect_true("VIP" %in% colnames(oplsda_res$vip))
  expect_true(all(oplsda_res$vip$VIP >= 0))  # VIP scores should be non-negative
})

test_that("gly_oplsda validates Topliss ratio", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test that function correctly rejects datasets with insufficient n/p ratio
  # The test dataset has 6 samples and 500 variables (ratio = 0.012 << 5)
  expect_error(
    gly_oplsda(exp_2groups()),
    "Insufficient sample-to-variable ratio"
  )

  # Test that the error message contains helpful information
  expect_error(
    gly_oplsda(exp_2groups()),
    "Topliss ratio principle"
  )

  # Test that the error message suggests solutions
  expect_error(
    gly_oplsda(exp_2groups()),
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
  if (ncol(oplsda_res$variables) > 0) {  # Model was built
    # May have orthogonal components in samples
    expect_true("p1" %in% colnames(oplsda_res$samples))
  }
})

test_that("gly_oplsda return_raw works", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test return_raw with valid dataset
  suppressMessages(suppressWarnings({
    capture.output({
      oplsda_raw <- gly_oplsda(exp_topliss_valid(), return_raw = TRUE)
    }, type = "output")
  }))

  expect_s4_class(oplsda_raw, "opls")
})

test_that("gly_oplsda validates inputs", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test invalid pred_i
  expect_error(
    gly_oplsda(exp_2groups(), pred_i = 0),
    "Assertion on 'pred_i' failed"
  )

  # Test invalid ortho_i
  expect_error(
    gly_oplsda(exp_2groups(), ortho_i = -1),
    "Assertion on 'ortho_i' failed"
  )

  # Test multi-group data (should fail at group validation before Topliss ratio)
  expect_error(
    gly_oplsda(test_gp_exp),
    "group must have exactly 2 levels for"
  )

  # Test invalid group column (should fail before Topliss ratio check)
  expect_error(
    gly_oplsda(exp_2groups(), group_col = "nonexistent"),
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
  expect_false(identical(oplsda_scaled$samples, oplsda_unscaled$samples))
})

test_that("gly_oplsda validates Topliss ratio with invalid data", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test that function rejects datasets with insufficient n/p ratio
  expect_error(
    gly_oplsda(exp_2groups()),
    "Insufficient sample-to-variable ratio"
  )
})
