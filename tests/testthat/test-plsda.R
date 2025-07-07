test_that("gly_plsda works with valid Topliss ratio", {
  # Skip test if mixOmics is not available
  skip_if_not_installed("mixOmics")

  # Test with a dataset that satisfies Topliss ratio (40 samples, 6 variables, ratio = 6.67)
  suppressMessages({
    plsda_res <- gly_plsda(exp_multigroup_valid())
  })

  expect_s3_class(plsda_res, c("glystats_plsda_res", "glystats_res"))
  expect_type(plsda_res, "list")
  expect_setequal(names(plsda_res), c("samples", "variables", "variance", "vip"))

  # Check samples tibble
  expect_s3_class(plsda_res$samples, "tbl_df")
  expect_true("group" %in% colnames(plsda_res$samples))
  expect_true("comp1" %in% colnames(plsda_res$samples))
  expect_true("comp2" %in% colnames(plsda_res$samples))

  # Check variables tibble
  expect_s3_class(plsda_res$variables, "tbl_df")
  expect_true("variable" %in% colnames(plsda_res$variables))
  expect_true("comp1" %in% colnames(plsda_res$variables))
  expect_true("comp2" %in% colnames(plsda_res$variables))

  # Check variance tibble
  expect_s3_class(plsda_res$variance, "tbl_df")
  expect_true("component" %in% colnames(plsda_res$variance))
  expect_true("prop_var_explained" %in% colnames(plsda_res$variance))
  expect_true("cumulative_prop_var" %in% colnames(plsda_res$variance))

  # Check VIP tibble
  expect_s3_class(plsda_res$vip, "tbl_df")
  expect_true("variable" %in% colnames(plsda_res$vip))
  expect_true("VIP" %in% colnames(plsda_res$vip))
  expect_true(all(plsda_res$vip$VIP >= 0))  # VIP scores should be non-negative
})

test_that("gly_plsda validates Topliss ratio", {
  # Skip test if mixOmics is not available
  skip_if_not_installed("mixOmics")

  # Test that function correctly rejects datasets with insufficient n/p ratio
  # The test dataset has 12 samples and 500 variables (ratio = 0.024 << 5)
  expect_error(
    suppressMessages(gly_plsda(test_gp_exp)),
    "Insufficient sample-to-variable ratio"
  )

  # Test that the error message contains helpful information
  expect_error(
    suppressMessages(gly_plsda(test_gp_exp)),
    "Topliss ratio principle"
  )

  # Test that the error message suggests solutions
  expect_error(
    suppressMessages(gly_plsda(test_gp_exp)),
    "Collecting more samples"
  )
})

test_that("gly_plsda return_raw works", {
  # Skip test if mixOmics is not available
  skip_if_not_installed("mixOmics")

  # Test return_raw with valid dataset
  suppressMessages({
    plsda_raw <- gly_plsda(exp_multigroup_valid(), return_raw = TRUE)
  })

  expect_s3_class(plsda_raw, "mixo_plsda")
})

test_that("gly_plsda validates inputs", {
  # Skip test if mixOmics is not available
  skip_if_not_installed("mixOmics")

  # Test invalid group column (should fail before Topliss ratio check)
  expect_error(
    suppressMessages(gly_plsda(test_gp_exp, group_col = "nonexistent")),
    "not found in sample information"
  )

  # Test invalid ncomp (should fail before Topliss ratio check)
  expect_error(
    suppressMessages(gly_plsda(test_gp_exp, ncomp = 0)),
    "Assertion on 'ncomp' failed"
  )
})
