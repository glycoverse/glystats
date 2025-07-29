test_that("gly_plsda works with valid Topliss ratio", {
  # Skip test if mixOmics is not available
  skip_if_not_installed("mixOmics")

  # Test with a dataset that satisfies Topliss ratio (40 samples, 6 variables, ratio = 6.67)
  suppressMessages({
    plsda_res <- gly_plsda(exp_multigroup_valid())
  })

  expect_s3_class(plsda_res, c("glystats_plsda_res", "glystats_res"))
  expect_type(plsda_res, "list")
  expect_setequal(names(plsda_res), c("tidy_result", "raw_result"))

  # Check tidy_result structure
  expect_type(plsda_res$tidy_result, "list")
  expect_setequal(names(plsda_res$tidy_result), c("samples", "variables", "variance", "vip"))

  # Check samples tibble
  expect_s3_class(plsda_res$tidy_result$samples, "tbl_df")
  expect_true("group" %in% colnames(plsda_res$tidy_result$samples))
  expect_true("comp1" %in% colnames(plsda_res$tidy_result$samples))
  expect_true("comp2" %in% colnames(plsda_res$tidy_result$samples))

  # Check variables tibble
  expect_s3_class(plsda_res$tidy_result$variables, "tbl_df")
  expect_true("variable" %in% colnames(plsda_res$tidy_result$variables))
  expect_true("comp1" %in% colnames(plsda_res$tidy_result$variables))
  expect_true("comp2" %in% colnames(plsda_res$tidy_result$variables))

  # Check variance tibble
  expect_s3_class(plsda_res$tidy_result$variance, "tbl_df")
  expect_true("component" %in% colnames(plsda_res$tidy_result$variance))
  expect_true("prop_var_explained" %in% colnames(plsda_res$tidy_result$variance))
  expect_true("cumulative_prop_var" %in% colnames(plsda_res$tidy_result$variance))

  # Check VIP tibble
  expect_s3_class(plsda_res$tidy_result$vip, "tbl_df")
  expect_true("variable" %in% colnames(plsda_res$tidy_result$vip))
  expect_true("VIP" %in% colnames(plsda_res$tidy_result$vip))
  expect_true(all(plsda_res$tidy_result$vip$VIP >= 0))  # VIP scores should be non-negative

  # Check raw_result
  expect_s3_class(plsda_res$raw_result, "mixo_plsda")
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

test_that("gly_plsda raw_result is accessible", {
  # Skip test if mixOmics is not available
  skip_if_not_installed("mixOmics")

  # Test raw_result access with valid dataset
  suppressMessages({
    plsda_res <- gly_plsda(exp_multigroup_valid())
  })

  expect_s3_class(plsda_res$raw_result, "mixo_plsda")
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

test_that("gly_plsda_ works correctly", {
  skip_if_not_installed("mixOmics")

  # Create test data with good Topliss ratio
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(60)) + 1, nrow = 6, ncol = 10)  # 6 vars, 10 samples
  rownames(expr_mat) <- paste0("var", 1:6)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B"), each = 5))

  # Test function execution
  suppressMessages({
    result <- gly_plsda_(expr_mat, groups)
  })

  # Verify results structure
  expect_s3_class(result, c("glystats_plsda_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))

  # Check tidy_result
  expect_type(result$tidy_result, "list")
  expect_setequal(names(result$tidy_result), c("samples", "variables", "variance", "vip"))
  expect_true(tibble::is_tibble(result$tidy_result$samples))
  expect_equal(nrow(result$tidy_result$samples), 10)

  # Check raw_result
  expect_s3_class(result$raw_result, "mixo_plsda")
})

test_that("gly_plsda add_info works", {
  skip_if_not_installed("mixOmics")

  # Test with add_info = TRUE (default)
  suppressMessages({
    plsda_with_info <- gly_plsda(exp_multigroup_valid(), add_info = TRUE)
  })

  # Test with add_info = FALSE
  suppressMessages({
    plsda_without_info <- gly_plsda(exp_multigroup_valid(), add_info = FALSE)
  })

  # Results should be different when add_info is different
  expect_false(identical(
    plsda_with_info$tidy_result$samples,
    plsda_without_info$tidy_result$samples
  ))
})
