test_that("gly_oplsda works with valid Topliss ratio", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test with a dataset that satisfies Topliss ratio (30 samples, 5 variables, ratio = 6)
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_res <- gly_oplsda(exp_topliss_valid())
      },
      type = "output"
    )
  }))

  expect_s3_class(oplsda_res, c("glystats_oplsda_res", "glystats_res"))
  expect_type(oplsda_res, "list")
  expect_setequal(
    names(oplsda_res),
    c("tidy_result", "raw_result", "meta_data")
  )

  # Check tidy_result structure
  expect_type(oplsda_res$tidy_result, "list")
  expect_setequal(
    names(oplsda_res$tidy_result),
    c("samples", "variables", "variance", "vip", "perm_test")
  )

  # Check samples tibble
  expect_s3_class(oplsda_res$tidy_result$samples, "tbl_df")
  expect_true("group" %in% colnames(oplsda_res$tidy_result$samples))
  expect_true("p1" %in% colnames(oplsda_res$tidy_result$samples)) # predictive component 1

  # Check variables tibble
  expect_s3_class(oplsda_res$tidy_result$variables, "tbl_df")
  expect_true("variable" %in% colnames(oplsda_res$tidy_result$variables))
  expect_true("p1" %in% colnames(oplsda_res$tidy_result$variables))
  expect_true("pcorr1" %in% colnames(oplsda_res$tidy_result$variables))

  # Check variance tibble
  expect_s3_class(oplsda_res$tidy_result$variance, "tbl_df")
  expect_true("component" %in% colnames(oplsda_res$tidy_result$variance))
  expect_true(
    "prop_var_explained" %in% colnames(oplsda_res$tidy_result$variance)
  )
  expect_true(
    "cumulative_prop_var" %in% colnames(oplsda_res$tidy_result$variance)
  )

  # Check VIP tibble
  expect_s3_class(oplsda_res$tidy_result$vip, "tbl_df")
  expect_true("variable" %in% colnames(oplsda_res$tidy_result$vip))
  expect_true("vip" %in% colnames(oplsda_res$tidy_result$vip))
  expect_true(all(oplsda_res$tidy_result$vip$vip >= 0)) # VIP scores should be non-negative

  # Check perm_test tibble
  expect_s3_class(oplsda_res$tidy_result$perm_test, "tbl_df")
  expect_true("model" %in% colnames(oplsda_res$tidy_result$perm_test))
  expect_true("perm_id" %in% colnames(oplsda_res$tidy_result$perm_test))

  # If permutation test was performed, check structure
  if (nrow(oplsda_res$tidy_result$perm_test) > 0) {
    expect_true("Original" %in% oplsda_res$tidy_result$perm_test$model)
    expect_true(0 %in% oplsda_res$tidy_result$perm_test$perm_id)
    expect_true(all(oplsda_res$tidy_result$perm_test$perm_id >= 0))
  }

  # Check raw_result
  expect_s4_class(oplsda_res$raw_result, "opls")
})

test_that("gly_oplsda works with orthogonal components", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test with valid dataset and orthogonal components
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_res <- gly_oplsda(exp_topliss_valid(), ortho_i = 1)
      },
      type = "output"
    )
  }))

  expect_s3_class(oplsda_res, c("glystats_oplsda_res", "glystats_res"))

  # Check that orthogonal components are present if model was built successfully
  if (ncol(oplsda_res$tidy_result$variables) > 0) {
    # Model was built
    # May have orthogonal components in samples
    expect_true("p1" %in% colnames(oplsda_res$tidy_result$samples))
  }
})

test_that("gly_oplsda correctly handles orthogonal loadings and scores", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Create a small dataset that allows orthogonal components
  # Use 2 variables with 12 samples to ensure we can force orthogonal components
  suppressMessages({
    small_exp <- test_gp_exp |>
      glyexp::slice_head_var(n = 2) |>
      glyexp::mutate_obs(
        group = dplyr::if_else(group %in% c("H", "M"), "control", "case")
      )
  })

  # Test without forced orthogonal components (automatic decision)
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_auto <- gly_oplsda(small_exp)
      },
      type = "output"
    )
  }))

  # Test with forced orthogonal components
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_forced <- gly_oplsda(small_exp, ortho_i = 1)
      },
      type = "output"
    )
  }))

  # Verify automatic case (should not have orthogonal components)
  expect_true("p1" %in% colnames(oplsda_auto$tidy_result$samples))
  expect_false("o1" %in% colnames(oplsda_auto$tidy_result$samples))
  expect_true("p1" %in% colnames(oplsda_auto$tidy_result$variables))
  expect_false("o1" %in% colnames(oplsda_auto$tidy_result$variables))

  # Verify forced case (should have orthogonal components)
  expect_true("p1" %in% colnames(oplsda_forced$tidy_result$samples))
  expect_true("o1" %in% colnames(oplsda_forced$tidy_result$samples))
  expect_true("p1" %in% colnames(oplsda_forced$tidy_result$variables))
  expect_true("o1" %in% colnames(oplsda_forced$tidy_result$variables))

  # Verify that orthogonal scores and loadings have correct dimensions
  expect_equal(nrow(oplsda_forced$tidy_result$samples), 12) # 12 samples
  expect_equal(nrow(oplsda_forced$tidy_result$variables), 2) # 2 variables

  # Verify that orthogonal values are numeric and not all NA
  expect_type(oplsda_forced$tidy_result$samples$o1, "double")
  expect_type(oplsda_forced$tidy_result$variables$o1, "double")
  expect_false(all(is.na(oplsda_forced$tidy_result$samples$o1)))
  expect_false(all(is.na(oplsda_forced$tidy_result$variables$o1)))

  # Verify that raw ropls object contains orthogonal components
  expect_equal(ncol(oplsda_forced$raw_result@orthoScoreMN), 1)
  expect_equal(ncol(oplsda_forced$raw_result@orthoLoadingMN), 1)
  expect_equal(nrow(oplsda_forced$raw_result@orthoScoreMN), 12)
  expect_equal(nrow(oplsda_forced$raw_result@orthoLoadingMN), 2)
})

test_that("gly_oplsda orthogonal components regression test", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # This test specifically prevents regression of the bug where orthogonal loadings
  # were not included in the variables tibble
  suppressMessages({
    small_exp <- test_gp_exp |>
      glyexp::slice_head_var(n = 2) |>
      glyexp::mutate_obs(
        group = dplyr::if_else(group %in% c("H", "M"), "control", "case")
      )
  })

  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_res <- gly_oplsda(small_exp, ortho_i = 1)
      },
      type = "output"
    )
  }))

  # Critical regression test: ensure both orthogonal scores AND loadings are present
  samples_cols <- colnames(oplsda_res$tidy_result$samples)
  variables_cols <- colnames(oplsda_res$tidy_result$variables)

  # Both samples and variables should have predictive component
  expect_true("p1" %in% samples_cols)
  expect_true("p1" %in% variables_cols)

  # Both samples and variables should have orthogonal component
  expect_true(
    "o1" %in% samples_cols,
    info = "Orthogonal scores (o1) missing from samples tibble"
  )
  expect_true(
    "o1" %in% variables_cols,
    info = "Orthogonal loadings (o1) missing from variables tibble - this was the bug!"
  )

  # Verify the values are meaningful (not all zeros or NAs)
  expect_false(all(oplsda_res$tidy_result$samples$o1 == 0))
  expect_false(all(oplsda_res$tidy_result$variables$o1 == 0))
  expect_false(all(is.na(oplsda_res$tidy_result$samples$o1)))
  expect_false(all(is.na(oplsda_res$tidy_result$variables$o1)))
})

test_that("gly_oplsda raw_result is accessible", {
  # Skip test if ropls is not available
  skip_if_not_installed("ropls")

  # Test raw_result access with valid dataset
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_res <- gly_oplsda(exp_topliss_valid())
      },
      type = "output"
    )
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
    capture.output(
      {
        oplsda_scaled <- gly_oplsda(exp_topliss_valid(), scale = TRUE)
      },
      type = "output"
    )
  }))
  expect_s3_class(oplsda_scaled, c("glystats_oplsda_res", "glystats_res"))

  # Test without scaling
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_unscaled <- gly_oplsda(exp_topliss_valid(), scale = FALSE)
      },
      type = "output"
    )
  }))
  expect_s3_class(oplsda_unscaled, c("glystats_oplsda_res", "glystats_res"))

  # Results should be different when scaling is different
  expect_false(identical(
    oplsda_scaled$tidy_result$samples,
    oplsda_unscaled$tidy_result$samples
  ))
})

test_that("gly_oplsda_ works correctly", {
  skip_if_not_installed("ropls")

  exp <- exp_topliss_valid()
  expr_mat <- glyexp::get_expr_mat(exp)
  groups <- factor(glyexp::get_sample_info(exp)$group)

  # Test function execution
  suppressMessages(suppressWarnings({
    capture.output(
      {
        result <- gly_oplsda_(expr_mat, groups)
      },
      type = "output"
    )
  }))

  # Verify results structure
  expect_s3_class(result, c("glystats_oplsda_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_type(result$tidy_result, "list")
  expect_s4_class(result$raw_result, "opls")
})

test_that("gly_oplsda includes pcorr columns for predictive components", {
  skip_if_not_installed("ropls")

  # Test with single predictive component
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_res1 <- gly_oplsda(exp_topliss_valid(), pred_i = 1)
      },
      type = "output"
    )
  }))

  # Check that pcorr1 is present
  expect_true("pcorr1" %in% colnames(oplsda_res1$tidy_result$variables))

  # Check that pcorr1 values are numeric and within valid correlation range
  pcorr1_values <- oplsda_res1$tidy_result$variables$pcorr1
  expect_type(pcorr1_values, "double")
  expect_true(all(pcorr1_values >= -1 & pcorr1_values <= 1, na.rm = TRUE))
  expect_false(all(is.na(pcorr1_values)))

  # Test with multiple predictive components
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_res2 <- gly_oplsda(exp_topliss_valid(), pred_i = 2)
      },
      type = "output"
    )
  }))

  # Check that both pcorr1 and pcorr2 are present
  variables_cols <- colnames(oplsda_res2$tidy_result$variables)
  expect_true("pcorr1" %in% variables_cols)
  expect_true("pcorr2" %in% variables_cols)

  # Check that both pcorr columns have valid values
  pcorr1_values2 <- oplsda_res2$tidy_result$variables$pcorr1
  pcorr2_values2 <- oplsda_res2$tidy_result$variables$pcorr2

  expect_type(pcorr1_values2, "double")
  expect_type(pcorr2_values2, "double")
  expect_true(all(pcorr1_values2 >= -1 & pcorr1_values2 <= 1, na.rm = TRUE))
  expect_true(all(pcorr2_values2 >= -1 & pcorr2_values2 <= 1, na.rm = TRUE))
  expect_false(all(is.na(pcorr1_values2)))
  expect_false(all(is.na(pcorr2_values2)))

  # Verify that pcorr values are different between components
  expect_false(identical(pcorr1_values2, pcorr2_values2))
})

test_that("gly_oplsda pcorr values match manual calculation", {
  skip_if_not_installed("ropls")

  # Use a simple test case
  suppressMessages(suppressWarnings({
    capture.output(
      {
        oplsda_res <- gly_oplsda(exp_topliss_valid(), pred_i = 1)
      },
      type = "output"
    )
  }))

  # Extract raw OPLS-DA object for manual calculation
  opls_obj <- oplsda_res$raw_result

  # Get the first predictive component scores
  t1 <- opls_obj@scoreMN[, "p1", drop = TRUE]

  # Get the modeling matrix X
  X <- opls_obj@suppLs$xModelMN

  # Ensure sample alignment
  if (!is.null(rownames(X)) && !is.null(rownames(opls_obj@scoreMN))) {
    X <- X[rownames(opls_obj@scoreMN), , drop = FALSE]
  }

  # Calculate manual pcorr1
  manual_pcorr1 <- apply(X, 2, function(x) {
    cor(x, t1, use = "pairwise.complete.obs")
  })

  # Get pcorr1 from our function
  function_pcorr1 <- oplsda_res$tidy_result$variables$pcorr1
  names(function_pcorr1) <- oplsda_res$tidy_result$variables$variable

  # Align by variable names
  common_vars <- intersect(names(manual_pcorr1), names(function_pcorr1))
  manual_aligned <- manual_pcorr1[common_vars]
  function_aligned <- function_pcorr1[common_vars]

  # Compare values (allowing for small numerical differences)
  expect_equal(manual_aligned, function_aligned, tolerance = 1e-10)
})
