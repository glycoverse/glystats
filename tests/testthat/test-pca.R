test_that("gly_pca works", {
  # Note: this integration test only makes sure the function runs,
  # it doesn't promise the result is correct.
  pca_res <- gly_pca(test_gp_exp)
  expect_s3_class(pca_res, c("glystats_pca_res", "glystats_res"))
  expect_type(pca_res, "list")
  expect_setequal(names(pca_res), c("tidy_result", "raw_result"))

  # Check tidy_result structure
  expect_type(pca_res$tidy_result, "list")
  expect_setequal(names(pca_res$tidy_result), c("samples", "variables", "eigenvalues"))

  # Check raw_result
  expect_s3_class(pca_res$raw_result, "prcomp")
})

test_that("gly_pca_ works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test function execution
  suppressMessages({
    result <- gly_pca_(expr_mat)
  })

  # Verify results structure
  expect_s3_class(result, c("glystats_pca_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))

  # Check tidy_result
  expect_type(result$tidy_result, "list")
  expect_setequal(names(result$tidy_result), c("samples", "variables", "eigenvalues"))
  expect_true(tibble::is_tibble(result$tidy_result$samples))
  expect_true(nrow(result$tidy_result$samples) > 0)  # Should have some samples

  # Check raw_result
  expect_s3_class(result$raw_result, "prcomp")
})

test_that("gly_pca add_info works", {
  # Test with add_info = TRUE (default)
  pca_with_info <- gly_pca(test_gp_exp, add_info = TRUE)

  # Test with add_info = FALSE
  pca_without_info <- gly_pca(test_gp_exp, add_info = FALSE)

  # Results should be different when add_info is different
  expect_false(identical(
    pca_with_info$tidy_result$samples,
    pca_without_info$tidy_result$samples
  ))
})

test_that("gly_pca_ handles constant columns with scale = TRUE", {
  # Create test data with a constant column
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Make the first variable constant (all same value)
  expr_mat[1, ] <- 1

  # With scale = TRUE (default), should warn and remove the constant column
  expect_warning(
    result <- gly_pca_(expr_mat, scale = TRUE),
    "Removed 1 constant column before PCA"
  )
  expect_s3_class(result, "glystats_pca_res")
  # The constant variable should not be in the loadings
  expect_false("var1" %in% result$tidy_result$variables$variable)

  # With scale = FALSE, no need to remove constant columns
  expect_no_warning(
    result_no_scale <- gly_pca_(expr_mat, scale = FALSE)
  )
  expect_s3_class(result_no_scale, "glystats_pca_res")
  # The constant variable should be in the loadings (with zero loading)
  expect_true("var1" %in% result_no_scale$tidy_result$variables$variable)
})

test_that("gly_pca_ handles all constant columns", {
  # Create test data where all columns become constant after log transformation
  expr_mat <- matrix(0, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Should error because no columns with variance remain
  expect_error(
    expect_warning(gly_pca_(expr_mat, scale = TRUE)),
    "No columns with non-zero variance remain"
  )
})