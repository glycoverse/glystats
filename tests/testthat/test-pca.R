test_that("gly_pca works", {
  # Note: this integration test only makes sure the function runs,
  # it doesn't promise the result is correct.
  pca_res <- gly_pca(test_gp_exp)
  expect_s3_class(pca_res, c("glystats_pca_res", "glystats_res"))
  expect_type(pca_res, "list")
  expect_setequal(names(pca_res), c("tidy_result", "raw_result"))

  # Check tidy_result structure
  expect_s3_class(pca_res$tidy_result, c("glystats_pca_res", "glystats_res"))
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
  expect_s3_class(result$tidy_result, c("glystats_pca_res", "glystats_res"))
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