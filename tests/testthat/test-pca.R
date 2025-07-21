test_that("gly_pca works", {
  # Note: this integration test only makes sure the function runs,
  # it doesn't promise the result is correct.
  pca_res <- gly_pca(test_gp_exp)
  expect_s3_class(pca_res, c("glystats_pca_res", "glystats_res"))
  expect_type(pca_res, "list")
  expect_setequal(names(pca_res), c("samples", "variables", "eigenvalues"))
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

  # Verify results
  expect_s3_class(result, "glystats_pca_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("samples", "variables", "eigenvalues"))
  expect_true(tibble::is_tibble(result$samples))
  expect_true(nrow(result$samples) > 0)  # Should have some samples
})