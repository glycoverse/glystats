test_that("gly_kmeans works with basic parameters", {
  # Use a subset of test data for faster testing
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10) |>
    glyexp::slice_sample_obs(n = 8)

  result <- suppressMessages(gly_kmeans(exp_subset))

  # Check that result is a tibble with expected structure
  expect_s3_class(result, "tbl_df")
  expect_s3_class(result, "glystats_kmeans_res")
  expect_s3_class(result, "glystats_res")

  # Check columns (should cluster variables by default)
  expect_true("variable" %in% colnames(result))
  expect_true("cluster" %in% colnames(result))

  # Check that we have the right number of rows
  expect_equal(nrow(result), 10)
})

test_that("gly_kmeans_ works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  
  # Test function execution
  suppressMessages({
    result <- gly_kmeans_(expr_mat)
  })
  
  # Verify results
  expect_s3_class(result, "glystats_kmeans_res")
  expect_true(tibble::is_tibble(result))
  expect_true("variable" %in% colnames(result))
  expect_true("cluster" %in% colnames(result))
  expect_equal(nrow(result), 10)
})
