test_that("gly_kmeans works with basic parameters", {
  # Use a subset of test data for faster testing
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10) |>
    glyexp::slice_sample_obs(n = 8) |>
    as_test_se()

  result <- suppressMessages(gly_kmeans(exp_subset))

  # Check that result is a list with expected structure
  expect_type(result, "list")
  expect_s3_class(result, "glystats_kmeans_res")
  expect_s3_class(result, "glystats_res")

  # Check list elements
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check tidy_result structure
  expect_s3_class(result$tidy_result, "tbl_df")
  expect_true("variable" %in% colnames(result$tidy_result))
  expect_true("cluster" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)

  # Check raw_result structure
  expect_s3_class(result$raw_result, "kmeans")
})

test_that(".analyze_kmeans works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test function execution
  suppressMessages({
    result <- .analyze_kmeans(expr_mat)
  })

  # Verify results
  expect_type(result, "list")
  expect_s3_class(result, "glystats_kmeans_res")
  expect_s3_class(result, "glystats_res")

  # Check list elements
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check tidy_result structure
  expect_s3_class(result$tidy_result, "tbl_df")
  expect_true("variable" %in% colnames(result$tidy_result))
  expect_true("cluster" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)

  # Check raw_result structure
  expect_s3_class(result$raw_result, "kmeans")
})

test_that("gly_kmeans works with on='sample'", {
  # Use a subset of test data for faster testing
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10) |>
    glyexp::slice_sample_obs(n = 8)

  result <- suppressMessages(gly_kmeans(exp_subset, on = "sample"))

  # Check that result is a list with expected structure
  expect_type(result, "list")
  expect_s3_class(result, "glystats_kmeans_res")
  expect_s3_class(result, "glystats_res")

  # Check tidy_result structure (should cluster samples)
  expect_s3_class(result$tidy_result, "tbl_df")
  expect_true("sample" %in% colnames(result$tidy_result))
  expect_true("cluster" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 8)

  # Check raw_result structure
  expect_s3_class(result$raw_result, "kmeans")
})

test_that("gly_kmeans add_info parameter works", {
  # Use a subset of test data for faster testing
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10) |>
    glyexp::slice_sample_obs(n = 8)

  # Test with add_info = TRUE (default)
  result_with_info <- suppressMessages(gly_kmeans(exp_subset, add_info = TRUE))

  # Test with add_info = FALSE
  result_without_info <- suppressMessages(gly_kmeans(
    exp_subset,
    add_info = FALSE
  ))

  # Both should have the same structure but different column counts
  expect_type(result_with_info, "list")
  expect_type(result_without_info, "list")

  # The result with add_info should have more columns
  expect_gte(
    ncol(result_with_info$tidy_result),
    ncol(result_without_info$tidy_result)
  )
})
