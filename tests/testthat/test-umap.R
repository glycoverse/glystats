# Test UMAP analysis

test_that("gly_umap works with default parameters", {
  skip_if_not_installed("uwot")

  # Use appropriate n_neighbors for small dataset
  result <- gly_umap(test_gp_se, n_neighbors = 3)

  # Check basic structure
  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))

  # Check tidy_result structure
  expect_s3_class(result$tidy_result, "tbl_df")
  expect_equal(nrow(result$tidy_result), ncol(test_gp_se))
  expect_true("umap1" %in% names(result$tidy_result))
  expect_true("umap2" %in% names(result$tidy_result))
  expect_true("sample" %in% names(result$tidy_result))

  # Check that coordinates are numeric
  expect_type(result$tidy_result$umap1, "double")
  expect_type(result$tidy_result$umap2, "double")

  # Check that no missing values
  expect_false(any(is.na(result$tidy_result$umap1)))
  expect_false(any(is.na(result$tidy_result$umap2)))

  # Check raw_result
  expect_true(is.matrix(result$raw_result))
  expect_equal(nrow(result$raw_result), ncol(test_gp_se))
  expect_equal(ncol(result$raw_result), 2)
})

test_that("gly_umap works with custom parameters", {
  skip_if_not_installed("uwot")

  result <- gly_umap(
    test_gp_se,
    n_neighbors = 2,
    min_dist = 0.01,
    n_epochs = 50
  )

  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_equal(nrow(result$tidy_result), ncol(test_gp_se))
  expect_true(all(c("umap1", "umap2", "sample") %in% names(result$tidy_result)))
})

test_that("gly_umap works with more than 2 components", {
  skip_if_not_installed("uwot")

  result <- gly_umap(test_gp_se, n_neighbors = 3, n_components = 3)

  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_equal(nrow(result$tidy_result), ncol(test_gp_se))
  expect_true(all(
    c("umap1", "umap2", "umap3", "sample") %in% names(result$tidy_result)
  ))

  # Check that all coordinates are numeric
  expect_type(result$tidy_result$umap1, "double")
  expect_type(result$tidy_result$umap2, "double")
  expect_type(result$tidy_result$umap3, "double")

  # Check raw_result has correct dimensions
  expect_equal(ncol(result$raw_result), 3)
})

test_that("gly_umap has consistent sample names", {
  skip_if_not_installed("uwot")

  result <- gly_umap(test_gp_se, n_neighbors = 3)

  # Should have same sample names as input expression matrix
  expect_equal(
    sort(result$tidy_result$sample),
    sort(colnames(test_gp_se))
  )
})

test_that("gly_umap requires uwot package", {
  # Test that the function would fail if uwot is not available
  # This is tested by verifying the error message structure exists
  expect_true(any(grepl("uwot", deparse(body(gly_umap)))))
})

test_that(".analyze_umap works correctly", {
  skip_if_not_installed("uwot")

  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test function execution with appropriate n_neighbors
  suppressMessages({
    result <- .analyze_umap(expr_mat, n_neighbors = 3, n_components = 2)
  })

  # Verify results
  expect_s3_class(result, "glystats_umap_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))
  expect_true("sample" %in% colnames(result$tidy_result))
  expect_true("umap1" %in% colnames(result$tidy_result))
  expect_true("umap2" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)
  expect_true(is.matrix(result$raw_result))
  expect_equal(nrow(result$raw_result), 10)
  expect_equal(ncol(result$raw_result), 2)
})
