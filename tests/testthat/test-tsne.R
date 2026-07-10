# Test t-SNE analysis

test_that("gly_tsne works with default parameters", {
  skip_if_not_installed("Rtsne")

  # Use appropriate perplexity for small dataset to avoid warnings
  result <- gly_tsne(test_gp_se, perplexity = 3)

  # Check basic structure
  expect_s3_class(result, c("glystats_tsne_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))

  # Check tidy_result structure
  expect_s3_class(result$tidy_result, "tbl_df")
  expect_equal(nrow(result$tidy_result), ncol(test_gp_se))
  expect_true("tsne1" %in% names(result$tidy_result))
  expect_true("tsne2" %in% names(result$tidy_result))
  expect_true("sample" %in% names(result$tidy_result))

  # Check that coordinates are numeric
  expect_type(result$tidy_result$tsne1, "double")
  expect_type(result$tidy_result$tsne2, "double")

  # Check that no missing values
  expect_false(any(is.na(result$tidy_result$tsne1)))
  expect_false(any(is.na(result$tidy_result$tsne2)))

  # Check raw_result
  expect_s3_class(result$raw_result, "Rtsne")
})

test_that("gly_tsne works with custom parameters", {
  skip_if_not_installed("Rtsne")

  result <- gly_tsne(test_gp_se, perplexity = 2, max_iter = 250)

  expect_s3_class(result, c("glystats_tsne_res", "glystats_res"))
  expect_equal(nrow(result$tidy_result), ncol(test_gp_se))
  expect_true(all(c("tsne1", "tsne2", "sample") %in% names(result$tidy_result)))
})

test_that("gly_tsne handles perplexity adjustment", {
  skip_if_not_installed("Rtsne")

  # Should work and adjust perplexity automatically when too large
  expect_warning(
    result <- suppressMessages(gly_tsne(test_gp_se, perplexity = 15)),
    "Perplexity should be smaller"
  )

  expect_s3_class(result, c("glystats_tsne_res", "glystats_res"))
  expect_equal(nrow(result$tidy_result), ncol(test_gp_se))
})

test_that("gly_tsne works with default perplexity", {
  skip_if_not_installed("Rtsne")

  # Test with default perplexity (30) - should trigger warning and adjustment
  expect_warning(
    result <- suppressMessages(gly_tsne(test_gp_se)),
    "Perplexity should be smaller"
  )

  expect_s3_class(result, c("glystats_tsne_res", "glystats_res"))
  expect_equal(nrow(result$tidy_result), ncol(test_gp_se))
})

test_that("gly_tsne has consistent sample names", {
  skip_if_not_installed("Rtsne")

  result <- gly_tsne(test_gp_se, perplexity = 3)

  # Should have same sample names as input expression matrix
  expect_equal(
    sort(result$tidy_result$sample),
    sort(colnames(test_gp_se))
  )
})

test_that(".analyze_tsne works correctly", {
  skip_if_not_installed("Rtsne")

  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test function execution with appropriate perplexity
  suppressMessages({
    result <- .analyze_tsne(expr_mat, perplexity = 3)
  })

  # Verify results
  expect_s3_class(result, "glystats_tsne_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))
  expect_true("sample" %in% colnames(result$tidy_result))
  expect_true("tsne1" %in% colnames(result$tidy_result))
  expect_true("tsne2" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)
  expect_s3_class(result$raw_result, "Rtsne")
})
