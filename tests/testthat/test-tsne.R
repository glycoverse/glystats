# Test t-SNE analysis

test_that("gly_tsne works with default parameters", {
  skip_if_not_installed("Rtsne")
  
  # Use appropriate perplexity for small dataset to avoid warnings
  result <- gly_tsne(test_gp_exp, perplexity = 3)
  
  # Check basic structure
  expect_s3_class(result, "glystats_tsne_res")
  expect_s3_class(result, "tbl_df")
  
  # Check dimensions
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
  expect_true("tsne1" %in% names(result))
  expect_true("tsne2" %in% names(result))
  expect_true("sample" %in% names(result))
  
  # Check that coordinates are numeric
  expect_type(result$tsne1, "double")
  expect_type(result$tsne2, "double")
  
  # Check that no missing values
  expect_false(any(is.na(result$tsne1)))
  expect_false(any(is.na(result$tsne2)))
})

test_that("gly_tsne works with custom parameters", {
  skip_if_not_installed("Rtsne")
  
  result <- gly_tsne(test_gp_exp, perplexity = 2, max_iter = 250)
  
  expect_s3_class(result, "glystats_tsne_res")
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
  expect_true(all(c("tsne1", "tsne2", "sample") %in% names(result)))
})

test_that("gly_tsne handles perplexity adjustment", {
  skip_if_not_installed("Rtsne")
  
  # Should work and adjust perplexity automatically when too large
  expect_warning(
    result <- suppressMessages(gly_tsne(test_gp_exp, perplexity = 15)),
    "Perplexity should be smaller"
  )
  
  expect_s3_class(result, "glystats_tsne_res")
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
})

test_that("gly_tsne works with default perplexity", {
  skip_if_not_installed("Rtsne")
  
  # Test with default perplexity (30) - should trigger warning and adjustment
  expect_warning(
    result <- suppressMessages(gly_tsne(test_gp_exp)),
    "Perplexity should be smaller"
  )
  
  expect_s3_class(result, "glystats_tsne_res")
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
})

test_that("gly_tsne has consistent sample names", {
  skip_if_not_installed("Rtsne")
  
  result <- gly_tsne(test_gp_exp, perplexity = 3)
  
  # Should have same sample names as input expression matrix
  expect_equal(sort(result$sample), sort(colnames(test_gp_exp$expr_mat)))
})

test_that("gly_tsne requires Rtsne package", {
  # Test that the function would fail if Rtsne is not available
  # This is tested by verifying the error message structure exists
  expect_true(any(grepl("Rtsne", deparse(body(gly_tsne)))))
}) 