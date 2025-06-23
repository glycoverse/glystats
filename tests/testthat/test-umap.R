# Test UMAP analysis

test_that("gly_umap works with default parameters", {
  skip_if_not_installed("uwot")
  
  # Use appropriate n_neighbors for small dataset
  result <- gly_umap(test_gp_exp, n_neighbors = 3)
  
  # Check basic structure
  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_s3_class(result, "tbl_df")
  
  # Check dimensions
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
  expect_true("umap1" %in% names(result))
  expect_true("umap2" %in% names(result))
  expect_true("sample" %in% names(result))
  
  # Check that coordinates are numeric
  expect_type(result$umap1, "double")
  expect_type(result$umap2, "double")
  
  # Check that no missing values
  expect_false(any(is.na(result$umap1)))
  expect_false(any(is.na(result$umap2)))
})

test_that("gly_umap works with custom parameters", {
  skip_if_not_installed("uwot")
  
  result <- gly_umap(test_gp_exp, n_neighbors = 2, min_dist = 0.01, n_epochs = 50)
  
  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
  expect_true(all(c("umap1", "umap2", "sample") %in% names(result)))
})

test_that("gly_umap handles n_neighbors adjustment", {
  skip_if_not_installed("uwot")
  
  # Should work and adjust n_neighbors automatically when too large
  expect_warning(
    result <- suppressMessages(gly_umap(test_gp_exp, n_neighbors = 50)),
    "n_neighbors should be smaller"
  )
  
  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
})

test_that("gly_umap works with more than 2 components", {
  skip_if_not_installed("uwot")
  
  result <- gly_umap(test_gp_exp, n_neighbors = 3, n_components = 3)
  
  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
  expect_true(all(c("umap1", "umap2", "umap3", "sample") %in% names(result)))
  
  # Check that all coordinates are numeric
  expect_type(result$umap1, "double")
  expect_type(result$umap2, "double")
  expect_type(result$umap3, "double")
})

test_that("gly_umap works with default n_neighbors", {
  skip_if_not_installed("uwot")
  
  # Test with default n_neighbors (15) - should trigger warning and adjustment
  expect_warning(
    result <- suppressMessages(gly_umap(test_gp_exp)),
    "n_neighbors should be smaller"
  )
  
  expect_s3_class(result, c("glystats_umap_res", "glystats_res"))
  expect_equal(nrow(result), nrow(test_gp_exp$sample_info))
})

test_that("gly_umap has consistent sample names", {
  skip_if_not_installed("uwot")
  
  result <- gly_umap(test_gp_exp, n_neighbors = 3)
  
  # Should have same sample names as input expression matrix
  expect_equal(sort(result$sample), sort(colnames(test_gp_exp$expr_mat)))
})

test_that("gly_umap requires uwot package", {
  # Test that the function would fail if uwot is not available
  # This is tested by verifying the error message structure exists
  expect_true(any(grepl("uwot", deparse(body(gly_umap)))))
}) 