test_that("gly_hclust does not leak graphics devices when ggdendro is available", {
  testthat::skip_if_not_installed("ggdendro")

  open_devices <- grDevices::dev.list()
  open_devices <- if (is.null(open_devices)) {
    integer(0)
  } else {
    as.integer(open_devices)
  }

  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 5) |>
    as_test_se()

  suppressMessages(gly_hclust(exp_subset, on = "sample", k_values = c(2)))

  after_devices <- grDevices::dev.list()
  after_devices <- if (is.null(after_devices)) {
    integer(0)
  } else {
    as.integer(after_devices)
  }

  expect_identical(after_devices, open_devices)
})

test_that("gly_hclust works with basic parameters (default: cluster variables)", {
  # Use a subset of test data for faster testing
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10) |>
    glyexp::slice_sample_obs(n = 8)

  result <- suppressMessages(gly_hclust(exp_subset))

  # Check that result is a list with expected components
  expect_type(result, "list")
  expect_s3_class(result, "glystats_hclust_res")
  expect_s3_class(result, "glystats_res")

  # Check that result has tidy_result and raw_result
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check raw_result is hclust object
  expect_s3_class(result$raw_result, "hclust")
  expect_true("height" %in% names(result$raw_result))
  expect_true("merge" %in% names(result$raw_result))
  expect_true("labels" %in% names(result$raw_result))

  # Check tidy_result structure
  expect_type(result$tidy_result, "list")

  # Check clusters component (should cluster variables by default)
  expect_true("clusters" %in% names(result$tidy_result))
  expect_s3_class(result$tidy_result$clusters, "tbl_df")
  expect_true("variable" %in% colnames(result$tidy_result$clusters))
  expect_true("cluster_k2" %in% colnames(result$tidy_result$clusters))
  expect_true("cluster_k3" %in% colnames(result$tidy_result$clusters))
  expect_true("cluster_k4" %in% colnames(result$tidy_result$clusters))
  expect_true("cluster_k5" %in% colnames(result$tidy_result$clusters))

  # Check heights component
  expect_true("heights" %in% names(result$tidy_result))
  expect_s3_class(result$tidy_result$heights, "tbl_df")
  expect_true("merge_step" %in% colnames(result$tidy_result$heights))
  expect_true("height" %in% colnames(result$tidy_result$heights))
  expect_true("n_clusters" %in% colnames(result$tidy_result$heights))

  # Check dendrogram and labels components exist
  expect_true("dendrogram" %in% names(result$tidy_result))
  expect_true("labels" %in% names(result$tidy_result))

  # Check that cluster assignments are valid
  expect_true(all(result$tidy_result$clusters$cluster_k2 %in% 1:2))
  expect_true(all(result$tidy_result$clusters$cluster_k3 %in% 1:3))
  expect_equal(nrow(result$tidy_result$clusters), 10) # Should match number of variables
})

test_that("gly_hclust works with on parameter", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 8) |>
    glyexp::slice_sample_obs(n = 6)

  # Test clustering variables (default)
  result_var <- suppressMessages(gly_hclust(
    exp_subset,
    on = "variable",
    k_values = c(2, 3)
  ))
  expect_s3_class(result_var, "glystats_hclust_res")
  expect_true("variable" %in% colnames(result_var$tidy_result$clusters))
  expect_equal(nrow(result_var$tidy_result$clusters), 8) # Number of variables

  # Test clustering samples
  result_sample <- suppressMessages(gly_hclust(
    exp_subset,
    on = "sample",
    k_values = c(2, 3)
  ))
  expect_s3_class(result_sample, "glystats_hclust_res")
  expect_true("sample" %in% colnames(result_sample$tidy_result$clusters))
  expect_equal(nrow(result_sample$tidy_result$clusters), 6) # Number of samples

  # Results should be different
  expect_false(identical(
    result_var$tidy_result$clusters,
    result_sample$tidy_result$clusters
  ))
})

test_that("gly_hclust works with different methods", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test different clustering methods
  methods <- c("complete", "average", "single", "ward.D2")

  for (method in methods) {
    # Pass hclust method via hclust.method parameter
    result <- suppressMessages(gly_hclust(exp_subset, hclust.method = method))
    expect_s3_class(result, "glystats_hclust_res")
    expect_true("clusters" %in% names(result$tidy_result))
  }
})

test_that("gly_hclust works with different distance methods", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test different distance methods
  dist_methods <- c("euclidean", "manhattan", "maximum")

  for (dist_method in dist_methods) {
    # Pass distance method via dist.method parameter
    result <- suppressMessages(gly_hclust(
      exp_subset,
      dist.method = dist_method
    ))
    expect_s3_class(result, "glystats_hclust_res")
    expect_true("clusters" %in% names(result$tidy_result))
  }
})

test_that("gly_hclust works with custom k_values", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 8)

  # Test with custom k values
  result <- suppressMessages(gly_hclust(exp_subset, k_values = c(2, 4)))

  expect_true("cluster_k2" %in% colnames(result$tidy_result$clusters))
  expect_true("cluster_k4" %in% colnames(result$tidy_result$clusters))
  expect_false("cluster_k3" %in% colnames(result$tidy_result$clusters))
  expect_false("cluster_k5" %in% colnames(result$tidy_result$clusters))
})

test_that("gly_hclust works with k_values = NULL", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test with no cluster assignments
  result <- suppressMessages(gly_hclust(exp_subset, k_values = NULL))

  expect_false("clusters" %in% names(result$tidy_result))
  expect_true("heights" %in% names(result$tidy_result))
  expect_true("dendrogram" %in% names(result$tidy_result))
  expect_true("labels" %in% names(result$tidy_result))
})


test_that("gly_hclust scale parameter works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test with and without scaling
  result_scaled <- suppressMessages(gly_hclust(exp_subset, scale = TRUE))
  result_unscaled <- suppressMessages(gly_hclust(exp_subset, scale = FALSE))

  expect_s3_class(result_scaled, "glystats_hclust_res")
  expect_s3_class(result_unscaled, "glystats_hclust_res")

  # Results should be different when scaling is applied
  expect_false(identical(
    result_scaled$tidy_result$heights,
    result_unscaled$tidy_result$heights
  ))
})

test_that("gly_hclust add_info parameter works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test add_info = FALSE (default clusters variables)
  result_no_info <- suppressMessages(gly_hclust(exp_subset, add_info = FALSE))
  result_with_info <- suppressMessages(gly_hclust(exp_subset, add_info = TRUE))

  # With add_info = FALSE, should have fewer columns in clusters tibble
  expect_true(
    ncol(result_no_info$tidy_result$clusters) <
      ncol(result_with_info$tidy_result$clusters)
  )
  expect_true("variable" %in% colnames(result_with_info$tidy_result$clusters))

  # Test with sample clustering
  result_sample_no_info <- suppressMessages(gly_hclust(
    exp_subset,
    on = "sample",
    add_info = FALSE
  ))
  result_sample_with_info <- suppressMessages(gly_hclust(
    exp_subset,
    on = "sample",
    add_info = TRUE
  ))

  expect_true(
    ncol(result_sample_no_info$tidy_result$clusters) <
      ncol(result_sample_with_info$tidy_result$clusters)
  )
  expect_true(
    "sample" %in% colnames(result_sample_with_info$tidy_result$clusters)
  )
})

test_that("gly_hclust handles edge cases", {
  # Test with minimal data (just 3 samples - minimum for clustering)
  exp_minimal <- test_gp_exp |>
    glyexp::slice_sample_var(n = 3) |>
    glyexp::slice_sample_obs(n = 3)

  result <- suppressMessages(gly_hclust(exp_minimal, k_values = c(2)))

  expect_s3_class(result, "glystats_hclust_res")
  expect_true("clusters" %in% names(result$tidy_result))
  expect_equal(nrow(result$tidy_result$clusters), 3)
  expect_true(all(result$tidy_result$clusters$cluster_k2 %in% 1:2))
})

test_that("gly_hclust input validation works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test invalid inputs
  expect_error(gly_hclust("not_an_experiment"))
  expect_error(gly_hclust(exp_subset, on = "invalid")) # Invalid on parameter
  expect_error(gly_hclust(exp_subset, method = 123))
  expect_error(gly_hclust(exp_subset, k_values = c(1))) # k must be >= 2
  expect_error(gly_hclust(exp_subset, scale = "yes"))
  expect_error(gly_hclust(exp_subset, add_info = "yes"))
})

test_that(".analyze_hclust works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test function execution
  suppressMessages({
    result <- .analyze_hclust(expr_mat)
  })

  # Verify results
  expect_s3_class(result, "glystats_hclust_res")
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))
  expect_true("clusters" %in% names(result$tidy_result))
  expect_true("heights" %in% names(result$tidy_result))
  expect_true(tibble::is_tibble(result$tidy_result$clusters))
  expect_equal(nrow(result$tidy_result$clusters), 10)
  expect_s3_class(result$raw_result, "hclust")
})
