test_that("gly_consensus_clustering works with basic parameters (default: cluster samples)", {
  # Use a subset of test data for faster testing
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10) |>
    glyexp::slice_sample_obs(n = 8)

  result <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    max_k = 4,
    reps = 50
  ))

  # Check that result is a list with expected structure
  expect_type(result, "list")
  expect_s3_class(result, "glystats_cc_res")
  expect_s3_class(result, "glystats_res")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check tidy_result is a tibble with expected structure
  tidy_result <- result$tidy_result
  expect_s3_class(tidy_result, "tbl_df")

  # Check columns (should cluster samples by default)
  expect_true("sample" %in% colnames(tidy_result))
  expect_true("k" %in% colnames(tidy_result))
  expect_true("cluster" %in% colnames(tidy_result))

  # Check that we have the right number of rows (number of samples * number of k values)
  # For max_k = 4, we have k = 2, 3, 4, so 3 k values * 8 samples = 24 rows
  expect_equal(nrow(tidy_result), 8 * 3)

  # Check that k values are correct
  expect_true(all(tidy_result$k %in% 2:4))
  # Check that cluster assignments are valid integers for each k
  for (k_val in 2:4) {
    k_subset <- tidy_result[tidy_result$k == k_val, ]
    expect_true(all(k_subset$cluster %in% 1:k_val))
  }

  # Check raw_result is a list (ConsensusClusterPlus result)
  expect_type(result$raw_result, "list")
})

test_that("gly_consensus_clustering works with on parameter", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 8) |>
    glyexp::slice_sample_obs(n = 6)

  # Test clustering samples (default)
  result_sample <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    on = "sample",
    max_k = 3,
    reps = 50
  ))
  expect_s3_class(result_sample, "glystats_cc_res")
  expect_true("sample" %in% colnames(result_sample$tidy_result))
  expect_equal(nrow(result_sample$tidy_result), 6 * 2) # 6 samples * 2 k values (k=2,3)

  # Test clustering variables
  result_var <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    on = "variable",
    max_k = 3,
    reps = 50
  ))
  expect_s3_class(result_var, "glystats_cc_res")
  expect_true("variable" %in% colnames(result_var$tidy_result))
  expect_equal(nrow(result_var$tidy_result), 8 * 2) # 8 variables * 2 k values (k=2,3)

  # Results should be different
  expect_false(identical(result_sample, result_var))
})

test_that("gly_consensus_clustering works with different parameters", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 8) # Increase sample size to allow max_k = 5

  # Test different max_k values
  result_max_k3 <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    max_k = 3,
    reps = 30
  ))
  result_max_k5 <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    max_k = 5,
    reps = 30
  ))

  expect_s3_class(result_max_k3, "glystats_cc_res")
  expect_s3_class(result_max_k5, "glystats_cc_res")

  # Different max_k should give different number of rows
  expect_equal(nrow(result_max_k3$tidy_result), 8 * 2) # k=2,3
  expect_equal(nrow(result_max_k5$tidy_result), 8 * 4) # k=2,3,4,5

  # Test different clustering algorithms
  result_hc <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    cluster_alg = "hc",
    max_k = 3,
    reps = 30
  ))
  result_km <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    cluster_alg = "km",
    max_k = 3,
    reps = 30
  ))

  expect_s3_class(result_hc, "glystats_cc_res")
  expect_s3_class(result_km, "glystats_cc_res")
})

test_that("gly_consensus_clustering scale parameter works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test with and without scaling
  result_scaled <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    scale = TRUE,
    max_k = 3,
    reps = 30
  ))
  result_unscaled <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    scale = FALSE,
    max_k = 3,
    reps = 30
  ))

  expect_s3_class(result_scaled, "glystats_cc_res")
  expect_s3_class(result_unscaled, "glystats_cc_res")

  # Both should have same structure
  expect_equal(
    nrow(result_scaled$tidy_result),
    nrow(result_unscaled$tidy_result)
  )
  expect_equal(
    colnames(result_scaled$tidy_result),
    colnames(result_unscaled$tidy_result)
  )
})

test_that("gly_consensus_clustering add_info parameter works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test add_info = FALSE (default clusters samples)
  result_no_info <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    add_info = FALSE,
    max_k = 3,
    reps = 30
  ))
  result_with_info <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    add_info = TRUE,
    max_k = 3,
    reps = 30
  ))

  # With add_info = FALSE, should have fewer columns
  expect_true(
    ncol(result_no_info$tidy_result) < ncol(result_with_info$tidy_result)
  )
  expect_true("sample" %in% colnames(result_with_info$tidy_result))

  # Test with variable clustering
  result_var_no_info <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    on = "variable",
    add_info = FALSE,
    max_k = 3,
    reps = 30
  ))
  result_var_with_info <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    on = "variable",
    add_info = TRUE,
    max_k = 3,
    reps = 30
  ))

  expect_true(
    ncol(result_var_no_info$tidy_result) <
      ncol(result_var_with_info$tidy_result)
  )
  expect_true("variable" %in% colnames(result_var_with_info$tidy_result))
})

test_that("gly_consensus_clustering returns both tidy and raw results", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test that function returns both tidy and raw results
  result <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    max_k = 3,
    reps = 30
  ))

  # Should return a list with both tidy_result and raw_result
  expect_type(result, "list")
  expect_s3_class(result, "glystats_cc_res")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # raw_result should have elements for each k value
  expect_true(length(result$raw_result) >= 3)
})

test_that("gly_consensus_clustering input validation works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  # Test invalid inputs
  expect_error(gly_consensus_clustering("not_an_experiment"))
  expect_error(gly_consensus_clustering(exp_subset, on = "invalid")) # Invalid on parameter
  expect_error(gly_consensus_clustering(exp_subset, max_k = 1)) # max_k must be >= 2
  expect_error(gly_consensus_clustering(exp_subset, reps = 0)) # reps must be >= 1
  expect_error(gly_consensus_clustering(exp_subset, p_item = 1.5)) # p_item must be <= 1
  expect_error(gly_consensus_clustering(exp_subset, cluster_alg = "invalid")) # Invalid algorithm
  expect_error(gly_consensus_clustering(exp_subset, scale = "yes"))
  expect_error(gly_consensus_clustering(exp_subset, add_info = "yes"))
})

test_that("gly_consensus_clustering_ works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test function execution (default: cluster samples)
  suppressMessages({
    result <- gly_consensus_clustering_(expr_mat, max_k = 3, reps = 30)
  })

  # Verify results
  expect_s3_class(result, "glystats_cc_res")
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  tidy_result <- result$tidy_result
  expect_true(tibble::is_tibble(tidy_result))
  expect_true("sample" %in% colnames(tidy_result))
  expect_true("k" %in% colnames(tidy_result))
  expect_true("cluster" %in% colnames(tidy_result))
  expect_equal(nrow(tidy_result), 10 * 2) # 10 samples * 2 k values (k=2,3)

  # Test clustering variables
  suppressMessages({
    result_var <- gly_consensus_clustering_(
      expr_mat,
      on = "variable",
      max_k = 3,
      reps = 30
    )
  })

  expect_true("variable" %in% colnames(result_var$tidy_result))
  expect_equal(nrow(result_var$tidy_result), 10 * 2) # 10 variables * 2 k values (k=2,3)
})

test_that("gly_consensus_clustering handles edge cases", {
  # Test with minimal data (5 samples - enough for clustering)
  exp_minimal <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 5)

  result <- suppressMessages(gly_consensus_clustering(
    exp_minimal,
    max_k = 3,
    reps = 20
  ))

  expect_s3_class(result, "glystats_cc_res")
  expect_equal(nrow(result$tidy_result), 5 * 2) # 5 samples * 2 k values (k=2,3)
  expect_true(all(result$tidy_result$k %in% 2:3))
})

test_that("gly_consensus_clustering long format structure is correct", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  result <- suppressMessages(gly_consensus_clustering(
    exp_subset,
    max_k = 4,
    reps = 30
  ))
  tidy_result <- result$tidy_result

  # Check that each sample appears for each k value
  samples <- unique(tidy_result$sample)
  k_values <- unique(tidy_result$k)

  expect_equal(length(samples), 6)
  expect_equal(length(k_values), 3) # k = 2, 3, 4
  expect_equal(nrow(tidy_result), 6 * 3) # 6 samples * 3 k values

  # Check that each sample-k combination appears exactly once
  for (sample in samples) {
    for (k in k_values) {
      sample_k_rows <- tidy_result[
        tidy_result$sample == sample & tidy_result$k == k,
      ]
      expect_equal(nrow(sample_k_rows), 1)
    }
  }
})

test_that("gly_consensus_clustering output_file parameter works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(50)) + 1, nrow = 5, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:5)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test output_file = NULL (should not generate any files)
  withr::with_tempdir({
    suppressMessages({
      result1 <- gly_consensus_clustering_(
        expr_mat,
        max_k = 3,
        reps = 10,
        output_file = NULL
      )
    })

    # Check that no PDF files are generated
    pdf_files <- list.files(pattern = "*.pdf")
    expect_equal(length(pdf_files), 0)

    # Check that result is still valid
    expect_s3_class(result1, "glystats_cc_res")
    expect_true("sample" %in% colnames(result1$tidy_result))
  })

  # Test output_file with specific path
  withr::with_tempdir({
    output_path <- "test_consensus.pdf"

    suppressMessages({
      result2 <- gly_consensus_clustering_(
        expr_mat,
        max_k = 3,
        reps = 10,
        output_file = output_path
      )
    })

    # Check that the specified file is created
    expect_true(file.exists(output_path))
    expect_true(file.size(output_path) > 0)

    # Check that no directories are created with the output file name
    expect_false(dir.exists(output_path))

    # Check that result is still valid
    expect_s3_class(result2, "glystats_cc_res")
    expect_true("sample" %in% colnames(result2$tidy_result))
  })

  # Test output_file with subdirectory
  withr::with_tempdir({
    output_path <- "results/consensus_plot.pdf"

    suppressMessages({
      result3 <- gly_consensus_clustering_(
        expr_mat,
        max_k = 3,
        reps = 10,
        output_file = output_path
      )
    })

    # Check that the file is created in the subdirectory
    expect_true(file.exists(output_path))
    expect_true(file.size(output_path) > 0)
    expect_true(dir.exists("results"))

    # Check that result is still valid
    expect_s3_class(result3, "glystats_cc_res")
    expect_true("sample" %in% colnames(result3$tidy_result))
  })
})
