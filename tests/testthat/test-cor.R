test_that("gly_cor basic functionality works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  result <- suppressMessages(gly_cor(exp_subset))

  # Check that result is a list with expected structure
  expect_type(result, "list")
  expect_s3_class(result, "glystats_cor_res")
  expect_s3_class(result, "glystats_res")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check tidy_result structure
  tidy_result <- result$tidy_result
  expect_s3_class(tidy_result, "tbl_df")
  expect_true("variable1" %in% colnames(tidy_result))
  expect_true("variable2" %in% colnames(tidy_result))
  expect_true("cor" %in% colnames(tidy_result))
  expect_true("p_val" %in% colnames(tidy_result))
  expect_true("p_adj" %in% colnames(tidy_result))

  # Check that we have the right number of pairs (n choose 2)
  n_vars <- 5
  expected_pairs <- n_vars * (n_vars - 1) / 2
  expect_equal(nrow(tidy_result), expected_pairs)

  # Check that correlation values are in valid range
  expect_true(all(tidy_result$cor >= -1 & tidy_result$cor <= 1))
  expect_true(all(tidy_result$p_val >= 0 & tidy_result$p_val <= 1))
  expect_true(all(tidy_result$p_adj >= 0 & tidy_result$p_adj <= 1))

  # Check raw_result structure
  raw_result <- result$raw_result
  expect_type(raw_result, "list")
  expect_true("r" %in% names(raw_result)) # correlation matrix in rcorr
  expect_true("P" %in% names(raw_result)) # p-value matrix in rcorr
  expect_true(is.matrix(raw_result$r))
  expect_true(is.matrix(raw_result$P))
})

test_that("gly_cor on parameter works correctly", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 6) |> # Need >4 variables for sample correlation
    glyexp::slice_sample_obs(n = 6) # Need >4 samples for variable correlation

  # Test correlating variables (default)
  result_var <- suppressMessages(gly_cor(exp_subset, on = "variable"))
  expect_s3_class(result_var, "glystats_cor_res")
  expect_true("variable1" %in% colnames(result_var$tidy_result))
  expect_true("variable2" %in% colnames(result_var$tidy_result))
  expect_equal(nrow(result_var$tidy_result), 6 * 5 / 2) # 6 choose 2

  # Test correlating samples
  result_sample <- suppressMessages(gly_cor(exp_subset, on = "sample"))
  expect_s3_class(result_sample, "glystats_cor_res")
  expect_true("sample1" %in% colnames(result_sample$tidy_result))
  expect_true("sample2" %in% colnames(result_sample$tidy_result))
  expect_equal(nrow(result_sample$tidy_result), 6 * 5 / 2) # 6 choose 2

  # Results should be different
  expect_false(identical(
    result_var$tidy_result$cor,
    result_sample$tidy_result$cor
  ))
})

test_that("gly_cor works with different correlation methods", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Test different correlation methods
  methods <- c("pearson", "spearman")

  for (method in methods) {
    result <- suppressMessages(gly_cor(exp_subset, method = method))
    expect_s3_class(result, "glystats_cor_res")
    expect_true("cor" %in% colnames(result$tidy_result))
    expect_true(all(result$tidy_result$cor >= -1 & result$tidy_result$cor <= 1))
  }
})

test_that("gly_cor p_adj_method parameter works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Test with p-value adjustment
  result_adj <- suppressMessages(gly_cor(exp_subset, p_adj_method = "BH"))
  expect_true("p_adj" %in% colnames(result_adj$tidy_result))

  # Test without p-value adjustment
  result_no_adj <- suppressMessages(gly_cor(exp_subset, p_adj_method = NULL))
  expect_false("p_adj" %in% colnames(result_no_adj$tidy_result))
  expect_true("p_val" %in% colnames(result_no_adj$tidy_result))
})

test_that("gly_cor raw_result contains expected structure", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Test that raw_result is included in the output
  result <- suppressMessages(gly_cor(exp_subset))

  expect_type(result$raw_result, "list")
  expect_true("r" %in% names(result$raw_result)) # correlation matrix in rcorr
  expect_true("P" %in% names(result$raw_result)) # p-value matrix in rcorr
  expect_true(is.matrix(result$raw_result$r))
  expect_true(is.matrix(result$raw_result$P))
  expect_equal(dim(result$raw_result$r), c(4, 4))
  expect_equal(dim(result$raw_result$P), c(4, 4))
})


test_that("gly_cor handles edge cases", {
  # Test with minimal data (3 variables, 6 observations - Hmisc::rcorr requires >4 observations)
  exp_minimal <- test_gp_exp |>
    glyexp::slice_sample_var(n = 3) |>
    glyexp::slice_sample_obs(n = 6)

  result <- suppressMessages(gly_cor(exp_minimal))

  expect_s3_class(result, "glystats_cor_res")
  expect_equal(nrow(result$tidy_result), 3) # 3 choose 2
  expect_true(all(result$tidy_result$cor >= -1 & result$tidy_result$cor <= 1))
})

test_that("gly_cor input validation works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Test invalid inputs
  expect_error(gly_cor("not_an_experiment"))
  expect_error(gly_cor(exp_subset, on = "invalid")) # Invalid on parameter
  expect_error(gly_cor(exp_subset, method = "invalid")) # Invalid method
  expect_error(gly_cor(exp_subset, p_adj_method = "invalid")) # Invalid p_adj_method
})

test_that("gly_cor produces consistent results", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Run the same analysis twice
  result1 <- suppressMessages(gly_cor(exp_subset))
  result2 <- suppressMessages(gly_cor(exp_subset))

  # Results should be identical
  expect_equal(result1$tidy_result$cor, result2$tidy_result$cor)
  expect_equal(result1$tidy_result$p_val, result2$tidy_result$p_val)
  expect_equal(result1$tidy_result$p_adj, result2$tidy_result$p_adj)
})

test_that(".analyze_cor works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Test function execution
  suppressMessages({
    result <- .analyze_cor(expr_mat)
  })

  # Verify results
  expect_s3_class(result, "glystats_cor_res")
  expect_true(nrow(result$tidy_result) > 0)
  expect_true("variable1" %in% colnames(result$tidy_result))
  expect_true("variable2" %in% colnames(result$tidy_result))
  expect_true("cor" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 45) # 10 choose 2 = 45 pairs
})
