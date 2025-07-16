test_that("gly_cor basic functionality works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::slice_sample_obs(n = 6)

  result <- suppressMessages(gly_cor(exp_subset))

  # Check that result is a tibble with expected structure
  expect_s3_class(result, "tbl_df")
  expect_s3_class(result, "glystats_cor_res")
  expect_s3_class(result, "glystats_res")

  # Check columns (should correlate variables by default)
  expect_true("variable1" %in% colnames(result))
  expect_true("variable2" %in% colnames(result))
  expect_true("cor" %in% colnames(result))
  expect_true("p_value" %in% colnames(result))
  expect_true("p_adj" %in% colnames(result))

  # Check that we have the right number of pairs (n choose 2)
  n_vars <- 5
  expected_pairs <- n_vars * (n_vars - 1) / 2
  expect_equal(nrow(result), expected_pairs)

  # Check that correlation values are in valid range
  expect_true(all(result$cor >= -1 & result$cor <= 1))
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
  expect_true(all(result$p_adj >= 0 & result$p_adj <= 1))
})

test_that("gly_cor on parameter works correctly", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 6) |>  # Need >4 variables for sample correlation
    glyexp::slice_sample_obs(n = 6)     # Need >4 samples for variable correlation

  # Test correlating variables (default)
  result_var <- suppressMessages(gly_cor(exp_subset, on = "variable"))
  expect_s3_class(result_var, "glystats_cor_res")
  expect_true("variable1" %in% colnames(result_var))
  expect_true("variable2" %in% colnames(result_var))
  expect_equal(nrow(result_var), 6 * 5 / 2)  # 6 choose 2

  # Test correlating samples
  result_sample <- suppressMessages(gly_cor(exp_subset, on = "sample"))
  expect_s3_class(result_sample, "glystats_cor_res")
  expect_true("sample1" %in% colnames(result_sample))
  expect_true("sample2" %in% colnames(result_sample))
  expect_equal(nrow(result_sample), 6 * 5 / 2)  # 6 choose 2

  # Results should be different
  expect_false(identical(result_var$cor, result_sample$cor))
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
    expect_true("cor" %in% colnames(result))
    expect_true(all(result$cor >= -1 & result$cor <= 1))
  }
})

test_that("gly_cor p_adj_method parameter works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Test with p-value adjustment
  result_adj <- suppressMessages(gly_cor(exp_subset, p_adj_method = "BH"))
  expect_true("p_adj" %in% colnames(result_adj))

  # Test without p-value adjustment
  result_no_adj <- suppressMessages(gly_cor(exp_subset, p_adj_method = NULL))
  expect_false("p_adj" %in% colnames(result_no_adj))
  expect_true("p_value" %in% colnames(result_no_adj))
})

test_that("gly_cor return_raw parameter works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Test return_raw = TRUE
  result_raw <- suppressMessages(gly_cor(exp_subset, return_raw = TRUE))

  expect_type(result_raw, "list")
  expect_true("r" %in% names(result_raw))  # correlation matrix in rcorr
  expect_true("P" %in% names(result_raw))  # p-value matrix in rcorr
  expect_true(is.matrix(result_raw$r))
  expect_true(is.matrix(result_raw$P))
  expect_equal(dim(result_raw$r), c(4, 4))
  expect_equal(dim(result_raw$P), c(4, 4))
})



test_that("gly_cor handles edge cases", {
  # Test with minimal data (3 variables, 6 observations - Hmisc::rcorr requires >4 observations)
  exp_minimal <- test_gp_exp |>
    glyexp::slice_sample_var(n = 3) |>
    glyexp::slice_sample_obs(n = 6)

  result <- suppressMessages(gly_cor(exp_minimal))

  expect_s3_class(result, "glystats_cor_res")
  expect_equal(nrow(result), 3)  # 3 choose 2
  expect_true(all(result$cor >= -1 & result$cor <= 1))
})

test_that("gly_cor input validation works", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Test invalid inputs
  expect_error(gly_cor("not_an_experiment"))
  expect_error(gly_cor(exp_subset, on = "invalid"))  # Invalid on parameter
  expect_error(gly_cor(exp_subset, method = "invalid"))  # Invalid method
  expect_error(gly_cor(exp_subset, p_adj_method = "invalid"))  # Invalid p_adj_method
  expect_error(gly_cor(exp_subset, return_raw = "yes"))
})

test_that("gly_cor produces consistent results", {
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 4) |>
    glyexp::slice_sample_obs(n = 6)

  # Run the same analysis twice
  result1 <- suppressMessages(gly_cor(exp_subset))
  result2 <- suppressMessages(gly_cor(exp_subset))

  # Results should be identical
  expect_equal(result1$cor, result2$cor)
  expect_equal(result1$p_value, result2$p_value)
  expect_equal(result1$p_adj, result2$p_adj)
})
