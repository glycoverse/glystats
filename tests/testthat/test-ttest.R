test_that("gly_ttest works with t-test method", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10) |>
    as_test_se() # Use smaller subset for faster testing

  # Run DEA with t-test
  result <- suppressMessages(gly_ttest(exp_2group))

  # Test core functionality
  expect_s3_class(result, c("glystats_ttest_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))
  expect_equal(nrow(result$tidy_result), 10)
  expect_true("p_adj" %in% colnames(result$tidy_result)) # p_adj should exist
  expect_true("log2fc" %in% colnames(result$tidy_result)) # log2fc should exist
  expect_true("effect_size" %in% colnames(result$tidy_result))
  expect_type(result$tidy_result$log2fc, "double") # log2fc should be numeric
  expect_type(result$tidy_result$effect_size, "double")

  # Test raw_result
  expect_type(result$raw_result, "list")
  expect_equal(length(result$raw_result), 10)
  expect_true(all(purrr::map_lgl(result$raw_result, ~ inherits(.x, "htest"))))
})

test_that("gly_ttest assigns NA for failed variables", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)
  exp_2group$expr_mat[1:3, ] <- NA # This will lead to stats::t.test() failing
  na_vars <- exp_2group$var_info$variable[1:3]

  # Run DEA with t-test
  expect_warning(result <- suppressMessages(gly_ttest(exp_2group)))

  # Test results
  main_test_raw <- result$raw_result
  expect_true(all(is.na(main_test_raw[na_vars])))
  main_test_tidy <- result$tidy_result
  p_values <- main_test_tidy |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(p_val)
  expect_true(all(is.na(p_values)))
})

test_that("gly_wilcox works with wilcoxon method", {
  # Use test_gp_exp and filter to 2 groups for wilcoxon test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 10) |>
    as_test_se() # Use smaller subset for faster testing

  # Run DEA with wilcoxon test
  result <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group)))

  # Test core functionality
  expect_s3_class(result, c("glystats_wilcox_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))
  expect_equal(nrow(result$tidy_result), 10)
  expect_true("log2fc" %in% colnames(result$tidy_result)) # Wilcoxon should now have log2fc
  expect_true("effect_size" %in% colnames(result$tidy_result))
  expect_type(result$tidy_result$log2fc, "double") # log2fc should be numeric
  expect_type(result$tidy_result$effect_size, "double")

  # Test raw_result
  expect_type(result$raw_result, "list")
  expect_equal(length(result$raw_result), 10)
  expect_true(all(purrr::map_lgl(result$raw_result, ~ inherits(.x, "htest"))))
})

test_that(".analyze_ttest works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B"), each = 5))

  # Test function execution
  suppressMessages({
    result <- .analyze_ttest(expr_mat, groups)
  })

  # Verify results
  expect_s3_class(result, "glystats_ttest_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_true(nrow(result$tidy_result) > 0)
  expect_true("log2fc" %in% colnames(result$tidy_result))
  expect_true("effect_size" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)
})

test_that(".analyze_ttest returns Cohen's d in effect_size", {
  expr_mat <- matrix(
    c(1, 2, 3, 8, 9, 10),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:6))
  )
  groups <- factor(c("A", "A", "A", "B", "B", "B"))

  result <- suppressMessages(.analyze_ttest(
    expr_mat,
    groups,
    p_adj_method = NULL
  ))

  log_values <- log2(c(1, 2, 3, 8, 9, 10) + 1e-6)
  group_a <- log_values[1:3]
  group_b <- log_values[4:6]
  pooled_sd <- sqrt(
    (((length(group_a) - 1) * stats::sd(group_a)^2) +
      ((length(group_b) - 1) * stats::sd(group_b)^2)) /
      (length(group_a) + length(group_b) - 2)
  )
  expected <- (mean(group_b) - mean(group_a)) / pooled_sd

  expect_equal(result$tidy_result$effect_size, expected, tolerance = 1e-10)
})

test_that(".log_transform_expr_mat uses a small pseudo-count", {
  expr_mat <- matrix(
    c(0, 1, 10),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:3))
  )

  result <- .log_transform_expr_mat(expr_mat)

  expect_equal(result, log2(expr_mat + 1e-6), tolerance = 1e-10)
})

test_that(".analyze_ttest direction matches log2fc and ref_group", {
  expr_mat <- matrix(
    c(10, 11, 12, 13, 14, 1, 2, 3, 4, 5),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:10))
  )
  groups <- factor(c(rep("A", 5), rep("B", 5)))

  result_default <- suppressMessages(.analyze_ttest(
    expr_mat,
    groups,
    p_adj_method = NULL
  ))
  result_ref_b <- suppressMessages(.analyze_ttest(
    expr_mat,
    groups,
    p_adj_method = NULL,
    ref_group = "B"
  ))

  expect_lt(result_default$tidy_result$log2fc, 0)
  expect_lt(result_default$tidy_result$estimate, 0)
  expect_lt(result_default$tidy_result$statistic, 0)
  expect_lt(result_default$tidy_result$conf_low, 0)
  expect_lt(result_default$tidy_result$conf_high, 0)

  expect_equal(
    result_default$tidy_result$estimate1,
    result_ref_b$tidy_result$estimate2,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$estimate2,
    result_ref_b$tidy_result$estimate1,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$log2fc,
    -result_ref_b$tidy_result$log2fc,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$estimate,
    -result_ref_b$tidy_result$estimate,
    tolerance = 1e-4
  )
  expect_equal(
    result_default$tidy_result$statistic,
    -result_ref_b$tidy_result$statistic,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$conf_low,
    -result_ref_b$tidy_result$conf_high,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$conf_high,
    -result_ref_b$tidy_result$conf_low,
    tolerance = 1e-10
  )
})

test_that(".analyze_ttest reports one-sided alternatives in the output direction", {
  expr_mat <- matrix(
    c(10, 11, 12, 13, 14, 1, 2, 3, 4, 5),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:10))
  )
  groups <- factor(c(rep("A", 5), rep("B", 5)))

  result_greater <- suppressMessages(.analyze_ttest(
    expr_mat,
    groups,
    p_adj_method = NULL,
    alternative = "greater"
  ))
  result_less <- suppressMessages(.analyze_ttest(
    expr_mat,
    groups,
    p_adj_method = NULL,
    alternative = "less"
  ))

  expect_identical(result_greater$tidy_result$alternative, "greater")
  expect_identical(result_less$tidy_result$alternative, "less")
  expect_gt(result_greater$tidy_result$p_val, 0.9)
  expect_lt(result_less$tidy_result$p_val, 0.1)
})

test_that(".analyze_wilcox direction matches log2fc and ref_group", {
  expr_mat <- matrix(
    c(10, 11, 12, 13, 14, 1, 2, 3, 4, 5),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:10))
  )
  groups <- factor(c(rep("A", 5), rep("B", 5)))

  result_default <- suppressMessages(suppressWarnings(.analyze_wilcox(
    expr_mat,
    groups,
    p_adj_method = NULL,
    conf.int = TRUE,
    exact = FALSE
  )))
  result_ref_b <- suppressMessages(suppressWarnings(.analyze_wilcox(
    expr_mat,
    groups,
    p_adj_method = NULL,
    ref_group = "B",
    conf.int = TRUE,
    exact = FALSE
  )))

  expect_lt(result_default$tidy_result$log2fc, 0)
  expect_lt(result_default$tidy_result$estimate, 0)
  expect_lt(
    result_default$tidy_result$statistic,
    result_ref_b$tidy_result$statistic
  )
  expect_lt(result_default$tidy_result$conf_low, 0)
  expect_lt(result_default$tidy_result$conf_high, 0)

  expect_equal(
    result_default$tidy_result$log2fc,
    -result_ref_b$tidy_result$log2fc,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$estimate,
    -result_ref_b$tidy_result$estimate,
    tolerance = 1e-4
  )
  expect_equal(
    result_default$tidy_result$conf_low,
    -result_ref_b$tidy_result$conf_high,
    tolerance = 1e-4
  )
  expect_equal(
    result_default$tidy_result$conf_high,
    -result_ref_b$tidy_result$conf_low,
    tolerance = 1e-4
  )
})

test_that(".analyze_wilcox preserves one-sided alternatives in the output direction", {
  expr_mat <- matrix(
    c(10, 11, 12, 13, 14, 1, 2, 3, 4, 5),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:10))
  )
  groups <- factor(c(rep("A", 5), rep("B", 5)))

  result_greater <- suppressMessages(suppressWarnings(.analyze_wilcox(
    expr_mat,
    groups,
    p_adj_method = NULL,
    alternative = "greater",
    conf.int = TRUE,
    exact = FALSE
  )))
  result_less <- suppressMessages(suppressWarnings(.analyze_wilcox(
    expr_mat,
    groups,
    p_adj_method = NULL,
    alternative = "less",
    conf.int = TRUE,
    exact = FALSE
  )))

  expect_identical(result_greater$tidy_result$alternative, "greater")
  expect_identical(result_less$tidy_result$alternative, "less")
  expect_gt(result_greater$tidy_result$p_val, 0.9)
  expect_lt(result_less$tidy_result$p_val, 0.1)
})

test_that(".analyze_wilcox works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B"), each = 5))

  # Test function execution
  suppressMessages({
    result <- .analyze_wilcox(expr_mat, groups)
  })

  # Verify results
  expect_s3_class(result, "glystats_wilcox_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_true(nrow(result$tidy_result) > 0)
  expect_true("log2fc" %in% colnames(result$tidy_result))
  expect_true("effect_size" %in% colnames(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)
})

test_that(".analyze_wilcox returns rank-biserial correlation in effect_size", {
  expr_mat <- matrix(
    c(1, 2, 3, 8, 9, 10),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:6))
  )
  groups <- factor(c("A", "A", "A", "B", "B", "B"))

  result <- suppressMessages(suppressWarnings(.analyze_wilcox(
    expr_mat,
    groups,
    p_adj_method = NULL,
    exact = FALSE
  )))

  log_values <- log2(c(1, 2, 3, 8, 9, 10) + 1e-6)
  ranks <- rank(log_values)
  n_ref <- 3
  n_test <- 3
  w_stat <- sum(ranks[4:6])
  u_stat <- w_stat - n_test * (n_test + 1) / 2
  expected <- (2 * u_stat / (n_ref * n_test)) - 1

  expect_equal(result$tidy_result$effect_size, expected, tolerance = 1e-10)
})

test_that("gly_ttest and gly_wilcox basic functionality works", {
  # Test both methods work with test_gp_exp
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5) # Use very small subset

  # 2-group methods
  exp_2group <- exp_small |> glyexp::filter_obs(group %in% c("C", "H"))
  expect_no_error(suppressMessages(gly_ttest(exp_2group)))
  expect_no_error(suppressMessages(suppressWarnings(gly_wilcox(exp_2group))))
})

test_that("gly_ttest and gly_wilcox error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)

  # Test various error conditions - group column not found
  expect_error(
    suppressMessages(gly_ttest(exp_small, group_col = "nonexistent")),
    "not found in sample information"
  )
  expect_error(
    suppressMessages(gly_wilcox(exp_small, group_col = "nonexistent")),
    "not found in sample information"
  )
})

test_that("gly_ttest and gly_wilcox group validation", {
  # Test with 3 groups (using C, H, M from test_gp_exp)
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with 1 group
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)

  # Test 3 groups with 2-group methods
  expect_error(suppressMessages(gly_ttest(exp_3group)), "exactly 2 levels")
  expect_error(suppressMessages(gly_wilcox(exp_3group)), "exactly 2 levels")

  # Test 1 group with 2-group methods
  expect_error(suppressMessages(gly_ttest(exp_1group)), "exactly 2 levels")
  expect_error(suppressMessages(gly_wilcox(exp_1group)), "exactly 2 levels")
})

test_that("gly_ttest ref_group parameter works", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with default (no ref_group)
  result_default <- suppressMessages(gly_ttest(exp_2group))

  # Test with ref_group = "H" (should reverse the comparison)
  result_ref_h <- suppressMessages(gly_ttest(exp_2group, ref_group = "H"))

  # Test with ref_group = "C" (should be same as default since C is first alphabetically)
  result_ref_c <- suppressMessages(gly_ttest(exp_2group, ref_group = "C"))

  # Check that log2fc values are negated when reference group changes
  expect_equal(
    result_default$tidy_result$log2fc,
    -result_ref_h$tidy_result$log2fc,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$log2fc,
    result_ref_c$tidy_result$log2fc,
    tolerance = 1e-10
  )

  # Test invalid ref_group
  expect_error(
    suppressMessages(gly_ttest(exp_2group, ref_group = "invalid")),
    "Must be element of set"
  )
})

test_that("gly_wilcox ref_group parameter works", {
  # Use test_gp_exp and filter to 2 groups for wilcoxon test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with default (no ref_group)
  result_default <- suppressMessages(suppressWarnings(gly_wilcox(exp_2group)))

  # Test with ref_group = "Y" (should reverse the comparison)
  result_ref_y <- suppressMessages(suppressWarnings(gly_wilcox(
    exp_2group,
    ref_group = "Y"
  )))

  # Test with ref_group = "M" (should be same as default since M is first alphabetically)
  result_ref_m <- suppressMessages(suppressWarnings(gly_wilcox(
    exp_2group,
    ref_group = "M"
  )))

  # Check that log2fc values are negated when reference group changes
  expect_equal(
    result_default$tidy_result$log2fc,
    -result_ref_y$tidy_result$log2fc,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$log2fc,
    result_ref_m$tidy_result$log2fc,
    tolerance = 1e-10
  )

  # Test invalid ref_group
  expect_error(
    suppressMessages(gly_wilcox(exp_2group, ref_group = "invalid")),
    "Must be element of set"
  )
})
