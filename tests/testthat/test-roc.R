test_that("gly_roc works with 2-group binary classification", {
  # Use test_gp_exp and filter to 2 groups for ROC analysis
  exp_2group <- exp_2groups() |>
    glyexp::slice_sample_var(n = 10) # Use smaller subset for faster testing

  # Run ROC analysis
  result <- suppressMessages(gly_roc(
    exp_2group,
    group_col = "group",
    pos_class = "H"
  ))

  # Test structure
  expect_s3_class(result, c("glystats_roc_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))

  # Test tidy_result structure
  expect_type(result$tidy_result, "list")
  expect_setequal(names(result$tidy_result), c("auc", "coords"))

  # Test AUC
  expect_s3_class(result$tidy_result$auc, "tbl_df")
  expect_true(all(
    c("variable", "auc", "auc_ci_low", "auc_ci_high") %in%
      colnames(result$tidy_result$auc)
  ))
  expect_true(all(
    result$tidy_result$auc$auc >= 0 & result$tidy_result$auc$auc <= 1
  )) # AUC should be between 0 and 1
  expect_true(all(
    result$tidy_result$auc$auc_ci_low >= 0 &
      result$tidy_result$auc$auc_ci_low <= 1
  ))
  expect_true(all(
    result$tidy_result$auc$auc_ci_high >= 0 &
      result$tidy_result$auc$auc_ci_high <= 1
  ))
  expect_true(all(
    result$tidy_result$auc$auc_ci_low <= result$tidy_result$auc$auc_ci_high
  ))

  # Test coords
  expect_s3_class(result$tidy_result$coords, "tbl_df")
  expect_true(all(
    c("variable", "threshold", "sensitivity", "specificity") %in%
      colnames(result$tidy_result$coords)
  ))
  expect_true(all(
    result$tidy_result$coords$sensitivity >= 0 &
      result$tidy_result$coords$sensitivity <= 1
  ))
  expect_true(all(
    result$tidy_result$coords$specificity >= 0 &
      result$tidy_result$coords$specificity <= 1
  ))
  expect_equal(length(unique(result$tidy_result$coords$variable)), 10)

  # Test raw_result
  expect_type(result$raw_result, "list")
  expect_equal(length(result$raw_result), 10)
  expect_true(all(purrr::map_lgl(result$raw_result, ~ inherits(.x, "roc"))))
})

test_that("gly_roc works with automatic pos_class detection", {
  # Use test_gp_exp and filter to 2 groups
  exp_2group <- exp_2groups() |>
    glyexp::slice_sample_var(n = 5)

  # Run ROC analysis without specifying pos_class
  result <- suppressMessages(gly_roc(exp_2group, group_col = "group"))

  # Test basic structure
  expect_s3_class(result, c("glystats_roc_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))
  expect_s3_class(result$tidy_result$coords, "tbl_df")
})

test_that("gly_roc error handling", {
  # Test with more than 2 groups
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)

  expect_error(suppressMessages(gly_roc(exp_3group)), "exactly 2 levels")

  # Test with 1 group
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)

  expect_error(suppressMessages(gly_roc(exp_1group)), "exactly 2 levels")

  # Test with non-existent group column
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)

  expect_error(
    suppressMessages(gly_roc(exp_2group, group_col = "nonexistent")),
    "not found in sample information"
  )

  # Test with invalid pos_class
  expect_error(
    suppressMessages(gly_roc(exp_2group, pos_class = "invalid")),
    "not found in group levels"
  )
})

test_that("gly_roc works with different group column names", {
  # Modify sample info to use different group column name
  exp_2group <- exp_2groups() |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::mutate_obs(condition = group)

  result <- suppressMessages(gly_roc(
    exp_2group,
    group_col = "condition",
    pos_class = "H"
  ))

  expect_s3_class(result, c("glystats_roc_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result", "meta_data"))
})

test_that("gly_roc assigns NA for failed variables", {
  # Use test_gp_exp and filter to 2 groups for ROC analysis
  exp_2group <- exp_2groups() |>
    glyexp::slice_sample_var(n = 10)
  exp_2group$expr_mat[1:3, ] <- NA # This will lead to pROC::roc() failing
  na_vars <- exp_2group$var_info$variable[1:3]

  # Run DEA with ROC analysis
  expect_warning(result <- suppressMessages(gly_roc(exp_2group)))

  # Test results
  expect_true(all(is.na(result$raw_result[na_vars])))
  auc_values <- result$tidy_result$auc |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(auc)
  expect_true(all(is.na(auc_values)))
  auc_ci_low_values <- result$tidy_result$auc |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(auc_ci_low)
  expect_true(all(is.na(auc_ci_low_values)))
  auc_ci_high_values <- result$tidy_result$auc |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(auc_ci_high)
  expect_true(all(is.na(auc_ci_high_values)))
  coords_values <- result$tidy_result$coords |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(sensitivity)
  expect_true(all(is.na(coords_values)))
})

test_that(".analyze_roc works correctly", {
  skip_if_not_installed("pROC")

  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B"), each = 5))

  # Test function execution
  suppressMessages({
    result <- .analyze_roc(expr_mat, groups, pos_class = "B")
  })

  # Verify results - .analyze_roc returns only tidy_result and raw_result (no meta_data)
  expect_s3_class(result, "glystats_roc_res")
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result$auc))
  expect_true(all(
    c("variable", "auc", "auc_ci_low", "auc_ci_high") %in%
      colnames(result$tidy_result$auc)
  ))
  expect_equal(nrow(result$tidy_result$auc), 10)
  expect_equal(length(result$raw_result), 10)
})
