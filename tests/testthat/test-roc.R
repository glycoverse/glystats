test_that("gly_roc works with 2-group binary classification", {
  # Use test_gp_exp and filter to 2 groups for ROC analysis
  exp_2group <- exp_2groups() |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing
  
  # Run ROC analysis
  result <- suppressMessages(gly_roc(exp_2group, group_col = "group", pos_class = "H"))
  
  # Test structure
  expect_s3_class(result, c("glystats_roc_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("auc", "coords"))
  
  # Test AUC
  expect_s3_class(result$auc, "tbl_df")
  expect_true(all(result$auc$auc >= 0 & result$auc$auc <= 1))  # AUC should be between 0 and 1

  # Test coords
  expect_s3_class(result$coords, "tbl_df")
  expect_true(all(c("variable", "threshold", "sensitivity", "specificity") %in% colnames(result$coords)))
  expect_true(all(result$coords$sensitivity >= 0 & result$coords$sensitivity <= 1))
  expect_true(all(result$coords$specificity >= 0 & result$coords$specificity <= 1))
  expect_equal(length(unique(result$coords$variable)), 10)
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
  expect_setequal(names(result), c("auc", "coords"))
  expect_s3_class(result$coords, "tbl_df")
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
  
  expect_error(suppressMessages(gly_roc(exp_2group, group_col = "nonexistent")), "not found in sample information")
  
  # Test with invalid pos_class
  expect_error(suppressMessages(gly_roc(exp_2group, pos_class = "invalid")), "not found in group levels")
})

test_that("gly_roc works with different group column names", {
  # Modify sample info to use different group column name
  exp_2group <- exp_2groups() |>
    glyexp::slice_sample_var(n = 5) |>
    glyexp::mutate_obs(condition = group)
  
  result <- suppressMessages(gly_roc(exp_2group, group_col = "condition", pos_class = "H"))
  
  expect_s3_class(result, c("glystats_roc_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("auc", "coords"))
})
