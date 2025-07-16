test_that("add_info parameter works correctly for functions returning tibbles with variable column", {
  # Use test_gp_exp and filter to 2 groups
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test gly_fold_change
  result_no_info <- suppressMessages(gly_fold_change(exp_2group, add_info = FALSE))
  result_with_info <- suppressMessages(gly_fold_change(exp_2group, add_info = TRUE))
  
  expect_equal(ncol(result_no_info), 2)
  expect_true(ncol(result_with_info) > 2)
  expect_true("variable" %in% colnames(result_with_info))
  expect_true("log2fc" %in% colnames(result_with_info))
  
  # Test gly_ttest
  result_ttest_no_info <- suppressMessages(gly_ttest(exp_2group, add_info = FALSE))
  result_ttest_with_info <- suppressMessages(gly_ttest(exp_2group, add_info = TRUE))
  
  expect_true(ncol(result_ttest_no_info) < ncol(result_ttest_with_info))
  expect_true("variable" %in% colnames(result_ttest_with_info))
})

test_that("add_info parameter works correctly for functions returning tibbles with sample column", {
  # Use test_gp_exp
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10)
  
  # Test gly_pca
  result_pca_no_info <- suppressMessages(gly_pca(exp_subset, add_info = FALSE))
  result_pca_with_info <- suppressMessages(gly_pca(exp_subset, add_info = TRUE))
  
  # Check samples tibble
  expect_true(ncol(result_pca_no_info$samples) < ncol(result_pca_with_info$samples))
  expect_true("sample" %in% colnames(result_pca_with_info$samples))
  
  # Check variables tibble
  expect_true(ncol(result_pca_no_info$variables) < ncol(result_pca_with_info$variables))
  expect_true("variable" %in% colnames(result_pca_with_info$variables))
  
  # Check eigenvalues tibble (should be the same since it doesn't have variable/sample columns)
  expect_equal(ncol(result_pca_no_info$eigenvalues), ncol(result_pca_with_info$eigenvalues))
})

test_that("add_info parameter works correctly for ROC analysis", {
  # Use test_gp_exp and filter to 2 groups
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test gly_roc
  result_roc_no_info <- suppressMessages(gly_roc(exp_2group, add_info = FALSE))
  result_roc_with_info <- suppressMessages(gly_roc(exp_2group, add_info = TRUE))
  
  # Check auc tibble
  expect_true(ncol(result_roc_no_info$auc) < ncol(result_roc_with_info$auc))
  expect_true("variable" %in% colnames(result_roc_with_info$auc))
  
  # Check coords tibble
  expect_true(ncol(result_roc_no_info$coords) < ncol(result_roc_with_info$coords))
  expect_true("variable" %in% colnames(result_roc_with_info$coords))
})
