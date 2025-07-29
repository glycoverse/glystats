# Test that all gly_xxx_() functions accept character groups

test_that("gly_xxx_() functions accept character groups", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  
  # Two-group data
  groups_char_2 <- rep(c("A", "B"), each = 5)
  groups_factor_2 <- factor(groups_char_2)
  
  # Test gly_fold_change_
  suppressMessages({
    result_fc_char <- gly_fold_change_(expr_mat, groups_char_2)
    result_fc_factor <- gly_fold_change_(expr_mat, groups_factor_2)
  })
  expect_equal(result_fc_char, result_fc_factor)
  
  # Test gly_ttest_
  suppressMessages({
    result_ttest_char <- gly_ttest_(expr_mat, groups_char_2)
    result_ttest_factor <- gly_ttest_(expr_mat, groups_factor_2)
  })
  expect_equal(result_ttest_char, result_ttest_factor)
  
  # Test gly_wilcox_
  suppressMessages({
    result_wilcox_char <- gly_wilcox_(expr_mat, groups_char_2)
    result_wilcox_factor <- gly_wilcox_(expr_mat, groups_factor_2)
  })
  expect_equal(result_wilcox_char, result_wilcox_factor)
  
  # Test gly_roc_ (requires pROC package)
  skip_if_not_installed("pROC")
  suppressMessages({
    result_roc_char <- gly_roc_(expr_mat, groups_char_2)
    result_roc_factor <- gly_roc_(expr_mat, groups_factor_2)
  })
  expect_equal(result_roc_char, result_roc_factor)
})

test_that("multi-group functions accept character groups", {
  # Create test data for multi-group analysis
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(150)) + 1, nrow = 10, ncol = 15)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:15)
  
  # Three-group data
  groups_char_3 <- rep(c("A", "B", "C"), each = 5)
  groups_factor_3 <- factor(groups_char_3)
  
  # Test gly_anova_
  suppressMessages({
    result_anova_char <- gly_anova_(expr_mat, groups_char_3)
    result_anova_factor <- gly_anova_(expr_mat, groups_factor_3)
  })
  # Results should have same structure even if not identical due to factor level ordering
  expect_equal(names(result_anova_char), names(result_anova_factor))
  expect_equal(colnames(result_anova_char$tidy_result$main_test), 
               colnames(result_anova_factor$tidy_result$main_test))
  
  # Test gly_kruskal_ (requires FSA package)
  skip_if_not_installed("FSA")
  suppressMessages({
    result_kruskal_char <- gly_kruskal_(expr_mat, groups_char_3)
    result_kruskal_factor <- gly_kruskal_(expr_mat, groups_factor_3)
  })
  expect_equal(names(result_kruskal_char), names(result_kruskal_factor))
  expect_equal(colnames(result_kruskal_char$tidy_result$main_test), 
               colnames(result_kruskal_factor$tidy_result$main_test))
  
  # Test gly_limma_ (requires limma package)
  skip_if_not_installed("limma")
  suppressMessages({
    result_limma_char <- gly_limma_(expr_mat, groups_char_3)
    result_limma_factor <- gly_limma_(expr_mat, groups_factor_3)
  })
  expect_equal(names(result_limma_char), names(result_limma_factor))
  expect_equal(colnames(result_limma_char$tidy_result), 
               colnames(result_limma_factor$tidy_result))
  
  # Test gly_plsda_ (requires ropls package)
  skip_if_not_installed("ropls")
  suppressMessages({
    result_plsda_char <- gly_plsda_(expr_mat, groups_char_3)
    result_plsda_factor <- gly_plsda_(expr_mat, groups_factor_3)
  })
  expect_equal(names(result_plsda_char), names(result_plsda_factor))
  expect_equal(names(result_plsda_char$tidy_result), 
               names(result_plsda_factor$tidy_result))
})

test_that("binary classification functions accept character groups", {
  # Create test data for binary classification
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  
  # Two-group data
  groups_char_2 <- rep(c("A", "B"), each = 5)
  groups_factor_2 <- factor(groups_char_2)
  
  # Test gly_oplsda_ (requires ropls package)
  skip_if_not_installed("ropls")
  suppressMessages({
    result_oplsda_char <- gly_oplsda_(expr_mat, groups_char_2)
    result_oplsda_factor <- gly_oplsda_(expr_mat, groups_factor_2)
  })
  expect_equal(names(result_oplsda_char), names(result_oplsda_factor))
  expect_equal(names(result_oplsda_char$tidy_result), 
               names(result_oplsda_factor$tidy_result))
})

test_that("invalid groups input throws error", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  
  # Test with numeric groups (should fail)
  groups_numeric <- rep(c(1, 2), each = 5)
  
  expect_error(gly_fold_change_(expr_mat, groups_numeric), 
               "groups must be a factor or character vector")
  expect_error(gly_ttest_(expr_mat, groups_numeric), 
               "groups must be a factor or character vector")
  expect_error(gly_anova_(expr_mat, groups_numeric), 
               "groups must be a factor or character vector")
})
