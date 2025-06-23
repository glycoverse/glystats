test_that("gly_dea works with t-test method", {
  # Use test_gp_exp and filter to 2 groups for t-test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing
  
  # Run DEA with t-test
  result <- suppressMessages(gly_dea(exp_2group, method = "t-test"))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_ttest", "glystats_dea_res", "glystats_res"))
  expect_equal(nrow(result), 10)
  expect_true("log2fc" %in% colnames(result))  # t-test should have log2fc
  expect_true("p_adj" %in% colnames(result))  # p_adj should exist
})

test_that("gly_dea works with wilcoxon method", {
  # Use test_gp_exp and filter to 2 groups for wilcoxon test
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("M", "Y")) |>
    glyexp::slice_sample_var(n = 10)  # Use smaller subset for faster testing
  
  # Run DEA with wilcoxon test
  result <- suppressMessages(suppressWarnings(gly_dea(exp_2group, method = "wilcoxon")))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_wilcoxon", "glystats_dea_res", "glystats_res"))
  expect_equal(nrow(result), 10)
  expect_false("log2fc" %in% colnames(result))  # Wilcoxon should NOT have log2fc
})

test_that("gly_dea basic functionality works", {
  # Test all methods work with test_gp_exp
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)  # Use very small subset
  
  # 2-group methods
  exp_2group <- exp_small |> glyexp::filter_obs(group %in% c("C", "H"))
  expect_no_error(suppressMessages(gly_dea(exp_2group, method = "t-test")))
  expect_no_error(suppressMessages(suppressWarnings(gly_dea(exp_2group, method = "wilcoxon"))))
  
  # Multi-group methods
  expect_no_error(suppressMessages(gly_dea(exp_small, method = "anova")))
  expect_no_error(suppressMessages(gly_dea(exp_small, method = "kruskal")))
})

test_that("gly_dea error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)
  
  # Test various error conditions
  expect_error(suppressMessages(gly_dea(exp_small, method = "invalid")), "Invalid method")
  expect_error(suppressMessages(gly_dea(exp_small, method = "t-test", group_col = "nonexistent")), 
               "not found in sample information")
})

test_that("gly_dea group validation", {
  # Test with 3 groups (using C, H, M from test_gp_exp)
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test with 1 group 
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)
  
  # Test 3 groups
  expect_error(suppressMessages(gly_dea(exp_3group, method = "t-test")), "exactly 2 levels")
  expect_error(suppressMessages(gly_dea(exp_3group, method = "wilcoxon")), "exactly 2 levels")
  expect_no_error(suppressMessages(gly_dea(exp_3group, method = "anova")))
  expect_no_error(suppressMessages(gly_dea(exp_3group, method = "kruskal")))
  
  # Test 1 group
  expect_error(suppressMessages(gly_dea(exp_1group, method = "t-test")), "exactly 2 levels")
  expect_error(suppressMessages(gly_dea(exp_1group, method = "wilcoxon")), "exactly 2 levels")
  expect_error(suppressMessages(gly_dea(exp_1group, method = "anova")), "at least 2 levels")
  expect_error(suppressMessages(gly_dea(exp_1group, method = "kruskal")), "at least 2 levels")
})

test_that("gly_dea works with anova method", {
  # Use test_gp_exp with 3 groups for ANOVA
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)
  
  # Run DEA with ANOVA
  result <- suppressMessages(gly_dea(exp_3group, method = "anova"))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_anova", "glystats_dea_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("main_test", "post_hoc"))
  
  # Check main_test structure
  expect_equal(nrow(result$main_test), 10)
  expect_true("p_adj" %in% colnames(result$main_test))
})

test_that("gly_dea works with kruskal method", {
  # Use test_gp_exp with 3 groups for Kruskal-Wallis test
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)
  
  # Run DEA with Kruskal-Wallis test
  result <- suppressMessages(gly_dea(exp_3group, method = "kruskal"))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_kruskal", "glystats_dea_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("main_test", "post_hoc"))
  
  # Check main_test structure
  expect_equal(nrow(result$main_test), 10)
  expect_true("method" %in% colnames(result$main_test))
  expect_true(all(result$main_test$df == 2))  # 3 groups - 1
  expect_false("log2fc" %in% colnames(result$main_test))
})

test_that("gly_dea works with real data", {
  # This test uses the full test_gp_exp to ensure integration works
  result <- gly_dea(test_gp_exp)
  expect_s3_class(result, c("glystats_dea_res_anova", "glystats_dea_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("main_test", "post_hoc"))
})
