test_that("gly_anova works with anova method", {
  # Use test_gp_exp with 3 groups for ANOVA
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)

  # Run DEA with ANOVA
  result <- suppressMessages(gly_anova(exp_3group))

  # Test core functionality - should return a list with two tibbles and S3 class
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_dea_res_anova", "glystats_dea_res", "glystats_res"))
  expect_named(result, c("main_test", "post_hoc_test"))

  # Test main_test tibble
  main_test <- result$main_test
  expect_true(tibble::is_tibble(main_test))
  expect_equal(nrow(main_test), 10)
  expect_true("p_adj" %in% colnames(main_test))
  expect_true("post_hoc" %in% colnames(main_test))

  # Test post_hoc_test tibble
  post_hoc_test <- result$post_hoc_test
  expect_true(tibble::is_tibble(post_hoc_test))
  expect_true(all(c("variable", "group1", "group2", "p_value", "p_adj") %in% colnames(post_hoc_test)))
})

test_that("gly_kruskal works with kruskal method", {
  # Use test_gp_exp with 3 groups for Kruskal-Wallis test
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)

  # Run DEA with Kruskal-Wallis test
  result <- suppressMessages(gly_kruskal(exp_3group))

  # Test core functionality - should return a list with two tibbles and S3 class
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_dea_res_kruskal", "glystats_dea_res", "glystats_res"))
  expect_named(result, c("main_test", "post_hoc_test"))

  # Test main_test tibble
  main_test <- result$main_test
  expect_true(tibble::is_tibble(main_test))
  expect_equal(nrow(main_test), 10)
  expect_true("method" %in% colnames(main_test))
  expect_true(all(main_test$df == 2))  # 3 groups - 1
  expect_true("post_hoc" %in% colnames(main_test))
  expect_false("log2fc" %in% colnames(main_test))

  # Test post_hoc_test tibble
  post_hoc_test <- result$post_hoc_test
  expect_true(tibble::is_tibble(post_hoc_test))
  expect_true(all(c("variable", "group1", "group2", "p_value", "p_adj") %in% colnames(post_hoc_test)))
})

test_that("gly_anova and gly_kruskal basic functionality works", {
  # Test multi-group methods work with test_gp_exp
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)  # Use very small subset
  
  # Multi-group methods
  expect_no_error(suppressMessages(gly_anova(exp_small)))
  expect_no_error(suppressMessages(gly_kruskal(exp_small)))
})

test_that("gly_anova and gly_kruskal error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)
  
  # Test various error conditions - group column not found
  expect_error(suppressMessages(gly_anova(exp_small, group_col = "nonexistent")),
               "not found in sample information")
  expect_error(suppressMessages(gly_kruskal(exp_small, group_col = "nonexistent")),
               "not found in sample information")
})

test_that("gly_anova and gly_kruskal group validation", {
  # Test with 3 groups (using C, H, M from test_gp_exp)
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)
  
  # Test with 1 group 
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)
  
  # Test 3 groups with multi-group methods (should work)
  expect_no_error(suppressMessages(gly_anova(exp_3group)))
  expect_no_error(suppressMessages(gly_kruskal(exp_3group)))
  
  # Test 1 group with multi-group methods
  expect_error(suppressMessages(gly_anova(exp_1group)), "at least 2 levels")
  expect_error(suppressMessages(gly_kruskal(exp_1group)), "at least 2 levels")
})

test_that("gly_anova works with real data", {
  # This test uses the full test_gp_exp to ensure integration works
  result <- suppressMessages(gly_anova(test_gp_exp))
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_dea_res_anova", "glystats_dea_res", "glystats_res"))
  expect_named(result, c("main_test", "post_hoc_test"))
  expect_true(tibble::is_tibble(result$main_test))
  expect_true("post_hoc" %in% colnames(result$main_test))
  expect_true(tibble::is_tibble(result$post_hoc_test))
})

test_that("gly_kruskal works with real data", {
  # This test uses the full test_gp_exp to ensure integration works
  result <- suppressMessages(gly_kruskal(test_gp_exp))
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_dea_res_kruskal", "glystats_dea_res", "glystats_res"))
  expect_named(result, c("main_test", "post_hoc_test"))
  expect_true(tibble::is_tibble(result$main_test))
  expect_true("post_hoc" %in% colnames(result$main_test))
  expect_true(tibble::is_tibble(result$post_hoc_test))
})
