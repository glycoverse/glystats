test_that("gly_anova works with anova method", {
  # Use test_gp_exp with 3 groups for ANOVA
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)

  # Run DEA with ANOVA
  result <- suppressMessages(gly_anova(exp_3group))

  # Test core functionality - should return a list with tidy_result and raw_result
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_anova_res", "glystats_res"))
  expect_named(result, c("tidy_result", "raw_result"))

  # Test tidy_result structure
  expect_type(result$tidy_result, "list")
  expect_named(result$tidy_result, c("main_test", "post_hoc_test"))

  # Test main_test tibble
  main_test <- result$tidy_result$main_test
  expect_true(tibble::is_tibble(main_test))
  expect_equal(nrow(main_test), 10)
  expect_true("p_adj" %in% colnames(main_test))
  expect_true("post_hoc" %in% colnames(main_test))

  # Test post_hoc_test tibble
  post_hoc_test <- result$tidy_result$post_hoc_test
  expect_true(tibble::is_tibble(post_hoc_test))
  expect_true(all(c("variable", "ref_group", "test_group", "p_val", "p_adj", "log2fc") %in% colnames(post_hoc_test)))

  # Test raw_result structure
  expect_type(result$raw_result, "list")
  expect_named(result$raw_result, c("main_test", "post_hoc_test"))
})

test_that("gly_anova comparison direction is correct for 2 groups", {
  # Create a test experiment with 2 groups
  # group B has higher mean than group A
  var_info <- tibble::tibble(variable = "V1")
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:20),
    group = factor(rep(c("A", "B"), each = 10), levels = c("A", "B"))
  )
  expr_mat <- matrix(
    c(rnorm(10, mean = 1, sd = 0.1), rnorm(10, mean = 2, sd = 0.1)),
    nrow = 1, byrow = TRUE
  )
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "glycomics", "N")

  # Call gly_anova
  result <- suppressMessages(gly_anova(exp))

  # Test post_hoc
  expect_identical(result$tidy_result$main_test$post_hoc, "A_vs_B")
  expect_true(all(result$tidy_result$post_hoc_test$log2fc > 0))
  expect_identical(result$tidy_result$post_hoc_test$ref_group, "A")
  expect_identical(result$tidy_result$post_hoc_test$test_group, "B")
})

test_that("gly_anova comparison direction is correct for 3 groups", {
  # Create a test experiment with 3 groups
  # Group means: A < B < C
  var_info <- tibble::tibble(variable = "V1")
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:30),
    group = factor(rep(c("A", "B", "C"), each = 10), levels = c("A", "B", "C"))
  )
  expr_mat <- matrix(
    c(rnorm(10, mean = 1, sd = 0.1), rnorm(10, mean = 2, sd = 0.1), rnorm(10, mean = 3, sd = 0.1)),
    nrow = 1, byrow = TRUE
  )
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "glycomics", "N")

  # Call gly_anova
  result <- suppressMessages(gly_anova(exp))

  # Test post_hoc
  expect_setequal(
    stringr::str_split_1(result$tidy_result$main_test$post_hoc, ";"),
    c("A_vs_B", "A_vs_C", "B_vs_C")
  )
  comparisons <- stringr::str_c(
    result$tidy_result$post_hoc_test$ref_group,
    "_vs_",
    result$tidy_result$post_hoc_test$test_group
  )
  expect_setequal(comparisons, c("A_vs_B", "A_vs_C", "B_vs_C"))
  expect_true(all(result$tidy_result$post_hoc_test$log2fc > 0))
})

test_that("gly_kruskal comparison direction is correct for 2 groups", {
  # Create a test experiment with 2 groups
  # group B has higher mean than group A
  var_info <- tibble::tibble(variable = "V1")
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:20),
    group = factor(rep(c("A", "B"), each = 10), levels = c("A", "B"))
  )
  expr_mat <- matrix(1:20, nrow = 1, byrow = TRUE)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "glycomics", "N")

  # Call gly_kruskal
  result <- suppressMessages(gly_kruskal(exp))

  # Test post_hoc
  expect_identical(result$tidy_result$main_test$post_hoc, "A_vs_B")
  expect_true(result$tidy_result$post_hoc_test$log2fc > 0)
  expect_identical(result$tidy_result$post_hoc_test$ref_group, "A")
  expect_identical(result$tidy_result$post_hoc_test$test_group, "B")
})

test_that("gly_kruskal comparison direction is correct for 3 groups", {
  # Create a test experiment with 3 groups
  # Group means: A < B < C
  var_info <- tibble::tibble(variable = "V1")
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:30),
    group = factor(rep(c("A", "B", "C"), each = 10), levels = c("A", "B", "C"))
  )
  expr_mat <- matrix(1:30, nrow = 1, byrow = TRUE)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "glycomics", "N")

  # Call gly_kruskal
  result <- suppressMessages(gly_kruskal(exp))

  # Test post_hoc
  expect_setequal(
    stringr::str_split_1(result$tidy_result$main_test$post_hoc, ";"),
    c("A_vs_B", "A_vs_C", "B_vs_C")
  )
  comparisons <- stringr::str_c(
    result$tidy_result$post_hoc_test$ref_group,
    "_vs_",
    result$tidy_result$post_hoc_test$test_group
  )
  expect_setequal(comparisons, c("A_vs_B", "A_vs_C", "B_vs_C"))
  expect_true(all(result$tidy_result$post_hoc_test$log2fc > 0))
})

test_that("gly_anova assigns NA for failed variables", {
  # Use test_gp_exp with 3 groups for ANOVA
  # The first three variables are set to NA
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)
  exp_3group$expr_mat[1:3, ] <- NA  # This will lead to stats::aov() failing
  na_vars <- exp_3group$var_info$variable[1:3]

  # Run DEA with ANOVA
  expect_warning(result <- suppressMessages(gly_anova(exp_3group)))

  # Test main test
  main_test_raw <- result$raw_result$main_test
  expect_true(all(is.na(main_test_raw[na_vars])))
  main_test_tidy <- result$tidy_result$main_test
  p_values <- main_test_tidy |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(p_val)
  expect_true(all(is.na(p_values)))
})

test_that("gly_kruskal works with kruskal method", {
  # Use test_gp_exp with 3 groups for Kruskal-Wallis test
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)

  # Run DEA with Kruskal-Wallis test
  result <- suppressMessages(gly_kruskal(exp_3group))

  # Test core functionality - should return a list with tidy_result and raw_result
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_kruskal_res", "glystats_res"))
  expect_named(result, c("tidy_result", "raw_result"))

  # Test tidy_result structure
  expect_type(result$tidy_result, "list")
  expect_named(result$tidy_result, c("main_test", "post_hoc_test"))

  # Test main_test tibble
  main_test <- result$tidy_result$main_test
  expect_true(tibble::is_tibble(main_test))
  expect_equal(nrow(main_test), 10)
  expect_true("method" %in% colnames(main_test))

  expect_true("post_hoc" %in% colnames(main_test))
  expect_false("log2fc" %in% colnames(main_test))

  # Test post_hoc_test tibble
  post_hoc_test <- result$tidy_result$post_hoc_test
  expect_true(tibble::is_tibble(post_hoc_test))
  expect_true(all(c("variable", "ref_group", "test_group", "p_val", "p_adj", "log2fc") %in% colnames(post_hoc_test)))

  # Test raw_result structure
  expect_type(result$raw_result, "list")
  expect_named(result$raw_result, c("main_test", "post_hoc_test"))
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

test_that("gly_anova_ works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B", "C"), c(3, 3, 4)))

  # Test function execution
  suppressMessages({
    result <- gly_anova_(expr_mat, groups)
  })

  # Verify results
  expect_s3_class(result, "glystats_anova_res")
  expect_type(result, "list")
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result$main_test))
  expect_equal(nrow(result$tidy_result$main_test), 10)
})

test_that("gly_kruskal_ works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B", "C"), c(3, 3, 4)))

  # Test function execution
  suppressMessages({
    result <- gly_kruskal_(expr_mat, groups)
  })

  # Verify results
  expect_s3_class(result, "glystats_kruskal_res")
  expect_type(result, "list")
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result$main_test))
  expect_equal(nrow(result$tidy_result$main_test), 10)
})
