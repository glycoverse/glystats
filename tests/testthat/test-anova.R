test_that("gly_anova works with anova method", {
  # Use test_gp_exp with 3 groups for ANOVA
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)

  # Run DEA with ANOVA
  result <- suppressMessages(gly_anova(exp_3group))

  # Test core functionality - should return a list with tidy_result, raw_result, and meta_data
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_anova_res", "glystats_res"))
  expect_named(result, c("tidy_result", "raw_result", "meta_data"))

  # Test tidy_result structure
  expect_type(result$tidy_result, "list")
  expect_named(result$tidy_result, c("main_test", "post_hoc_test"))

  # Test main_test tibble
  main_test <- result$tidy_result$main_test
  expect_true(tibble::is_tibble(main_test))
  expect_equal(nrow(main_test), 10)
  expect_true("p_adj" %in% colnames(main_test))
  expect_true("post_hoc" %in% colnames(main_test))
  expect_true("effect_size" %in% colnames(main_test))
  expect_type(main_test$effect_size, "double")

  # Test post_hoc_test tibble
  post_hoc_test <- result$tidy_result$post_hoc_test
  expect_true(tibble::is_tibble(post_hoc_test))
  expect_true(all(
    c("variable", "ref_group", "test_group", "p_val", "p_adj", "log2fc") %in%
      colnames(post_hoc_test)
  ))

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
    nrow = 1,
    byrow = TRUE
  )
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

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
    c(
      rnorm(10, mean = 1, sd = 0.1),
      rnorm(10, mean = 2, sd = 0.1),
      rnorm(10, mean = 3, sd = 0.1)
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

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
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

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
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

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

test_that("gly_anova raw post-hoc uses test-minus-reference direction", {
  var_info <- tibble::tibble(variable = "V1")
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:30),
    group = factor(rep(c("A", "B", "C"), each = 10), levels = c("A", "B", "C"))
  )
  expr_mat <- matrix(1:30, nrow = 1, byrow = TRUE)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

  result <- suppressMessages(gly_anova(exp))
  raw_post_hoc <- result$raw_result$post_hoc_test[[1]]$group

  expect_setequal(rownames(raw_post_hoc), c("B-A", "C-A", "C-B"))
  expect_true(all(raw_post_hoc[, "diff"] > 0))
  expect_true(all(raw_post_hoc[, "lwr"] > 0))
  expect_true(all(raw_post_hoc[, "upr"] > 0))
})

test_that("gly_kruskal raw post-hoc uses test-minus-reference direction", {
  var_info <- tibble::tibble(variable = "V1")

  sample_info_2 <- tibble::tibble(
    sample = paste0("S", 1:20),
    group = factor(rep(c("A", "B"), each = 10), levels = c("A", "B"))
  )
  expr_mat_2 <- matrix(1:20, nrow = 1, byrow = TRUE)
  colnames(expr_mat_2) <- sample_info_2$sample
  rownames(expr_mat_2) <- var_info$variable
  exp_2 <- glyexp::experiment(expr_mat_2, sample_info_2, var_info, "others")

  result_2 <- suppressMessages(gly_kruskal(exp_2))
  raw_post_hoc_2 <- result_2$raw_result$post_hoc_test[[1]]$res

  expect_identical(raw_post_hoc_2$Comparison, "B - A")

  sample_info_3 <- tibble::tibble(
    sample = paste0("S", 1:30),
    group = factor(rep(c("A", "B", "C"), each = 10), levels = c("A", "B", "C"))
  )
  expr_mat_3 <- matrix(1:30, nrow = 1, byrow = TRUE)
  colnames(expr_mat_3) <- sample_info_3$sample
  rownames(expr_mat_3) <- var_info$variable
  exp_3 <- glyexp::experiment(expr_mat_3, sample_info_3, var_info, "others")

  result_3 <- suppressMessages(gly_kruskal(exp_3))
  raw_post_hoc_3 <- result_3$raw_result$post_hoc_test[[1]]$res

  expect_setequal(raw_post_hoc_3$Comparison, c("B - A", "C - A", "C - B"))
  expect_true(all(raw_post_hoc_3$Z > 0))
})

test_that("gly_anova assigns NA for failed variables", {
  # Use test_gp_exp with 3 groups for ANOVA
  # The first three variables are set to NA
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10)
  exp_3group$expr_mat[1:3, ] <- NA # This will lead to stats::aov() failing
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

  # Test core functionality - should return a list with tidy_result, raw_result, and meta_data
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_kruskal_res", "glystats_res"))
  expect_named(result, c("tidy_result", "raw_result", "meta_data"))

  # Test tidy_result structure
  expect_type(result$tidy_result, "list")
  expect_named(result$tidy_result, c("main_test", "post_hoc_test"))

  # Test main_test tibble
  main_test <- result$tidy_result$main_test
  expect_true(tibble::is_tibble(main_test))
  expect_equal(nrow(main_test), 10)
  expect_true("method" %in% colnames(main_test))

  expect_true("post_hoc" %in% colnames(main_test))
  expect_true("effect_size" %in% colnames(main_test))
  expect_type(main_test$effect_size, "double")
  expect_false("log2fc" %in% colnames(main_test))

  # Test post_hoc_test tibble
  post_hoc_test <- result$tidy_result$post_hoc_test
  expect_true(tibble::is_tibble(post_hoc_test))
  expect_true(all(
    c("variable", "ref_group", "test_group", "p_val", "p_adj", "log2fc") %in%
      colnames(post_hoc_test)
  ))

  # Test raw_result structure
  expect_type(result$raw_result, "list")
  expect_named(result$raw_result, c("main_test", "post_hoc_test"))
})

test_that("gly_anova and gly_kruskal basic functionality works", {
  # Test multi-group methods work with test_gp_exp
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5) # Use very small subset

  # Multi-group methods
  expect_no_error(suppressMessages(gly_anova(exp_small)))
  expect_no_error(suppressMessages(gly_kruskal(exp_small)))
})

test_that("gly_anova and gly_kruskal error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)

  # Test various error conditions - group column not found
  expect_error(
    suppressMessages(gly_anova(exp_small, group_col = "nonexistent")),
    "not found in sample information"
  )
  expect_error(
    suppressMessages(gly_kruskal(exp_small, group_col = "nonexistent")),
    "not found in sample information"
  )
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
  expect_true("effect_size" %in% colnames(result$tidy_result$main_test))
})

test_that("gly_anova_ returns eta-squared in effect_size", {
  expr_mat <- matrix(
    c(1, 2, 3, 5, 6, 7, 9, 10, 11),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:9))
  )
  groups <- factor(rep(c("A", "B", "C"), each = 3))

  result <- suppressMessages(gly_anova_(
    expr_mat,
    groups,
    p_adj_method = NULL
  ))

  log_values <- log2(c(1, 2, 3, 5, 6, 7, 9, 10, 11) + 1)
  grand_mean <- mean(log_values)
  group_means <- purrr::map_dbl(levels(groups), function(group) {
    mean(log_values[groups == group])
  })
  group_sizes <- purrr::map_int(levels(groups), function(group) {
    sum(groups == group)
  })
  ss_between <- sum(group_sizes * (group_means - grand_mean)^2)
  ss_total <- sum((log_values - grand_mean)^2)
  expected <- ss_between / ss_total

  expect_equal(
    result$tidy_result$main_test$effect_size,
    expected,
    tolerance = 1e-10
  )
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
  expect_true("effect_size" %in% colnames(result$tidy_result$main_test))
})

test_that("gly_kruskal_ returns epsilon-squared in effect_size", {
  expr_mat <- matrix(
    c(1, 2, 3, 5, 6, 7, 9, 10, 11),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:9))
  )
  groups <- factor(rep(c("A", "B", "C"), each = 3))

  result <- suppressMessages(gly_kruskal_(
    expr_mat,
    groups,
    p_adj_method = NULL
  ))

  log_values <- log2(c(1, 2, 3, 5, 6, 7, 9, 10, 11) + 1)
  h_stat <- unname(stats::kruskal.test(log_values ~ groups)$statistic)
  n_obs <- length(log_values)
  n_groups <- nlevels(groups)
  expected <- (h_stat - n_groups + 1) / (n_obs - n_groups)

  expect_equal(
    result$tidy_result$main_test$effect_size,
    expected,
    tolerance = 1e-10
  )
})

test_that(".tibblify_main_test_results validates effect-size inputs", {
  expr_mat <- matrix(
    c(1, 2, 3, 5, 6, 7, 9, 10, 11),
    nrow = 1,
    dimnames = list("var1", paste0("sample", 1:9))
  )
  groups <- factor(rep(c("A", "B", "C"), each = 3))
  log_values <- log2(expr_mat["var1", ] + 1)
  main_test_raw <- list(var1 = stats::kruskal.test(log_values ~ groups))

  expect_error(
    .tibblify_main_test_results(
      main_test_raw,
      stats::kruskal.test,
      p_adj_method = NULL,
      effect_size_method = "epsilon_squared"
    ),
    "must be supplied"
  )

  expect_error(
    .add_effect_size_to_main_test(
      tibble::tibble(variable = "var1", statistic = unname(main_test_raw$var1$statistic)),
      effect_size_method = "bad_method",
      expr_mat = expr_mat,
      groups = groups
    ),
    "Must be element of set"
  )
})

test_that("post_hoc_test should not contain NA rows when add_info is TRUE", {
  # Run ANOVA with add_info = TRUE (default)
  result <- suppressMessages(gly_anova(test_gp_exp))

  post_hoc_test <- result$tidy_result$post_hoc_test

  # Count total NA values in ref_group
  n_na_ref_group <- sum(is.na(post_hoc_test$ref_group))

  # There should be no NA values in ref_group
  # because post_hoc_test should only contain variables with significant main effects
  expect_equal(
    n_na_ref_group,
    0,
    info = "post_hoc_test should not have NA values in ref_group column"
  )

  # Verify that all rows have complete data
  expect_equal(sum(is.na(post_hoc_test$test_group)), 0)
  expect_equal(sum(is.na(post_hoc_test$p_val)), 0)
  expect_equal(sum(is.na(post_hoc_test$p_adj)), 0)

  # Verify that post_hoc_test only contains variables from main_test that are significant
  main_test <- result$tidy_result$main_test
  sig_vars <- main_test$variable[main_test$p_adj < 0.05]

  # All variables in post_hoc_test should be significant in main_test
  expect_true(
    all(post_hoc_test$variable %in% sig_vars),
    info = "All variables in post_hoc_test should have significant main effects"
  )
})
