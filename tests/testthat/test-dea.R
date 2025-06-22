test_that("gly_dea works with t-test method", {
  set.seed(123)
  
  # Create test data with clear group differences
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:6),
    group = rep(c("A", "B"), each = 3)
  )
  
  var_info <- tibble::tibble(
    variable = paste0("V", 1:3)
  )
  
  expr_mat <- matrix(c(
    15, 16, 14, 10, 11, 9,   # V1: A higher than B
    12, 13, 11, 8, 9, 7,    # V2: A higher than B  
    20, 21, 19, 15, 16, 14  # V3: A higher than B
  ), nrow = 3, ncol = 6, byrow = TRUE)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Run DEA with t-test
  result <- suppressMessages(gly_dea(exp, method = "t-test"))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_ttest", "glystats_dea_res"))
  expect_equal(nrow(result), 3)
  expect_true("log2fc" %in% colnames(result))  # t-test should have log2fc
  expect_true("p_adj" %in% colnames(result))  # p_adj should exist
  expect_true(all(result$log2fc > 0))  # All should be positive (A > B)
})

test_that("gly_dea works with wilcoxon method", {
  # Create simple test data
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:6),
    group = rep(c("A", "B"), each = 3)
  )
  
  var_info <- tibble::tibble(
    variable = paste0("V", 1:2)
  )
  
  expr_mat <- matrix(c(
    1, 2, 3, 10, 11, 12,   # V1: clear difference
    5, 6, 7, 15, 16, 17    # V2: clear difference
  ), nrow = 2, ncol = 6, byrow = TRUE)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Run DEA with wilcoxon test
  result <- suppressMessages(suppressWarnings(gly_dea(exp, method = "wilcoxon")))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_wilcoxon", "glystats_dea_res"))
  expect_equal(nrow(result), 2)
  expect_false("log2fc" %in% colnames(result))  # Wilcoxon should NOT have log2fc
})

test_that("gly_dea basic functionality works", {
  # Simple test for basic functionality
  sample_info <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    condition = c("X", "Y", "X", "Y")
  )
  
  var_info <- tibble::tibble(variable = c("V1", "V2"))
  
  expr_mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Test all methods work
  expect_no_error(suppressMessages(gly_dea(exp, method = "t-test")))
  expect_no_error(suppressMessages(suppressWarnings(gly_dea(exp, method = "wilcoxon"))))
  expect_no_error(suppressMessages(gly_dea(exp, method = "anova")))
  expect_no_error(suppressMessages(gly_dea(exp, method = "kruskal")))
  
  # Test custom group column
  expect_no_error(suppressMessages(gly_dea(exp, method = "t-test", group_col = "condition")))
})

test_that("gly_dea error handling", {
  # Create simple test data
  sample_info <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B")
  )
  
  var_info <- tibble::tibble(variable = c("V1", "V2"))
  
  expr_mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Test various error conditions
  expect_error(suppressMessages(gly_dea(exp, method = "invalid")), "Invalid method")
  expect_error(suppressMessages(gly_dea(exp, method = "t-test", group_col = "nonexistent")), 
               "not found in sample information")
})

test_that("gly_dea group validation", {
  # Test with 3 groups - should work for multi-group methods only
  sample_info_3g <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "B", "C", "A")
  )
  
  # Test with 1 group - should fail for all methods
  sample_info_1g <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = rep("A", 4)
  )
  
  var_info <- tibble::tibble(variable = c("V1", "V2"))
  expr_mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4)
  rownames(expr_mat) <- var_info$variable
  
  # Test 3 groups
  colnames(expr_mat) <- sample_info_3g$sample
  exp_3g <- glyexp::experiment(expr_mat, sample_info_3g, var_info, "glycomics", "N")
  
  expect_error(suppressMessages(gly_dea(exp_3g, method = "t-test")), "exactly 2 levels")
  expect_error(suppressMessages(gly_dea(exp_3g, method = "wilcoxon")), "exactly 2 levels")
  expect_no_error(suppressMessages(gly_dea(exp_3g, method = "anova")))
  expect_no_error(suppressMessages(gly_dea(exp_3g, method = "kruskal")))
  
  # Test 1 group
  colnames(expr_mat) <- sample_info_1g$sample
  exp_1g <- glyexp::experiment(expr_mat, sample_info_1g, var_info, "glycomics", "N")
  
  expect_error(suppressMessages(gly_dea(exp_1g, method = "t-test")), "exactly 2 levels")
  expect_error(suppressMessages(gly_dea(exp_1g, method = "wilcoxon")), "exactly 2 levels")
  expect_error(suppressMessages(gly_dea(exp_1g, method = "anova")), "at least 2 levels")
  expect_error(suppressMessages(gly_dea(exp_1g, method = "kruskal")), "at least 2 levels")
})

test_that("gly_dea works with anova method", {
  # Create test data with 3 groups
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:9),
    group = rep(c("A", "B", "C"), each = 3)
  )
  
  var_info <- tibble::tibble(
    variable = paste0("V", 1:2)
  )
  
  expr_mat <- matrix(c(
    10, 11, 12, 15, 16, 17, 20, 21, 22,  # V1: clear group differences
    5, 6, 7, 10, 11, 12, 15, 16, 17      # V2: clear group differences
  ), nrow = 2, ncol = 9, byrow = TRUE)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Run DEA with ANOVA
  result <- suppressMessages(gly_dea(exp, method = "anova"))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_anova", "glystats_dea_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("main_test", "post_hoc"))
  
  # Check main_test only group rows
  expect_equal(nrow(result$main_test), 2)
  expect_true("p_adj" %in% colnames(result$main_test))
})

test_that("gly_dea works with kruskal method", {
  # Create test data with 3 groups
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:9),
    group = rep(c("A", "B", "C"), each = 3)
  )
  
  var_info <- tibble::tibble(
    variable = paste0("V", 1:2)
  )
  
  expr_mat <- matrix(c(
    1, 2, 3, 10, 11, 12, 20, 21, 22,  # V1: clear group differences
    5, 6, 7, 15, 16, 17, 25, 26, 27   # V2: clear group differences
  ), nrow = 2, ncol = 9, byrow = TRUE)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Run DEA with Kruskal-Wallis test
  result <- suppressMessages(gly_dea(exp, method = "kruskal"))
  
  # Test core functionality
  expect_s3_class(result, c("glystats_dea_res_kruskal", "glystats_dea_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("main_test", "post_hoc"))
  
  # Check main_test structure
  expect_equal(nrow(result$main_test), 2)
  expect_true("method" %in% colnames(result$main_test))
  expect_true(all(result$main_test$df == 2))  # 3 groups - 1
  expect_false("log2fc" %in% colnames(result$main_test))
})

test_that("gly_dea works with real data", {
  result <- gly_dea(test_gp_exp)
  expect_s3_class(result, c("glystats_dea_res_anova", "glystats_dea_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("main_test", "post_hoc"))
})
