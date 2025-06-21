test_that("gly_dea works with t-test method", {
  # Set seed for reproducible results
  set.seed(123)
  
  # Create test data
  n_vars <- 5
  n_samples_per_group <- 6
  total_samples <- n_samples_per_group * 2
  
  # Create sample information with two groups
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:total_samples),
    group = rep(c("A", "B"), each = n_samples_per_group),
    batch = rep(1:2, times = n_samples_per_group)
  )
  
  # Create variable information
  var_info <- tibble::tibble(
    variable = paste0("V", 1:n_vars),
    protein = paste0("PRO", rep(1:3, length.out = n_vars)),
    glycan_composition = paste0("H", 3:7, "N2")
  )
  
  # Create expression matrix with some difference between groups
  # Group A: mean around 10, Group B: mean around 15
  expr_mat <- matrix(nrow = n_vars, ncol = total_samples)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  for (i in 1:n_vars) {
    # Group A samples (higher expression)
    expr_mat[i, 1:n_samples_per_group] <- rnorm(n_samples_per_group, mean = 15, sd = 2)
    # Group B samples (lower expression)  
    expr_mat[i, (n_samples_per_group + 1):total_samples] <- rnorm(n_samples_per_group, mean = 10, sd = 2)
  }
  
  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Run DEA with t-test
  result <- suppressMessages(gly_dea(exp, method = "t-test", group_col = "group"))
  
  # Check result structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), n_vars)
  
  # Check required columns exist
  expect_true("variable" %in% colnames(result))
  expect_true("p" %in% colnames(result))
  expect_true("log2fc" %in% colnames(result))
  expect_true("mean_group1" %in% colnames(result))
  expect_true("mean_group2" %in% colnames(result))
  
  # Check variable names match
  expect_equal(sort(result$variable), sort(var_info$variable))
  
  # Check p-values are numeric and in valid range
  expect_true(all(is.numeric(result$p)))
  expect_true(all(result$p >= 0 & result$p <= 1))
  
  # Check log2fc is numeric
  expect_true(all(is.numeric(result$log2fc)))
  
  # Since we set group A to have higher expression than group B,
  # most log2fc values should be positive (mean_group1 - mean_group2)
  expect_true(mean(result$log2fc > 0) > 0.5)
})

test_that("gly_dea works with wilcoxon method", {
  # Set seed for reproducible results
  set.seed(456)
  
  # Create test data
  n_vars <- 4
  n_samples_per_group <- 5
  total_samples <- n_samples_per_group * 2
  
  # Create sample information with two groups
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:total_samples),
    group = rep(c("Control", "Treatment"), each = n_samples_per_group)
  )
  
  # Create variable information
  var_info <- tibble::tibble(
    variable = paste0("Var", 1:n_vars),
    annotation = paste0("Feature_", 1:n_vars)
  )
  
  # Create expression matrix with clear differences
  expr_mat <- matrix(nrow = n_vars, ncol = total_samples)
  rownames(expr_mat) <- var_info$variable
  colnames(expr_mat) <- sample_info$sample
  
  for (i in 1:n_vars) {
    # Control group (lower values)
    expr_mat[i, 1:n_samples_per_group] <- rpois(n_samples_per_group, lambda = 5)
    # Treatment group (higher values)
    expr_mat[i, (n_samples_per_group + 1):total_samples] <- rpois(n_samples_per_group, lambda = 12)
  }
  
  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Run DEA with wilcoxon test
  result <- suppressMessages(suppressWarnings(gly_dea(exp, method = "wilcoxon", group_col = "group")))
  
  # Check result structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), n_vars)
  
  # Check required columns exist
  expect_true("variable" %in% colnames(result))
  expect_true("p" %in% colnames(result))
  
  # Wilcoxon test should NOT have log2fc column
  expect_false("log2fc" %in% colnames(result))
  
  # Check variable names match
  expect_equal(sort(result$variable), sort(var_info$variable))
  
  # Check p-values are numeric and in valid range
  expect_true(all(is.numeric(result$p)))
  expect_true(all(result$p >= 0 & result$p <= 1))
  
  # Given the clear difference we set up, most p-values should be small
  expect_true(mean(result$p < 0.05) > 0.5)
})

test_that("gly_dea works with valid inputs", {
  # Create a simple test experiment
  sample_info <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    condition = c("X", "Y", "X", "Y")
  )
  
  var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    annotation = c("Ann1", "Ann2")
  )
  
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
  
  # Test with valid inputs
  expect_no_error(suppressMessages(gly_dea(exp, method = "t-test")))
  expect_no_error(suppressMessages(suppressWarnings(gly_dea(exp, method = "wilcoxon"))))
  
  # Test with custom group column
  expect_no_error(suppressMessages(gly_dea(exp, method = "t-test", group_col = "condition")))
})

test_that("gly_dea throws error with invalid method", {
  # Create a simple test experiment
  sample_info <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B")
  )
  
  var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    annotation = c("Ann1", "Ann2")
  )
  
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
  
  # Test error: invalid method
  expect_error(suppressMessages(gly_dea(exp, method = "invalid")), "Invalid method")
})

test_that("gly_dea throws error with non-existent group column", {
  # Create a simple test experiment
  sample_info <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B")
  )
  
  var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    annotation = c("Ann1", "Ann2")
  )
  
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
  
  # Test error: non-existent group column
  expect_error(suppressMessages(gly_dea(exp, method = "t-test", group_col = "nonexistent")), 
               "not found in sample information")
})

test_that("gly_dea throws error with more than 2 groups", {
  # Create test experiment with 3 groups
  sample_info <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "B", "C", "A")
  )
  
  var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    annotation = c("Ann1", "Ann2")
  )
  
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
  
  # Test error: more than 2 groups
  expect_error(suppressMessages(gly_dea(exp, method = "t-test")), 
               "must be a factor with exactly 2 levels")
})

test_that("gly_dea throws error with only 1 group", {
  # Create test experiment with only 1 group
  sample_info <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    group = rep("A", 4)
  )
  
  var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    annotation = c("Ann1", "Ann2")
  )
  
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
  
  # Test error: only 1 group
  expect_error(suppressMessages(gly_dea(exp, method = "t-test")),
               "must be a factor with exactly 2 levels")
})
