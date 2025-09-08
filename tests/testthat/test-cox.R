test_that("gly_cox works with basic survival data", {
  # Create test survival data
  set.seed(123)
  
  # Create expression matrix (10 variables, 20 samples)
  expr_mat <- matrix(rnorm(200), nrow = 10, ncol = 20)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:20)
  
  # Create sample info with survival data
  sample_info <- tibble::tibble(
    sample = paste0("sample", 1:20),
    time = rexp(20, rate = 0.1),  # Exponential survival times
    event = rbinom(20, 1, 0.7),  # 70% event rate
    group = rep(c("A", "B"), each = 10)
  )
  
  # Create variable info
  var_info <- tibble::tibble(
    variable = paste0("var", 1:10),
    type = rep("biomarker", 10)
  )
  
  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Test gly_cox function
  result <- suppressMessages(gly_cox(exp))

  # Check result structure
  expect_s3_class(result, c("glystats_cox_res", "glystats_res"))
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check tidy_result structure
  expect_true(tibble::is_tibble(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 10)

  # Check required columns in tidy_result
  expect_true("variable" %in% colnames(result$tidy_result))
  expect_true("coefficient" %in% colnames(result$tidy_result))
  expect_true("hr" %in% colnames(result$tidy_result))
  expect_true("p" %in% colnames(result$tidy_result))
  expect_true("p_adj" %in% colnames(result$tidy_result))

  # Check that add_info worked (variable info should be joined)
  expect_true("type" %in% colnames(result$tidy_result))

  # Check raw_result structure
  expect_type(result$raw_result, "list")
  expect_length(result$raw_result, 10)
  expect_true(all(purrr::map_lgl(result$raw_result, ~ inherits(.x, "coxph"))))
})

test_that("gly_cox assigns NA for failed variables", {
  # Create test survival data
  set.seed(123)

  # Create expression matrix (10 variables, 20 samples)
  expr_mat <- matrix(rnorm(200), nrow = 10, ncol = 20)
  expr_mat[1:3, ] <- NA  # This will lead to stats::coxph() failing
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:20)
  na_vars <- paste0("var", 1:3)

  # Create sample info with survival data
  sample_info <- tibble::tibble(
    sample = paste0("sample", 1:20),
    time = rexp(20, rate = 0.1),  # Exponential survival times
    event = rbinom(20, 1, 0.7),  # 70% event rate
    group = rep(c("A", "B"), each = 10)
  )

  # Create variable info
  var_info <- tibble::tibble(
    variable = paste0("var", 1:10),
    type = rep("biomarker", 10)
  )

  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )

  # Run DEA with cox test
  expect_warning(result <- suppressMessages(gly_cox(exp)))

  # Test results
  expect_true(all(is.na(result$raw_result[na_vars])))
  p_values <- result$tidy_result |>
    dplyr::filter(variable %in% na_vars) |>
    dplyr::pull(p)
  expect_true(all(is.na(p_values)))
})

test_that("gly_cox_ works with matrix input", {
  # Create test data
  set.seed(456)
  expr_mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("gene", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  
  time <- rexp(10, rate = 0.2)
  event <- rbinom(10, 1, 0.6)
  
  # Test function execution
  result <- suppressMessages(gly_cox_(expr_mat, time, event))

  # Verify results
  expect_s3_class(result, c("glystats_cox_res", "glystats_res"))
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check tidy_result
  expect_equal(nrow(result$tidy_result), 10)
  expect_true("coefficient" %in% colnames(result$tidy_result))
  expect_true("hr" %in% colnames(result$tidy_result))
  expect_true("p" %in% colnames(result$tidy_result))
  expect_true("p_adj" %in% colnames(result$tidy_result))

  # Check raw_result
  expect_type(result$raw_result, "list")
  expect_length(result$raw_result, 10)
})

test_that("gly_cox returns list with tidy and raw results", {
  # Create test survival data
  set.seed(789)

  # Create expression matrix (5 variables, 15 samples)
  expr_mat <- matrix(rnorm(75), nrow = 5, ncol = 15)
  rownames(expr_mat) <- paste0("var", 1:5)
  colnames(expr_mat) <- paste0("sample", 1:15)

  # Create sample info with survival data
  sample_info <- tibble::tibble(
    sample = paste0("sample", 1:15),
    time = rexp(15, rate = 0.15),
    event = rbinom(15, 1, 0.8)
  )

  # Create variable info
  var_info <- tibble::tibble(
    variable = paste0("var", 1:5),
    type = rep("marker", 5)
  )

  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )

  # Test gly_cox function
  result <- suppressMessages(gly_cox(exp))

  # Check result structure
  expect_s3_class(result, c("glystats_cox_res", "glystats_res"))
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))

  # Check tidy_result
  expect_true(tibble::is_tibble(result$tidy_result))
  expect_equal(nrow(result$tidy_result), 5)
  expect_true("type" %in% colnames(result$tidy_result))  # add_info worked

  # Check raw_result
  expect_type(result$raw_result, "list")
  expect_length(result$raw_result, 5)
  expect_true(all(purrr::map_lgl(result$raw_result, ~ inherits(.x, "coxph"))))
})

test_that("gly_cox add_info parameter works", {
  # Create test survival data
  set.seed(321)
  
  # Create expression matrix (3 variables, 12 samples)
  expr_mat <- matrix(rnorm(36), nrow = 3, ncol = 12)
  rownames(expr_mat) <- paste0("var", 1:3)
  colnames(expr_mat) <- paste0("sample", 1:12)
  
  # Create sample info with survival data
  sample_info <- tibble::tibble(
    sample = paste0("sample", 1:12),
    time = rexp(12, rate = 0.1),
    event = rbinom(12, 1, 0.75)
  )
  
  # Create variable info
  var_info <- tibble::tibble(
    variable = paste0("var", 1:3),
    category = c("A", "B", "C"),
    importance = c(1, 2, 3)
  )
  
  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )
  
  # Test add_info = TRUE (default)
  result_with_info <- suppressMessages(gly_cox(exp, add_info = TRUE))
  expect_true("category" %in% colnames(result_with_info$tidy_result))
  expect_true("importance" %in% colnames(result_with_info$tidy_result))

  # Test add_info = FALSE
  result_without_info <- suppressMessages(gly_cox(exp, add_info = FALSE))
  expect_false("category" %in% colnames(result_without_info$tidy_result))
  expect_false("importance" %in% colnames(result_without_info$tidy_result))
})

test_that("gly_cox custom time and event columns work", {
  # Create test survival data with custom column names
  set.seed(654)

  # Create expression matrix (5 variables, 10 samples)
  expr_mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:5)
  colnames(expr_mat) <- paste0("sample", 1:10)

  # Create sample info with custom survival column names
  sample_info <- tibble::tibble(
    sample = paste0("sample", 1:10),
    survival_time = rexp(10, rate = 0.2),  # Custom time column
    death_event = rbinom(10, 1, 0.6)       # Custom event column
  )

  # Create variable info
  var_info <- tibble::tibble(
    variable = paste0("var", 1:5),
    type = rep("gene", 5)
  )

  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )

  # Test with custom column names
  result <- suppressMessages(gly_cox(exp, time_col = "survival_time", event_col = "death_event"))

  # Check result structure
  expect_s3_class(result, c("glystats_cox_res", "glystats_res"))
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))
  expect_equal(nrow(result$tidy_result), 5)
  expect_true("coefficient" %in% colnames(result$tidy_result))
  expect_true("hr" %in% colnames(result$tidy_result))
  expect_true("p" %in% colnames(result$tidy_result))
})

test_that("gly_cox p-value adjustment methods work", {
  # Create test survival data
  set.seed(987)

  # Create expression matrix (8 variables, 15 samples)
  expr_mat <- matrix(rnorm(120), nrow = 8, ncol = 15)
  rownames(expr_mat) <- paste0("var", 1:8)
  colnames(expr_mat) <- paste0("sample", 1:15)

  # Create sample info with survival data
  sample_info <- tibble::tibble(
    sample = paste0("sample", 1:15),
    time = rexp(15, rate = 0.1),
    event = rbinom(15, 1, 0.7)
  )

  # Create variable info
  var_info <- tibble::tibble(
    variable = paste0("var", 1:8),
    type = rep("biomarker", 8)
  )

  # Create experiment object
  exp <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )

  # Test different p-value adjustment methods
  result_bh <- suppressMessages(gly_cox(exp, p_adj_method = "BH"))
  result_bonf <- suppressMessages(gly_cox(exp, p_adj_method = "bonferroni"))
  result_none <- suppressMessages(gly_cox(exp, p_adj_method = NULL))

  # Check that p_adj column exists for BH and bonferroni
  expect_true("p_adj" %in% colnames(result_bh$tidy_result))
  expect_true("p_adj" %in% colnames(result_bonf$tidy_result))

  # Check that p_adj column doesn't exist when p_adj_method is NULL
  expect_false("p_adj" %in% colnames(result_none$tidy_result))

  # Check that adjusted p-values are different between methods
  expect_false(identical(result_bh$tidy_result$p_adj, result_bonf$tidy_result$p_adj))
})

test_that("gly_cox error handling works", {
  # Create test data for error testing
  set.seed(111)

  # Create expression matrix
  expr_mat <- matrix(rnorm(30), nrow = 5, ncol = 6)
  rownames(expr_mat) <- paste0("var", 1:5)
  colnames(expr_mat) <- paste0("sample", 1:6)

  # Create sample info without survival columns
  sample_info_no_surv <- tibble::tibble(
    sample = paste0("sample", 1:6),
    group = rep(c("A", "B"), each = 3)
  )

  # Create sample info with survival columns
  sample_info_with_surv <- tibble::tibble(
    sample = paste0("sample", 1:6),
    time = rexp(6, rate = 0.2),
    event = rbinom(6, 1, 0.5)
  )

  # Create variable info
  var_info <- tibble::tibble(
    variable = paste0("var", 1:5),
    type = rep("gene", 5)
  )

  # Test error when time column doesn't exist
  exp_no_time <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info_no_surv,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )

  expect_error(suppressMessages(gly_cox(exp_no_time)), "Must be of type 'numeric'")

  # Test error when event column doesn't exist
  expect_error(suppressMessages(gly_cox(exp_no_time, event_col = "nonexistent")), "Must be of type 'numeric'")

  # Test error with invalid p_adj_method
  exp_valid <- glyexp::experiment(
    expr_mat = expr_mat,
    sample_info = sample_info_with_surv,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = "N"
  )

  expect_error(suppressMessages(gly_cox_(expr_mat, sample_info_with_surv$time, sample_info_with_surv$event, p_adj_method = "invalid")),
               "Must be element of set")
})

test_that("gly_cox works with test_gp_exp data", {
  # Create survival data for test_gp_exp
  set.seed(42)

  # Use a small subset for faster testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 10)

  # Add survival data to sample_info
  sample_info <- glyexp::get_sample_info(exp_small)
  n_samples <- nrow(sample_info)

  # Create realistic survival data based on groups
  survival_data <- sample_info |>
    dplyr::mutate(
      # Different baseline hazards for different groups
      base_hazard = dplyr::case_when(
        group == "H" ~ 0.05,
        group == "M" ~ 0.08,
        group == "Y" ~ 0.12,
        group == "C" ~ 0.15,
        TRUE ~ 0.1
      ),
      # Generate survival times
      time = rexp(n_samples, rate = base_hazard),
      # Generate events (80% event rate)
      event = rbinom(n_samples, 1, 0.8)
    ) |>
    dplyr::select(-base_hazard)

  # Update experiment with survival data
  exp_with_survival <- glyexp::experiment(
    expr_mat = glyexp::get_expr_mat(exp_small),
    sample_info = survival_data,
    var_info = glyexp::get_var_info(exp_small),
    exp_type = glyexp::get_exp_type(exp_small),
    glycan_type = glyexp::get_glycan_type(exp_small)
  )

  # Test gly_cox function
  result <- suppressMessages(gly_cox(exp_with_survival))

  # Check result structure
  expect_s3_class(result, c("glystats_cox_res", "glystats_res"))
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))
  expect_equal(nrow(result$tidy_result), 10)

  # Check that variable info was joined (add_info = TRUE by default)
  expect_true("protein" %in% colnames(result$tidy_result))
  expect_true("gene" %in% colnames(result$tidy_result))

  # Check that results are reasonable
  expect_true(all(is.finite(result$tidy_result$coefficient)))
  expect_true(all(result$tidy_result$p >= 0 & result$tidy_result$p <= 1))
  expect_true(all(result$tidy_result$p_adj >= 0 & result$tidy_result$p_adj <= 1))
})

test_that("gly_cox input validation works", {
  # Test with non-experiment object
  expect_error(gly_cox("not_an_experiment"), "Assertion on 'exp' failed")

  # Test gly_cox_ with invalid inputs
  expr_mat <- matrix(rnorm(20), nrow = 4, ncol = 5)
  time <- rexp(5, rate = 0.1)
  event <- rbinom(5, 1, 0.7)

  # Test with mismatched dimensions
  expect_error(gly_cox_(expr_mat, time[1:3], event), "must match length of time vector")
  expect_error(gly_cox_(expr_mat, time, event[1:3]), "must match length of event vector")

  # Test with empty matrix
  expect_error(gly_cox_(matrix(nrow = 0, ncol = 0), numeric(0), numeric(0)),
               "Assertion on 'expr_mat' failed")
})

test_that("gly_cox handles edge cases", {
  # Test with normal case (some events, some censoring) using larger sample size for stability
  set.seed(999)
  expr_mat <- matrix(rnorm(80), nrow = 4, ncol = 20)  # Increased sample size
  rownames(expr_mat) <- paste0("var", 1:4)
  colnames(expr_mat) <- paste0("sample", 1:20)

  # Create more realistic survival times and events
  time <- rexp(20, rate = 0.1)
  event <- rbinom(20, 1, 0.7)  # 70% event rate

  # This should work normally
  result <- suppressMessages(gly_cox_(expr_mat, time, event))
  expect_s3_class(result, c("glystats_cox_res", "glystats_res"))
  expect_type(result, "list")
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))
  expect_equal(nrow(result$tidy_result), 4)
  expect_true("variable" %in% colnames(result$tidy_result))

  # Test with all events = 1 (no censoring) - also use larger sample
  event_all <- rep(1, 20)
  result_all_events <- suppressMessages(gly_cox_(expr_mat, time, event_all))
  expect_s3_class(result_all_events, c("glystats_cox_res", "glystats_res"))
  expect_equal(nrow(result_all_events$tidy_result), 4)
})
