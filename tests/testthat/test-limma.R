skip_if_not_installed("limma")

test_that("gly_limma works with 2 groups", {
  # Use test_gp_exp and filter to 2 groups for limma
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10) |>
    as_test_se() # Use smaller subset for faster testing

  # Run DEA with limma
  result <- suppressMessages(gly_limma(exp_2group))

  # Test core functionality
  expect_s3_class(result, c("glystats_limma_res", "glystats_res"))
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))
  expect_equal(nrow(result$tidy_result), 10)
  expect_true("p_adj" %in% colnames(result$tidy_result)) # p_adj should exist
  expect_true("log2fc" %in% colnames(result$tidy_result)) # log2fc should exist
  expect_type(result$tidy_result$log2fc, "double") # log2fc should be numeric
  expect_true("t" %in% colnames(result$tidy_result)) # t-statistic
  expect_true("b" %in% colnames(result$tidy_result)) # log-odds (B-statistic)
})

test_that("gly_limma direction is correct for 2 groups", {
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

  # Call gly_limma
  result <- suppressMessages(gly_limma(exp))

  # Test post_hoc
  expect_true(result$tidy_result$log2fc > 0)
  expect_true(result$tidy_result$t > 0)
})

test_that("gly_limma respects custom contrasts for 2 groups", {
  # Create a test experiment with 2 groups
  # group B has higher mean than group A
  set.seed(123)
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

  # Call gly_limma with reversed contrast
  result <- suppressMessages(gly_limma(exp, contrasts = "B_vs_A"))

  # Test post_hoc (A vs B should be negative)
  expect_true(result$tidy_result$log2fc < 0)
  expect_true(result$tidy_result$t < 0)
})

test_that("gly_limma error handling", {
  # Use test_gp_exp for error testing
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)

  # Test various error conditions - group column not found
  expect_error(
    suppressMessages(gly_limma(exp_small, group_col = "nonexistent")),
    "not found in sample information"
  )
})

test_that("gly_limma group validation", {
  # Test with 1 group
  exp_1group <- test_gp_exp |>
    glyexp::filter_obs(group == "C") |>
    glyexp::slice_sample_var(n = 5)

  # Test 1 group should still fail (need at least 2 groups)
  expect_error(suppressMessages(gly_limma(exp_1group)), "at least 2 levels")
})

test_that("gly_limma ref_group parameter works", {
  # Use test_gp_exp and filter to 2 groups for limma
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)

  # Test with default (no ref_group)
  result_default <- suppressMessages(gly_limma(exp_2group))

  # Test with ref_group = "H" (should reverse the comparison)
  result_ref_h <- suppressMessages(gly_limma(exp_2group, ref_group = "H"))

  # Test with ref_group = "C" (should be same as default since C is first alphabetically)
  result_ref_c <- suppressMessages(gly_limma(exp_2group, ref_group = "C"))

  # Check that log2fc values are negated when reference group changes
  expect_equal(
    result_default$tidy_result$log2fc,
    -result_ref_h$tidy_result$log2fc,
    tolerance = 1e-10
  )
  expect_equal(
    result_default$tidy_result$log2fc,
    result_ref_c$tidy_result$log2fc,
    tolerance = 1e-10
  )

  # Test invalid ref_group
  expect_error(
    suppressMessages(gly_limma(exp_2group, ref_group = "invalid")),
    "Must be element of set"
  )
})

test_that("gly_limma works with multi-group data", {
  # Test with 3 groups (using C, H, M from test_gp_exp)
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::mutate_obs(group = factor(group, levels = c("H", "M", "C"))) |>
    glyexp::slice_sample_var(n = 10) # Use smaller subset for faster testing

  result <- suppressMessages(gly_limma(exp_3group))

  # Test core functionality
  expect_s3_class(result, c("glystats_limma_res", "glystats_res"))
  expect_true(tibble::is_tibble(result$tidy_result))
  expect_true("ref_group" %in% colnames(result$tidy_result)) # Should have contrast column
  expect_true("test_group" %in% colnames(result$tidy_result)) # Should have contrast column
  expect_true("log2fc" %in% colnames(result$tidy_result))
  expect_true("p_adj" %in% colnames(result$tidy_result))

  # Should have 3 pairwise comparisons: H_vs_C, H_vs_M, M_vs_C
  contrasts <- stringr::str_c(
    result$tidy_result$ref_group,
    "_vs_",
    result$tidy_result$test_group
  )
  expect_setequal(unique(contrasts), c("H_vs_C", "H_vs_M", "M_vs_C"))

  # Each contrast should have the same number of variables
  expect_equal(nrow(result$tidy_result), 10 * 3) # 10 variables * 3 contrasts
})

test_that("gly_limma direction is correct for 3 groups", {
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

  # Call gly_limma
  result <- suppressMessages(gly_limma(exp))

  # Test post_hoc
  expect_true(all(result$tidy_result$log2fc > 0))
  expect_true(all(result$tidy_result$t > 0))
})

test_that("gly_limma direction is correct for 3 groups with custom contrasts", {
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

  # Call gly_limma
  custom_contrasts <- c("A_vs_B", "A_vs_C")
  result <- suppressMessages(gly_limma(exp, contrasts = custom_contrasts))

  # Test post_hoc
  expect_true(all(result$tidy_result$log2fc > 0))
  expect_true(all(result$tidy_result$t > 0))
})

test_that("gly_limma custom contrasts work with hyphen format", {
  # Test with 4 groups using custom contrasts (hyphen format)
  exp_4group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M", "Y")) |>
    glyexp::slice_sample_var(n = 5)

  # Test custom contrasts with hyphen format
  custom_contrasts <- c("H-C", "H-M", "H-Y")
  result <- suppressMessages(gly_limma(
    exp_4group,
    contrasts = custom_contrasts
  ))

  # Should have exactly 3 contrasts as specified
  contrasts <- stringr::str_c(
    result$tidy_result$ref_group,
    "_vs_",
    result$tidy_result$test_group
  )
  expect_setequal(unique(contrasts), c("H_vs_C", "H_vs_M", "H_vs_Y"))
  expect_equal(nrow(result$tidy_result), 5 * 3) # 5 variables * 3 contrasts
})

test_that("gly_limma custom contrasts work with _vs_ format", {
  # Test with 4 groups using custom contrasts (_vs_ format)
  exp_4group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M", "Y")) |>
    glyexp::slice_sample_var(n = 5)

  # Test custom contrasts with _vs_ format
  custom_contrasts <- c("H_vs_C", "M_vs_C", "Y_vs_C")
  result <- suppressMessages(gly_limma(
    exp_4group,
    contrasts = custom_contrasts
  ))

  # Should have exactly 3 contrasts as specified
  contrasts <- stringr::str_c(
    result$tidy_result$ref_group,
    "_vs_",
    result$tidy_result$test_group
  )
  expect_setequal(contrasts, c("H_vs_C", "M_vs_C", "Y_vs_C"))
  expect_equal(nrow(result$tidy_result), 5 * 3) # 5 variables * 3 contrasts
})

test_that("gly_limma works with covariate_cols", {
  set.seed(101)
  var_info <- tibble::tibble(variable = paste0("V", 1:5))
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:20),
    group = factor(rep(c("A", "B"), each = 10), levels = c("A", "B")),
    batch = rep(c("X", "Y"), times = 10),
    age = seq(30, 49)
  )
  expr_mat <- matrix(rnorm(5 * 20, mean = 1, sd = 0.2), nrow = 5)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

  result <- suppressMessages(gly_limma(exp, covariate_cols = c("batch", "age")))

  expect_s3_class(result, c("glystats_limma_res", "glystats_res"))
  expect_equal(nrow(result$tidy_result), 5)
  expect_true("log2fc" %in% colnames(result$tidy_result))
})

test_that("gly_limma handles all-NA variables without covariates", {
  set.seed(303)
  var_info <- tibble::tibble(variable = c("V1", "V2"))
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:6),
    group = factor(rep(c("A", "B"), each = 3), levels = c("A", "B"))
  )
  expr_mat <- matrix(rnorm(2 * 6, mean = 1, sd = 0.2), nrow = 2)
  expr_mat[1, ] <- NA_real_
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

  result <- suppressMessages(gly_limma(exp))

  expect_s3_class(result, c("glystats_limma_res", "glystats_res"))
  expect_setequal(result$tidy_result$variable, var_info$variable)
})

test_that("gly_limma falls back when eBayes trend fails", {
  set.seed(404)
  expr_mat <- matrix(abs(rnorm(20)) + 1, nrow = 4)
  rownames(expr_mat) <- paste0("V", 1:4)
  colnames(expr_mat) <- paste0("S", 1:5)
  groups <- factor(rep(c("A", "B"), length.out = 5), levels = c("A", "B"))

  result <- NULL
  trend_calls <- logical(0)
  real_ebayes <- limma::eBayes

  with_mocked_bindings(
    {
      expect_warning(
        {
          result <- suppressMessages(.analyze_limma(expr_mat, groups))
        },
        "trend = FALSE"
      )
    },
    eBayes = function(fit, trend = FALSE, ...) {
      trend_calls <<- c(trend_calls, trend)
      if (isTRUE(trend)) {
        stop("Problem with covariate")
      }
      real_ebayes(fit, trend = trend, ...)
    },
    .package = "limma"
  )

  expect_s3_class(result, c("glystats_limma_res", "glystats_res"))
  expect_true(any(trend_calls))
  expect_true(any(!trend_calls))
})

test_that(".analyze_limma works with covariates", {
  set.seed(202)
  var_info <- tibble::tibble(variable = paste0("V", 1:4))
  sample_names <- paste0("S", 1:16)
  groups <- factor(rep(c("A", "B"), each = 8), levels = c("A", "B"))
  expr_mat <- matrix(rnorm(4 * 16, mean = 1.5, sd = 0.3), nrow = 4)
  colnames(expr_mat) <- sample_names
  rownames(expr_mat) <- var_info$variable
  covariates <- data.frame(
    batch = rep(c("B1", "B2"), times = 8),
    age = seq(40, 55),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(.analyze_limma(
    expr_mat,
    groups,
    covariates = covariates
  ))

  expect_s3_class(result, c("glystats_limma_res", "glystats_res"))
  expect_equal(nrow(result$tidy_result), 4)
})

test_that(".analyze_limma validates covariate rows", {
  expr_mat <- matrix(rnorm(3 * 10), nrow = 3)
  groups <- factor(rep(c("A", "B"), each = 5), levels = c("A", "B"))
  covariates <- data.frame(batch = rep(c("B1", "B2"), each = 4))

  expect_error(
    suppressMessages(.analyze_limma(expr_mat, groups, covariates = covariates)),
    "covariates must have"
  )
})

test_that("gly_limma contrasts error handling works", {
  # Test with 3 groups
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)

  # Test invalid group names
  expect_error(
    suppressMessages(gly_limma(exp_3group, contrasts = c("H-C", "X-Y"))),
    "Group.*not found in data"
  )

  # Test invalid contrast format
  expect_error(
    suppressMessages(gly_limma(exp_3group, contrasts = c("H_C"))),
    "Invalid contrast format"
  )
})

test_that("gly_limma handles group names with hyphens correctly", {
  # Create test data with group names containing hyphens
  exp_hyphen <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 5)

  # Modify group names to include hyphens
  sample_info_modified <- glyexp::get_sample_info(exp_hyphen)
  sample_info_modified$group <- factor(
    ifelse(
      sample_info_modified$group == "C",
      "Control-1",
      ifelse(sample_info_modified$group == "H", "High-dose", "Medium-dose")
    ),
    levels = c("Control-1", "High-dose", "Medium-dose")
  )

  # Create new experiment with modified sample info
  exp_hyphen_modified <- glyexp::experiment(
    expr_mat = glyexp::get_expr_mat(exp_hyphen),
    sample_info = sample_info_modified,
    var_info = glyexp::get_var_info(exp_hyphen),
    exp_type = "glycoproteomics",
    glycan_type = "N"
  )

  # Test that hyphen format fails with helpful error message
  expect_error(
    suppressMessages(gly_limma(
      exp_hyphen_modified,
      contrasts = c("High-dose-Control-1")
    )),
    "Group names contain hyphens.*Use the format.*_vs_"
  )

  # Test that _vs_ format works
  result <- suppressMessages(gly_limma(
    exp_hyphen_modified,
    contrasts = c("High-dose_vs_Control-1")
  ))
  contrasts <- stringr::str_c(
    result$tidy_result$ref_group,
    "_vs_",
    result$tidy_result$test_group
  )
  expect_setequal(contrasts, c("High-dose_vs_Control-1"))
})

test_that(".analyze_limma works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B"), each = 5))

  # Test function execution
  suppressMessages({
    result <- .analyze_limma(expr_mat, groups)
  })

  # Verify results
  expect_s3_class(result, c("glystats_limma_res", "glystats_res"))
  expect_true("tidy_result" %in% names(result))
  expect_true("raw_result" %in% names(result))
  expect_true(tibble::is_tibble(result$tidy_result))
  expect_true("log2fc" %in% colnames(result$tidy_result))
  expect_true("p_val" %in% colnames(result$tidy_result)) # Now uses p_val
  expect_equal(nrow(result$tidy_result), 10)
})

test_that("gly_limma supports subject_cols for paired design", {
  set.seed(42)
  n_subjects <- 6
  subjects <- factor(rep(paste0("S", seq_len(n_subjects)), each = 2))
  groups <- factor(rep(c("A", "B"), times = n_subjects), levels = c("A", "B"))
  subject_effect <- rep(seq_len(n_subjects) * 2, each = 2)
  group_effect <- ifelse(groups == "B", 1.5, 0)
  n_genes <- 4
  n_samples <- length(groups)
  noise <- matrix(rnorm(n_genes * n_samples, sd = 0.05), nrow = n_genes)
  expr_mat <- matrix(10, nrow = n_genes, ncol = n_samples) +
    matrix(rep(subject_effect + group_effect, each = n_genes), nrow = n_genes) +
    noise
  rownames(expr_mat) <- paste0("V", seq_len(n_genes))
  sample_info <- tibble::tibble(
    sample = paste0("sample", seq_len(n_samples)),
    group = groups,
    subject = subjects
  )
  colnames(expr_mat) <- sample_info$sample
  exp <- glyexp::experiment(
    expr_mat,
    sample_info,
    tibble::tibble(variable = rownames(expr_mat)),
    "others"
  )

  result_paired <- suppressMessages(gly_limma(
    exp,
    subject_col = "subject",
    add_info = FALSE
  ))
  result_unpaired <- suppressMessages(gly_limma(exp, add_info = FALSE))

  expect_s3_class(result_paired, c("glystats_limma_res", "glystats_res"))
  expect_true(
    mean(abs(result_paired$tidy_result$t)) >
      mean(abs(result_unpaired$tidy_result$t))
  )
})

test_that(".analyze_limma supports subjects for paired design", {
  set.seed(101)
  n_subjects <- 5
  subjects <- factor(rep(paste0("S", seq_len(n_subjects)), each = 2))
  groups <- factor(rep(c("A", "B"), times = n_subjects), levels = c("A", "B"))
  subject_effect <- rep(seq_len(n_subjects) * 1.5, each = 2)
  group_effect <- ifelse(groups == "B", 1.2, 0)
  n_genes <- 3
  n_samples <- length(groups)
  noise <- matrix(rnorm(n_genes * n_samples, sd = 0.05), nrow = n_genes)
  expr_mat <- matrix(8, nrow = n_genes, ncol = n_samples) +
    matrix(rep(subject_effect + group_effect, each = n_genes), nrow = n_genes) +
    noise
  rownames(expr_mat) <- paste0("V", seq_len(n_genes))
  colnames(expr_mat) <- paste0("sample", seq_len(n_samples))

  result_paired <- suppressMessages(.analyze_limma(
    expr_mat,
    groups,
    subjects = subjects
  ))
  result_unpaired <- suppressMessages(.analyze_limma(expr_mat, groups))

  expect_s3_class(result_paired, c("glystats_limma_res", "glystats_res"))
  expect_true(
    mean(abs(result_paired$tidy_result$t)) >
      mean(abs(result_unpaired$tidy_result$t))
  )
})

test_that(".analyze_limma validates subjects length", {
  expr_mat <- matrix(rnorm(3 * 6), nrow = 3)
  groups <- factor(rep(c("A", "B"), each = 3), levels = c("A", "B"))
  subjects <- factor(rep(c("S1", "S2"), each = 2))

  expect_error(
    suppressMessages(.analyze_limma(expr_mat, groups, subjects = subjects)),
    "subjects must have"
  )
})

test_that("gly_limma returns meta_data from experiment", {
  # Use test_gp_exp and filter to 2 groups for limma
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 10)

  result <- suppressMessages(gly_limma(exp_2group))

  # Check that meta_data exists and contains expected values
  expect_true("meta_data" %in% names(result))
  expect_equal(result$meta_data$exp_type, "glycoproteomics")
  expect_equal(result$meta_data$glycan_type, "N")
})

test_that(".analyze_limma does not return meta_data", {
  # Use test_gp_exp and extract expression matrix
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H"))
  expr_mat <- glyexp::get_expr_mat(exp_2group)
  groups <- factor(rep(c("C", "H"), each = 3))

  result <- suppressMessages(.analyze_limma(expr_mat, groups))

  # Check that meta_data does NOT exist
  expect_false("meta_data" %in% names(result))
})
