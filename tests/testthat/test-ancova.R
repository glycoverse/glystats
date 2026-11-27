skip_if_not_installed("emmeans")

test_that("gly_ancova works with ancova method", {
  # Use test_gp_exp with 3 groups for ANCOVA
  exp_3group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H", "M")) |>
    glyexp::slice_sample_var(n = 10) |>
    glyexp::mutate_obs(covar = seq_along(.data$group)) |>
    as_test_se()

  # Run DEA with ANCOVA
  result <- suppressMessages(gly_ancova(exp_3group, covariate_cols = "covar"))

  # Test core functionality - should return a list with tidy_result and raw_result
  expect_type(result, "list")
  expect_s3_class(result, c("glystats_ancova_res", "glystats_res"))
  expect_named(
    result,
    c("tidy_result", "raw_result", "meta_data", "preprocessing", "identifiers")
  )

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
  expect_true(all(
    c(".variable", "ref_group", "test_group", "p_val", "p_adj", "log2fc") %in%
      colnames(post_hoc_test)
  ))

  # Test raw_result structure
  expect_type(result$raw_result, "list")
  expect_named(result$raw_result, c("main_test", "post_hoc_test"))
})

test_that("gly_ancova accepts a covariate named subject", {
  set.seed(456)
  subjects <- seq_len(24)
  groups <- factor(rep(c("A", "B"), times = 12), levels = c("A", "B"))
  log_expression <-
    3 +
    1.5 * (groups == "B") +
    0.02 * subjects +
    stats::rnorm(24, sd = 0.05)
  expression <- rbind(V1 = 2^log_expression - 1e-6)
  colnames(expression) <- paste0("sample", subjects)
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expression),
    colData = S4Vectors::DataFrame(group = groups, subject = subjects)
  )

  result <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = "subject",
    add_info = FALSE,
    p_adj_method = NULL
  ))

  expect_s3_class(result$raw_result$main_test$V1, "aov")
  expect_true(any(is.finite(result$tidy_result$main_test$p_val)))
})

test_that("gly_ancova comparison direction is correct for 2 groups", {
  # Create a test experiment with 2 groups
  # group B has higher mean than group A
  set.seed(123)
  var_info <- tibble::tibble(variable = "V1")
  sample_info <- tibble::tibble(
    sample = paste0("S", 1:20),
    group = factor(rep(c("A", "B"), each = 10), levels = c("A", "B")),
    covar = rnorm(20)
  )
  values <- c(
    rnorm(10, mean = 5, sd = 0.2) + 0.1 * sample_info$covar[1:10],
    rnorm(10, mean = 8, sd = 0.2) + 0.1 * sample_info$covar[11:20]
  )
  expr_mat <- matrix(values, nrow = 1, byrow = TRUE)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, "others")

  # Call gly_ancova
  result <- suppressMessages(gly_ancova(exp, covariate_cols = "covar"))

  # Test post_hoc
  expect_identical(result$tidy_result$main_test$post_hoc, "A_vs_B")
  expect_true(all(result$tidy_result$post_hoc_test$log2fc > 0))
  expect_identical(result$tidy_result$post_hoc_test$ref_group, "A")
  expect_identical(result$tidy_result$post_hoc_test$test_group, "B")
})

test_that(".analyze_ancova works correctly", {
  # Create test data
  set.seed(123)
  expr_mat <- matrix(abs(rnorm(100)) + 1, nrow = 10, ncol = 10)
  rownames(expr_mat) <- paste0("var", 1:10)
  colnames(expr_mat) <- paste0("sample", 1:10)
  groups <- factor(rep(c("A", "B", "C"), c(3, 3, 4)))
  covariates <- data.frame(covar = rnorm(10), stringsAsFactors = FALSE)

  # Test function execution
  suppressMessages({
    result <- .analyze_ancova(expr_mat, groups, covariates)
  })

  # Verify results
  expect_s3_class(result, "glystats_ancova_res")
  expect_type(result, "list")
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result$main_test))
  expect_equal(nrow(result$tidy_result$main_test), 10)
})

test_that("gly_ancova error handling", {
  exp_small <- test_gp_exp |> glyexp::slice_sample_var(n = 5)

  # Missing covariate columns
  expect_error(
    suppressMessages(gly_ancova(exp_small, covariate_cols = "nonexistent")),
    "covariate_cols not found"
  )

  # Missing covariates for .analyze_ancova
  expr_mat <- glyexp::get_expr_mat(exp_small)
  groups <- factor(glyexp::get_sample_info(exp_small)$group)
  expect_error(
    suppressMessages(.analyze_ancova(expr_mat, groups, covariates = NULL)),
    "covariates must be provided"
  )
})

test_that("ANCOVA screens features on the covariate-adjusted group effect", {
  group <- factor(rep(c("A", "B"), each = 12))
  age <- rep(1:12, 2) + 10 * (group == "B")
  confounded <- age + 0.3 * rep(c(1, -1, -1, 1), 6)
  values <- rbind(
    confounded = confounded,
    signal = confounded + 4 * (group == "B")
  )
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(group = group, age = age)
  )

  res <- gly_ancova(exp, covariate_cols = "age", log_transform = FALSE)
  expected <- apply(values, 1, function(y) {
    stats::coef(summary(stats::lm(y ~ group + age)))["groupB", "Pr(>|t|)"]
  })
  main <- res$tidy_result$main_test
  expect_equal(main$p_val, unname(expected[main$.variable]))
  expect_equal(
    main$p_adj,
    unname(stats::p.adjust(expected, "BH")[main$.variable])
  )
  expect_identical(unique(res$tidy_result$post_hoc_test$.variable), "signal")

  exp$batch <- group
  aliased <- gly_ancova(exp, covariate_cols = "batch", log_transform = FALSE)
  expect_identical(aliased$tidy_result$main_test$.variable, rownames(exp))
  expect_true(all(is.na(aliased$tidy_result$main_test$p_val)))
  expect_equal(nrow(aliased$tidy_result$post_hoc_test), 0L)
})

test_that("ANCOVA models covariate names literally", {
  exp <- model_argument_experiment()
  exp$age <- factor(exp$age > 0)
  expected <- gly_ancova(exp, covariate_cols = "age")
  for (name in c(
    "variable",
    "sample",
    "log_value",
    "batch score",
    "age+group"
  )) {
    test_exp <- exp
    names(SummarizedExperiment::colData(test_exp))[3] <- name
    result <- gly_ancova(
      test_exp,
      covariate_cols = name,
      contrasts = stats::setNames(list("contr.sum"), name)
    )
    expect_equal(
      result$tidy_result$main_test$p_val,
      expected$tidy_result$main_test$p_val
    )
    expect_equal(
      result$tidy_result$post_hoc_test$p_adj,
      expected$tidy_result$post_hoc_test$p_adj
    )
    expect_true("contr.sum" %in% result$raw_result$main_test[[1]]$contrasts)
  }
})
