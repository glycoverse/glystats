make_correlated_set_exp <- function() {
  set.seed(20260817)
  n <- 60
  target <- matrix(
    c(
      1,
      0.92,
      0.70,
      0.92,
      1,
      0.90,
      0.70,
      0.90,
      1
    ),
    nrow = 3,
    byrow = TRUE
  )
  basis <- scale(matrix(stats::rnorm(n * 3), nrow = n), scale = FALSE)
  basis <- qr.Q(qr(basis))
  profiles <- basis %*% chol(target) * sqrt(n - 1)
  colnames(profiles) <- c("A", "B", "C")

  log_expression <- rbind(
    A = profiles[, "A"],
    B = profiles[, "B"],
    C = profiles[, "C"],
    A_alias = profiles[, "A"],
    constant = rep(1, n),
    missing = rep(NA_real_, n)
  )
  expression <- 2^log_expression - 1e-6
  colnames(expression) <- paste0("S", seq_len(n))

  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expression),
    rowData = S4Vectors::DataFrame(
      site = c("one", "one", "one", "one", "one", "one"),
      row.names = rownames(expression)
    ),
    colData = S4Vectors::DataFrame(
      group = factor(
        rep(c("control", "case"), each = n / 2),
        levels = c("control", "case")
      ),
      row.names = colnames(expression)
    ),
    metadata = list(study = "correlated-set-test")
  )
}

make_hotelling_exp <- function(paired = FALSE) {
  set.seed(42)
  if (paired) {
    n <- 14
    control <- cbind(
      A = stats::rnorm(n),
      B = stats::rnorm(n),
      C = stats::rnorm(n)
    )
    differences <- cbind(
      A = stats::rnorm(n, 0.5, 0.3),
      B = stats::rnorm(n, -0.25, 0.4),
      C = stats::rnorm(n, 0.1, 0.5)
    )
    case <- control + differences
    log_expression <- t(rbind(control, case))
    groups <- factor(
      rep(c("control", "case"), each = n),
      levels = c("control", "case")
    )
    subjects <- rep(paste0("P", seq_len(n)), 2)
  } else {
    n <- 24
    control <- cbind(
      A = stats::rnorm(n),
      B = 0.4 * stats::rnorm(n) + stats::rnorm(n),
      C = stats::rnorm(n)
    )
    case <- cbind(
      A = stats::rnorm(n, 0.6),
      B = 0.4 * stats::rnorm(n) + stats::rnorm(n, -0.3),
      C = stats::rnorm(n, 0.1)
    )
    log_expression <- t(rbind(control, case))
    groups <- factor(
      rep(c("control", "case"), each = n),
      levels = c("control", "case")
    )
    subjects <- paste0("S", seq_len(n * 2))
  }

  log_expression <- rbind(
    log_expression,
    A_alias = log_expression["A", ],
    sum_ab = log_expression["A", ] + log_expression["B", ],
    missing1 = NA_real_,
    missing2 = NA_real_
  )
  expression <- 2^log_expression - 1e-6
  colnames(expression) <- paste0("S", seq_len(ncol(expression)))

  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expression),
    colData = S4Vectors::DataFrame(
      group = groups,
      subject = subjects,
      row.names = colnames(expression)
    ),
    metadata = list(study = "hotelling-test")
  )
}

make_matrix_set_exp <- function(ref, test, paired = FALSE) {
  stopifnot(identical(colnames(ref), colnames(test)))
  n_ref <- nrow(ref)
  n_test <- nrow(test)
  log_expression <- t(rbind(ref, test))
  expression <- 2^log_expression - 1e-6
  colnames(expression) <- paste0("S", seq_len(ncol(expression)))
  sample_info <- S4Vectors::DataFrame(
    group = factor(
      rep(c("control", "case"), c(n_ref, n_test)),
      levels = c("control", "case")
    ),
    row.names = colnames(expression)
  )
  if (paired) {
    stopifnot(n_ref == n_test)
    sample_info$subject <- rep(paste0("P", seq_len(n_ref)), 2)
  }

  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expression),
    colData = sample_info
  )
}

test_that("gly_set_test constructs complete and connected sets", {
  exp <- make_correlated_set_exp()

  complete <- gly_set_test(
    exp,
    threshold = 0.85,
    clustering = "complete",
    p_adj_method = NULL
  )
  connected <- gly_set_test(
    exp,
    threshold = 0.85,
    p_adj_method = NULL
  )

  complete_sets <- complete$raw_result$set_construction
  connected_sets <- connected$raw_result$set_construction
  expect_s3_class(complete, c("glystats_set_test_res", "glystats_res"))
  expect_identical(complete_sets$sets$set_1, c("A", "B", "A_alias"))
  expect_identical(complete_sets$sets$set_2, "C")
  expect_identical(connected_sets$sets$set_1, c("A", "B", "C", "A_alias"))
  expect_named(
    complete_sets,
    c(
      "sets",
      "membership",
      "correlation_matrix",
      "excluded_variables",
      "threshold",
      "correlation",
      "clustering",
      "within"
    )
  )
  expect_equal(complete$meta_data, S4Vectors::metadata(exp))
})

test_that("gly_set_test records exclusions and honors within strata", {
  exp <- make_correlated_set_exp()
  SummarizedExperiment::rowData(exp)$site <- c(
    "one",
    "one",
    "two",
    "one",
    "one",
    "one"
  )

  result <- gly_set_test(
    exp,
    threshold = 0.85,
    clustering = "connected",
    within = "site",
    p_adj_method = NULL
  )
  construction <- result$raw_result$set_construction

  expect_identical(construction$sets$set_1, c("A", "B", "A_alias"))
  expect_identical(construction$sets$set_2, "C")
  expect_identical(
    is.na(construction$correlation_matrix["A", "C"]),
    TRUE
  )
  expect_equal(
    construction$excluded_variables,
    tibble::tibble(
      variable = c("constant", "missing"),
      reason = c("zero_variance", "all_missing")
    )
  )
  expect_identical(
    construction$membership,
    tibble::tibble(
      set_id = c("set_1", "set_1", "set_1", "set_2"),
      variable = c("A", "B", "A_alias", "C")
    )
  )
})

test_that("gly_set_test matches an independent Hotelling calculation", {
  exp <- make_hotelling_exp()
  result <- gly_set_test(exp, list(signal = c("A", "B")))

  log_expression <- log2(SummarizedExperiment::assay(exp) + 1e-6)
  groups <- SummarizedExperiment::colData(exp)$group
  ref <- t(log_expression[c("A", "B"), groups == "control"])
  test <- t(log_expression[c("A", "B"), groups == "case"])
  n_ref <- nrow(ref)
  n_test <- nrow(test)
  p <- ncol(ref)
  covariance <- ((n_ref - 1) *
    stats::cov(ref) +
    (n_test - 1) * stats::cov(test)) /
    (n_ref + n_test - 2)
  estimate <- colMeans(test) - colMeans(ref)
  distance_squared <- as.numeric(
    t(estimate) %*% solve(covariance) %*% estimate
  )
  expected_t2 <- n_ref * n_test / (n_ref + n_test) * distance_squared
  expected_df2 <- n_ref + n_test - p - 1
  expected_f <- expected_df2 * expected_t2 / ((n_ref + n_test - 2) * p)

  expect_s3_class(result, c("glystats_set_test_res", "glystats_res"))
  expect_named(result, c("tidy_result", "raw_result", "meta_data"))
  expect_named(result$tidy_result, c("sets", "members"))
  expect_equal(result$tidy_result$sets$statistic, expected_t2)
  expect_equal(result$tidy_result$sets$df1, p)
  expect_equal(result$tidy_result$sets$df2, expected_df2)
  expect_equal(result$tidy_result$sets$effect_size, sqrt(distance_squared))
  expect_equal(
    result$tidy_result$sets$p_val,
    stats::pf(expected_f, p, expected_df2, lower.tail = FALSE)
  )
  expect_identical(result$tidy_result$sets$status, "ok")
  expect_identical(
    get_tidy_result(result, "members"),
    result$tidy_result$members
  )
  expect_equal(result$meta_data, S4Vectors::metadata(exp))
})

test_that("gly_set_test matches a paired Hotelling calculation", {
  exp <- make_hotelling_exp(paired = TRUE)
  result <- gly_set_test(
    exp,
    list(signal = c("A", "B")),
    subject_col = "subject"
  )

  log_expression <- log2(SummarizedExperiment::assay(exp) + 1e-6)
  groups <- SummarizedExperiment::colData(exp)$group
  differences <- t(
    log_expression[c("A", "B"), groups == "case"] -
      log_expression[c("A", "B"), groups == "control"]
  )
  n <- nrow(differences)
  p <- ncol(differences)
  estimate <- colMeans(differences)
  covariance <- stats::cov(differences)
  distance_squared <- as.numeric(
    t(estimate) %*% solve(covariance) %*% estimate
  )
  expected_t2 <- n * distance_squared
  expected_f <- (n - p) * expected_t2 / (p * (n - 1))

  expect_equal(result$tidy_result$sets$statistic, expected_t2)
  expect_equal(result$tidy_result$sets$df1, p)
  expect_equal(result$tidy_result$sets$df2, n - p)
  expect_equal(
    result$tidy_result$sets$p_val,
    stats::pf(expected_f, p, n - p, lower.tail = FALSE)
  )
  expect_identical(result$tidy_result$sets$paired, TRUE)
  expect_identical(
    result$raw_result$tests$signal$included_subjects,
    paste0("P", seq_len(n))
  )
})

test_that("gly_set_test tests identical profiles in their effective subspace", {
  exp <- make_hotelling_exp()
  redundant <- gly_set_test(
    exp,
    list(signal = c("A", "A_alias", "B")),
    p_adj_method = NULL
  )
  nonredundant <- gly_set_test(
    exp,
    list(signal = c("A", "B")),
    p_adj_method = NULL
  )

  expect_identical(redundant$tidy_result$sets$n_variables, 3L)
  expect_identical(redundant$tidy_result$sets$test_dimension, 2L)
  expect_equal(
    redundant$tidy_result$sets$statistic,
    nonredundant$tidy_result$sets$statistic
  )
  expect_identical(
    names(redundant$tidy_result$sets$estimate[[1]]),
    c("A", "A_alias", "B")
  )
  expect_equal(
    unname(redundant$tidy_result$sets$estimate[[1]][c("A", "A_alias")]),
    rep(redundant$tidy_result$sets$estimate[[1]][["A"]], 2)
  )
  expect_setequal(
    redundant$tidy_result$members$variable,
    c("A", "A_alias", "B")
  )
  expect_disjoint(
    names(redundant$raw_result$tests$signal),
    c("representatives", "aliases")
  )
  expect_disjoint(colnames(redundant$tidy_result$sets), "p_adj")
})

test_that("gly_set_test determines effective rank after complete-case selection", {
  for (paired in c(FALSE, TRUE)) {
    exp <- make_hotelling_exp(paired = paired)
    expression <- SummarizedExperiment::assay(exp)
    partial_alias <- expression["A", ]
    partial_alias[[1]] <- NA_real_
    expression <- rbind(expression, A_partial = partial_alias)
    exp <- SummarizedExperiment::SummarizedExperiment(
      assays = list(expression = expression),
      colData = SummarizedExperiment::colData(exp),
      metadata = S4Vectors::metadata(exp)
    )

    result <- gly_set_test(
      exp,
      list(signal = c("A", "A_partial", "B")),
      subject_col = if (paired) "subject" else NULL
    )

    expect_identical(result$tidy_result$sets$status, "ok")
    expect_identical(result$tidy_result$sets$test_dimension, 2L)
    expect_identical(
      names(result$tidy_result$sets$estimate[[1]]),
      c("A", "A_partial", "B")
    )
  }

  exp <- make_correlated_set_exp()
  expression <- SummarizedExperiment::assay(exp)
  partial_alias <- expression["A", ]
  partial_alias[[1]] <- NA_real_
  expression <- rbind(expression, A_partial = partial_alias)
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expression),
    colData = SummarizedExperiment::colData(exp),
    metadata = S4Vectors::metadata(exp)
  )
  result <- gly_set_test(
    exp,
    threshold = 0.85,
    clustering = "complete"
  )

  expect_contains(
    result$raw_result$set_construction$sets$set_1,
    "A_partial"
  )
  expect_identical(result$tidy_result$sets$status[[1]], "ok")
  expect_identical(result$tidy_result$sets$test_dimension[[1]], 2L)
})

test_that("gly_set_test tests general linear dependencies by effective rank", {
  exp <- make_hotelling_exp()
  dependent <- gly_set_test(
    exp,
    list(signal = c("A", "B", "sum_ab")),
    p_adj_method = NULL
  )
  nonredundant <- gly_set_test(
    exp,
    list(signal = c("A", "B")),
    p_adj_method = NULL
  )

  expect_identical(dependent$tidy_result$sets$status, "ok")
  expect_identical(dependent$tidy_result$sets$test_dimension, 2L)
  expect_equal(
    dependent$tidy_result$sets$statistic,
    nonredundant$tidy_result$sets$statistic
  )
})

test_that("gly_set_test rejects contrasts outside the covariance subspace", {
  basis <- matrix(c(1, 1) / sqrt(2), ncol = 1)
  expect_true(.estimate_in_covariance_subspace(c(1e-12, 1e-12), basis))
  expect_false(.estimate_in_covariance_subspace(c(0, 1e-12), basis))

  profile <- seq(-2.5, 2.5)

  for (paired in c(FALSE, TRUE)) {
    if (paired) {
      ref <- cbind(A = profile, B = -profile)
      test <- ref + cbind(A = profile, B = profile + 1)
    } else {
      ref <- cbind(A = profile, B = profile)
      test <- cbind(A = profile, B = profile + 1)
    }
    exp <- make_matrix_set_exp(ref, test, paired)
    result <- gly_set_test(
      exp,
      list(signal = c("A", "B")),
      subject_col = if (paired) "subject" else NULL
    )
    sets <- result$tidy_result$sets

    expect_identical(sets$status, "failed")
    expect_identical(sets$test_dimension, 1L)
    expect_identical(sets$statistic, NA_real_)
    expect_identical(sets$p_val, NA_real_)
    expect_equal(unname(sets$estimate[[1]]), c(0, 1), tolerance = 1e-12)
    expect_identical(
      sets$failure_reason,
      "The mean difference lies outside the estimable covariance subspace."
    )
  }
})

test_that("gly_set_test rejects sample-limited covariance rank", {
  set.seed(20260825)
  independent <- make_matrix_set_exp(
    matrix(stats::rnorm(15), nrow = 3, dimnames = list(NULL, LETTERS[1:5])),
    matrix(stats::rnorm(15), nrow = 3, dimnames = list(NULL, LETTERS[1:5]))
  )
  independent_result <- gly_set_test(
    independent,
    list(signal = LETTERS[1:5])
  )

  ref <- matrix(
    stats::rnorm(16),
    nrow = 4,
    dimnames = list(NULL, LETTERS[1:4])
  )
  differences <- matrix(
    stats::rnorm(16),
    nrow = 4,
    dimnames = list(NULL, LETTERS[1:4])
  )
  paired <- make_matrix_set_exp(ref, ref + differences, paired = TRUE)
  paired_result <- gly_set_test(
    paired,
    list(signal = LETTERS[1:4]),
    subject_col = "subject"
  )

  expect_identical(independent_result$tidy_result$sets$status, "failed")
  expect_identical(independent_result$tidy_result$sets$test_dimension, 4L)
  expect_identical(independent_result$tidy_result$sets$df2, 0)
  expect_identical(
    independent_result$tidy_result$sets$failure_reason,
    "The number of variables is too large for the complete group sample sizes."
  )
  expect_identical(paired_result$tidy_result$sets$status, "failed")
  expect_identical(paired_result$tidy_result$sets$test_dimension, 3L)
  expect_identical(paired_result$tidy_result$sets$df2, 0)
  expect_identical(
    paired_result$tidy_result$sets$failure_reason,
    "The number of variables must be smaller than the number of complete pairs."
  )
})

test_that("gly_set_test retains the classical sample-size boundary", {
  set.seed(20260826)
  independent <- make_matrix_set_exp(
    matrix(stats::rnorm(12), nrow = 3, dimnames = list(NULL, LETTERS[1:4])),
    matrix(stats::rnorm(12), nrow = 3, dimnames = list(NULL, LETTERS[1:4]))
  )
  independent_result <- gly_set_test(
    independent,
    list(signal = LETTERS[1:4])
  )

  ref <- matrix(
    stats::rnorm(20),
    nrow = 5,
    dimnames = list(NULL, LETTERS[1:4])
  )
  differences <- matrix(
    stats::rnorm(20),
    nrow = 5,
    dimnames = list(NULL, LETTERS[1:4])
  )
  paired <- make_matrix_set_exp(ref, ref + differences, paired = TRUE)
  paired_result <- gly_set_test(
    paired,
    list(signal = LETTERS[1:4]),
    subject_col = "subject"
  )

  expect_identical(independent_result$tidy_result$sets$status, "ok")
  expect_identical(independent_result$tidy_result$sets$test_dimension, 4L)
  expect_identical(independent_result$tidy_result$sets$df2, 1)
  expect_identical(paired_result$tidy_result$sets$status, "ok")
  expect_identical(paired_result$tidy_result$sets$test_dimension, 4L)
  expect_identical(paired_result$tidy_result$sets$df2, 1L)
})

test_that("gly_set_test supports automatic and custom singleton sets", {
  exp <- make_correlated_set_exp()
  automatic <- gly_set_test(exp, threshold = 1)
  custom <- gly_set_test(exp, list(singleton = "A"))

  expect_identical(
    unname(lengths(automatic$raw_result$definitions)),
    rep(1L, 4)
  )
  expect_identical(automatic$tidy_result$sets$test_dimension, rep(1L, 4))
  expect_identical(
    automatic$tidy_result$members$correlation_summary,
    rep(NA_real_, 4)
  )
  expect_identical(custom$tidy_result$sets$status, "ok")
  expect_identical(custom$tidy_result$sets$test_dimension, 1L)
  expect_identical(custom$tidy_result$members$correlation_summary, NA_real_)
})

test_that("gly_set_test retains failed sets and accepts no usable variables", {
  exp <- make_hotelling_exp()
  failed <- gly_set_test(exp, list(all_missing = c("missing1", "missing2")))

  paired <- make_hotelling_exp(paired = TRUE)
  expression <- SummarizedExperiment::assay(paired)
  n_pairs <- ncol(expression) / 2
  expression <- rbind(
    expression,
    zero_difference = rep(expression["A", seq_len(n_pairs)], 2)
  )
  paired <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expression),
    colData = SummarizedExperiment::colData(paired)
  )
  rank_zero <- gly_set_test(
    paired,
    list(rank_zero = "zero_difference"),
    subject_col = "subject"
  )

  unusable <- make_correlated_set_exp()[c("constant", "missing"), ]
  empty <- gly_set_test(unusable)

  expect_identical(failed$tidy_result$sets$status, "failed")
  expect_identical(failed$tidy_result$sets$statistic, NA_real_)
  expect_identical(failed$tidy_result$sets$p_adj, NA_real_)
  expect_identical(
    failed$tidy_result$members$marginal_estimate,
    rep(NA_real_, 2)
  )
  expect_identical(rank_zero$tidy_result$sets$status, "failed")
  expect_identical(rank_zero$tidy_result$sets$test_dimension, 0L)
  expect_match(
    rank_zero$tidy_result$sets$failure_reason,
    "no estimable dimensions",
    fixed = TRUE
  )
  expect_equal(nrow(empty$tidy_result$sets), 0)
  expect_equal(nrow(empty$tidy_result$members), 0)
  expect_length(empty$raw_result$tests, 0)
})

test_that("set construction and testing validate metadata contracts", {
  exp <- make_correlated_set_exp()

  expect_snapshot(
    error = TRUE,
    gly_set_test(exp, list(empty = character()))
  )
  expect_snapshot(
    error = TRUE,
    gly_set_test(exp, within = "unknown")
  )
  expect_snapshot(
    error = TRUE,
    gly_set_test(exp, list(signal = c("A", "unknown")))
  )

  paired <- make_hotelling_exp(paired = TRUE)
  SummarizedExperiment::colData(paired)$subject[2] <-
    SummarizedExperiment::colData(paired)$subject[1]
  expect_snapshot(
    error = TRUE,
    gly_set_test(
      paired,
      list(signal = c("A", "B")),
      subject_col = "subject"
    )
  )
})
