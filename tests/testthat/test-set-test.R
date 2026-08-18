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
    clustering = "connected",
    p_adj_method = NULL
  )

  complete_sets <- complete$raw_result$set_construction
  connected_sets <- connected$raw_result$set_construction
  expect_s3_class(complete, c("glystats_set_test_res", "glystats_res"))
  expect_identical(complete_sets$sets$set_1, c("A", "B", "A_alias"))
  expect_identical(complete_sets$representatives$set_1, c("A", "B"))
  expect_identical(connected_sets$sets$set_1, c("A", "B", "C", "A_alias"))
  expect_identical(connected_sets$representatives$set_1, c("A", "B", "C"))
  expect_identical(complete_sets$aliases$A, "A_alias")
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
  expect_true(is.na(construction$correlation_matrix["A", "C"]))
  expect_equal(
    construction$excluded_variables,
    tibble::tibble(
      variable = c("constant", "missing"),
      reason = c("zero_variance", "all_missing")
    )
  )
  expect_identical(
    construction$membership$is_alias,
    c(FALSE, FALSE, TRUE)
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

test_that("gly_set_test collapses aliases without hiding members", {
  exp <- make_hotelling_exp()
  result <- gly_set_test(
    exp,
    list(signal = c("A", "A_alias", "B")),
    p_adj_method = NULL
  )

  expect_identical(result$tidy_result$sets$n_variables, 3L)
  expect_identical(result$tidy_result$sets$test_dimension, 2L)
  expect_identical(
    result$raw_result$tests$signal$representatives,
    c("A", "B")
  )
  expect_identical(result$raw_result$tests$signal$aliases$A, "A_alias")
  expect_setequal(
    result$tidy_result$members$variable,
    c("A", "A_alias", "B")
  )
  expect_false("p_adj" %in% colnames(result$tidy_result$sets))
})

test_that("gly_set_test collapses aliases after complete-case selection", {
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
      result$raw_result$tests$signal$representatives,
      c("A", "B")
    )
    expect_identical(
      result$raw_result$tests$signal$aliases$A,
      "A_partial"
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
  expect_identical(result$tidy_result$sets$status, "ok")
  expect_identical(result$tidy_result$sets$test_dimension, 2L)
})

test_that("gly_set_test reports unfit sets without dropping them", {
  exp <- make_hotelling_exp()
  result <- gly_set_test(
    exp,
    list(
      all_missing = c("missing1", "missing2"),
      singular = c("A", "B", "sum_ab")
    )
  )

  expect_identical(result$tidy_result$sets$status, c("failed", "failed"))
  expect_true(all(is.na(result$tidy_result$sets$statistic)))
  expect_true(all(is.na(result$tidy_result$sets$p_adj)))
  expect_true(all(is.na(result$tidy_result$members$marginal_estimate[1:2])))
})

test_that("gly_set_test accepts an empty correlated-set result", {
  exp <- make_correlated_set_exp()
  result <- gly_set_test(exp, threshold = 1)

  expect_equal(nrow(result$tidy_result$sets), 0)
  expect_equal(nrow(result$tidy_result$members), 0)
  expect_length(result$raw_result$tests, 0)
})

test_that("set construction and testing validate metadata contracts", {
  exp <- make_correlated_set_exp()

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
