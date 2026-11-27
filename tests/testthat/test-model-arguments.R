test_that("sample selection precedes main tests, effect sizes and post-hoc fits", {
  exp <- model_argument_experiment()
  keep <- seq_len(ncol(exp)) %% 4 != 0
  compare <- function(analyze, ...) {
    got <- analyze(exp, ..., subset = keep)
    expected <- analyze(exp[, keep], ...)
    expect_equal(got$tidy_result, expected$tidy_result)
    expect_gt(nrow(got$tidy_result$post_hoc_test), 0L)
  }
  compare(gly_anova)
  compare(gly_anova, subject_col = "subject")
  compare(gly_kruskal, subject_col = "subject")
  if (requireNamespace("FSA", quietly = TRUE)) {
    compare(gly_kruskal)
  }
  if (requireNamespace("emmeans", quietly = TRUE)) {
    compare(gly_ancova, covariate_cols = "age")
  }
})

test_that("paired sample arguments track exclusions and model order", {
  exp <- model_argument_experiment()
  values <- SummarizedExperiment::assay(exp)
  values[1, 1] <- NA_real_
  SummarizedExperiment::assay(exp) <- values
  weights <- seq_len(ncol(exp))
  order <- c(seq(32, 2, -2), seq(1, 31, 2))
  exp <- exp[, order]
  keep <- exp$subject != 8
  r <- gly_anova(
    exp,
    subset = keep,
    weights = weights[order],
    subject_col = "subject"
  )
  fit <- r$raw_result$main_test[[1]]
  expected_weights <- weights[order][match(
    paste(fit$model$subject, fit$model$group),
    paste(exp$subject, exp$group)
  )]
  expect_equal(unname(stats::weights(fit)), unname(expected_weights))
  for (analyze in list(gly_ttest, gly_wilcox)) {
    aligned <- analyze(exp[, order(order)], subject_col = "subject")
    positional <- analyze(exp[, order(order)], paired = TRUE)
    expect_equal(positional$tidy_result, aligned$tidy_result)
    expect_error(
      analyze(exp, subject_col = "subject", paired = FALSE),
      "conflicts"
    )
  }
  expect_error(gly_ttest(exp, var.euqal = TRUE), "var.euqal")
})

test_that("weight filtering in both limma wrappers preserves observation and sample weights", {
  skip_if_not_installed("limma")
  exp <- model_argument_experiment()[, 1:8]
  exp$group <- factor(rep(c("A", "B"), 4))
  SummarizedExperiment::assay(exp)[c(1, 4), ] <- NA_real_
  keep <- !seq_len(nrow(exp)) %in% c(1, 4)
  w <- matrix(seq_len(length(SummarizedExperiment::assay(exp))), nrow(exp))
  compare <- function(analyze, ...) {
    got <- analyze(exp, ..., weights = w)$tidy_result
    expected <- analyze(exp[keep, ], ..., weights = w[keep, ])$tidy_result
    expect_equal(
      got$p_val[match(expected$.variable, got$.variable)],
      expected$p_val
    )
    sample_weights <- seq_len(ncol(exp))
    got <- analyze(exp, ..., weights = sample_weights)$tidy_result
    expected <- analyze(
      exp[keep, ],
      ...,
      weights = rbind(sample_weights)
    )$tidy_result
    expect_equal(
      got$p_val[match(expected$.variable, got$.variable)],
      expected$p_val
    )
  }
  compare(gly_limma)
  compare(gly_linear_model, ~group)
})

test_that("Cox weights and backend controls reach the fitted model", {
  exp <- model_argument_experiment()
  w <- seq_len(ncol(exp))
  got <- gly_cox(
    exp,
    time_col = "time",
    event_col = "event",
    weights = w,
    ties = "breslow",
    iter.max = 30,
    log_transform = FALSE
  )
  x <- SummarizedExperiment::assay(exp)[1, ]
  expected <- survival::coxph(
    survival::Surv(exp$time, exp$event) ~ x,
    weights = w,
    ties = "breslow",
    iter.max = 30
  )
  fit <- got$raw_result[[1]]
  expect_equal(stats::coef(fit), stats::coef(expected))
  expect_equal(fit$weights, w)
  expect_error(
    gly_cox(exp, time_col = "time", event_col = "event", tiees = "exact"),
    "tiees"
  )
})

test_that("filter thresholds are taken from the caller when annotations share their names", {
  exp <- model_argument_experiment()
  SummarizedExperiment::rowData(exp)$p_adj_cutoff <- 0
  SummarizedExperiment::rowData(exp)$p_val_cutoff <- 0
  SummarizedExperiment::rowData(exp)$fc_cutoff <- Inf
  fit <- gly_ttest(exp)
  expect_equal(
    nrow(filter_sig_vars(exp, fit, p_adj_cutoff = 1, fc_cutoff = 1)),
    nrow(exp)
  )
  expect_equal(
    nrow(filter_sig_vars(exp, fit, p_adj_cutoff = NULL, p_val_cutoff = 1)),
    nrow(exp)
  )
})

test_that("discarding raw QR storage preserves ANOVA statistics", {
  exp <- model_argument_experiment()
  for (analyze in list(gly_anova, function(...) {
    gly_ancova(..., covariate_cols = "age")
  })) {
    if (!identical(analyze, gly_anova)) {
      skip_if_not_installed("emmeans")
    }
    got <- analyze(exp, qr = FALSE)
    expect_null(got$raw_result$main_test[[1]]$qr)
    expect_equal(got$tidy_result, analyze(exp)$tidy_result)
  }
})

test_that("correlation backend names agree with the wrapper and preserve auxiliary data", {
  skip_if_not_installed("Hmisc")
  exp <- model_argument_experiment()
  expected <- gly_cor(exp, method = "spearman")
  got <- gly_cor(exp, type = "spearman")
  expect_equal(got$raw_result$r, expected$raw_result$r)
  expect_error(
    gly_cor(exp, method = "pearson", type = "spearman"),
    "different correlations"
  )
  y <- cbind(a = seq_len(ncol(exp)), b = rev(seq_len(ncol(exp))))
  SummarizedExperiment::assay(exp)[] <- NA_real_
  expect_no_warning(fit <- gly_cor(exp, y = y))
  expected <- Hmisc::rcorr(y)
  for (field in c("r", "n", "P")) {
    expect_equal(
      fit$raw_result[[field]][colnames(y), colnames(y)],
      expected[[field]]
    )
  }
  expect_equal(fit$raw_result$n[1, "a"], 0L)
})
