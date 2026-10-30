make_log_transform_differential_exp <- function(n_groups = 3L) {
  set.seed(371)
  groups <- factor(rep(LETTERS[seq_len(n_groups)], each = 6L))
  covariate <- rep(c(-2, -1, 0, 1, 2, 3), n_groups)
  log_expression <- outer(seq(2, 5, length.out = 8), rep(1, length(groups))) +
    outer(seq(0.4, 1.4, length.out = 8), as.integer(groups) - 1) +
    outer(seq(0.02, 0.12, length.out = 8), covariate) +
    matrix(stats::rnorm(8 * length(groups), sd = 0.12), nrow = 8)
  expression <- 2^log_expression
  dimnames(expression) <- list(
    paste0("V", seq_len(nrow(expression))),
    paste0("S", seq_len(ncol(expression)))
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = expression),
    colData = S4Vectors::DataFrame(
      group = groups,
      subject = factor(rep(seq_len(6), n_groups)),
      covariate = covariate,
      row.names = colnames(expression)
    )
  )
}

test_that("the log-transform helper preserves values when disabled", {
  expression <- matrix(
    c(-4, 0, NA_real_, 2, 8, 16),
    nrow = 2,
    dimnames = list(c("V1", "V2"), c("S1", "S2", "S3"))
  )
  expect_identical(.log_transform_expr_mat(expression, FALSE), expression)
  positive <- abs(expression)
  expect_equal(.log_transform_expr_mat(positive), log2(positive + 1e-6))
  expect_equal(.log_transform_expr_mat(positive, TRUE), log2(positive + 1e-6))
})

test_that("differential functions accept automatic mode or a scalar logical flag", {
  exp <- make_log_transform_differential_exp(2)
  functions <- list(
    gly_ttest,
    gly_wilcox,
    gly_anova,
    gly_ancova,
    gly_kruskal,
    gly_limma,
    gly_linear_model
  )
  invalid <- list(NA, logical(), c(TRUE, FALSE), 0, 1, "TRUE")
  for (fun in functions) {
    expect_null(formals(fun)$log_transform)
    for (flag in invalid) {
      expect_error(fun(exp, log_transform = flag), "log_transform")
    }
  }
  for (flag in c(invalid, list(NULL))) {
    expect_error(.log_transform_expr_mat(matrix(1), flag), "log_transform")
  }
})

test_that("t-tests use the chosen scale and preserve raw fold changes", {
  exp <- make_log_transform_differential_exp(2)
  expression <- SummarizedExperiment::assay(exp)
  groups <- SummarizedExperiment::colData(exp)$group
  default <- gly_ttest(exp, add_info = FALSE)
  logged <- gly_ttest(exp, add_info = FALSE, log_transform = TRUE)
  expect_equal(default$tidy_result, logged$tidy_result)

  for (flag in c(TRUE, FALSE)) {
    result <- gly_ttest(
      exp,
      add_info = FALSE,
      var.equal = TRUE,
      alternative = "greater",
      log_transform = flag
    )
    values <- if (flag) log2(expression["V1", ] + 1e-6) else expression["V1", ]
    ref <- values[groups == "A"]
    test <- values[groups == "B"]
    reference <- stats::t.test(
      test,
      ref,
      var.equal = TRUE,
      alternative = "greater"
    )
    row <- result$tidy_result[result$tidy_result$variable == "V1", ]
    pooled_sd <- sqrt((stats::var(ref) + stats::var(test)) / 2)
    expect_equal(row$p_val, reference$p.value)
    expect_equal(unname(row$statistic), unname(reference$statistic))
    expect_equal(row$effect_size, (mean(test) - mean(ref)) / pooled_sd)
    expect_equal(result$tidy_result$log2fc, logged$tidy_result$log2fc)
  }
})

test_that("paired t-tests and Wilcoxon effects use the chosen scale", {
  exp <- make_log_transform_differential_exp(2)[1, ]
  SummarizedExperiment::assay(exp)[1, ] <- c(
    1,
    10,
    3,
    8,
    5,
    20,
    2,
    8,
    5,
    6,
    8,
    17
  )
  values <- as.numeric(SummarizedExperiment::assay(exp))
  for (flag in c(TRUE, FALSE)) {
    analysis <- if (flag) log2(values + 1e-6) else values
    ref <- analysis[seq_len(6)]
    test <- analysis[6 + seq_len(6)]
    differences <- test - ref
    result <- gly_ttest(
      exp,
      subject_col = "subject",
      add_info = FALSE,
      log_transform = flag
    )
    reference <- stats::t.test(test, ref, paired = TRUE)
    expect_equal(result$tidy_result$p_val, reference$p.value)
    expect_equal(unname(result$tidy_result$estimate), mean(differences))
    expect_equal(
      result$tidy_result$effect_size,
      mean(differences) / stats::sd(differences)
    )

    wilcox <- gly_wilcox(
      exp,
      subject_col = "subject",
      exact = FALSE,
      add_info = FALSE,
      log_transform = flag
    )
    wilcox_reference <- stats::wilcox.test(
      test,
      ref,
      paired = TRUE,
      exact = FALSE
    )
    difference_ranks <- rank(abs(differences))
    expect_equal(wilcox$tidy_result$p_val, wilcox_reference$p.value)
    expect_equal(
      unname(wilcox$tidy_result$statistic),
      unname(wilcox_reference$statistic)
    )
    expect_equal(
      wilcox$tidy_result$effect_size,
      sum(sign(differences) * difference_ranks) / sum(difference_ranks)
    )
  }
  default <- gly_wilcox(exp, exact = FALSE, add_info = FALSE)
  logged <- gly_wilcox(
    exp,
    exact = FALSE,
    add_info = FALSE,
    log_transform = TRUE
  )
  expect_equal(default$tidy_result, logged$tidy_result)
})

test_that("ANOVA and its post-hoc models use the chosen scale", {
  exp <- make_log_transform_differential_exp()
  expression <- SummarizedExperiment::assay(exp)
  group <- SummarizedExperiment::colData(exp)$group
  default <- gly_anova(exp, add_info = FALSE)
  logged <- gly_anova(exp, add_info = FALSE, log_transform = TRUE)
  expect_equal(default$tidy_result, logged$tidy_result)

  for (flag in c(TRUE, FALSE)) {
    result <- gly_anova(exp, add_info = FALSE, log_transform = flag)
    value <- if (flag) log2(expression["V1", ] + 1e-6) else expression["V1", ]
    reference <- stats::aov(value ~ group)
    reference_table <- summary(reference)[[1]]
    row <- result$tidy_result$main_test[
      result$tidy_result$main_test$variable == "V1",
    ]
    expect_equal(row$p_val, reference_table[1, "Pr(>F)"])
    expect_equal(
      row$effect_size,
      reference_table[1, "Sum Sq"] / sum(reference_table[, "Sum Sq"])
    )
    expect_equal(
      unname(result$raw_result$post_hoc_test$V1$group[, "diff"]),
      unname(stats::TukeyHSD(reference)$group[, "diff"])
    )
    posthoc <- result$tidy_result$post_hoc_test
    expect_gt(nrow(posthoc), 0)
    expected_fc <- vapply(
      seq_len(nrow(posthoc)),
      function(i) {
        values <- expression[posthoc$variable[i], ]
        log2(
          mean(values[group == posthoc$test_group[i]]) /
            mean(values[group == posthoc$ref_group[i]])
        )
      },
      numeric(1)
    )
    expect_equal(unname(posthoc$log2fc), expected_fc)
  }
})

test_that("paired ANOVA and Friedman tests can use untransformed values", {
  exp <- make_log_transform_differential_exp()
  expression <- SummarizedExperiment::assay(exp)
  sample_info <- as.data.frame(SummarizedExperiment::colData(exp))
  sample_info$value <- as.numeric(expression["V1", ])
  result <- gly_anova(
    exp,
    subject_col = "subject",
    add_info = FALSE,
    log_transform = FALSE
  )
  reference <- stats::aov(value ~ subject + group, data = sample_info)
  expect_equal(
    unname(stats::model.response(stats::model.frame(
      result$raw_result$main_test$V1
    ))),
    sample_info$value
  )
  expect_equal(
    result$tidy_result$main_test$p_val[1],
    summary(reference)[[1]]["group", "Pr(>F)"]
  )

  # Negative supplied values must not be transformed by the rank-test path.
  SummarizedExperiment::assay(exp) <- expression - max(expression) - 1
  sample_info$value <- as.numeric(SummarizedExperiment::assay(exp)["V1", ])
  friedman <- gly_kruskal(
    exp,
    subject_col = "subject",
    add_info = FALSE,
    log_transform = FALSE
  )
  friedman_reference <- stats::friedman.test(
    value ~ group | subject,
    sample_info
  )
  expect_equal(
    friedman$tidy_result$main_test$p_val[1],
    friedman_reference$p.value
  )
  expect_equal(
    unname(friedman$tidy_result$main_test$statistic[1]),
    unname(friedman_reference$statistic)
  )
})

test_that("ANCOVA forwards log_transform without changing model arguments", {
  skip_if_not_installed("emmeans")
  exp <- make_log_transform_differential_exp()
  sample_info <- as.data.frame(SummarizedExperiment::colData(exp))
  expression <- SummarizedExperiment::assay(exp)
  default <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = "covariate",
    add_info = FALSE
  ))
  logged <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = "covariate",
    add_info = FALSE,
    log_transform = TRUE
  ))
  expect_equal(default$tidy_result, logged$tidy_result)
  for (flag in c(TRUE, FALSE)) {
    result <- suppressMessages(gly_ancova(
      exp,
      covariate_cols = "covariate",
      add_info = FALSE,
      log_transform = flag
    ))
    sample_info$value <- if (flag) {
      log2(expression["V1", ] + 1e-6)
    } else {
      expression["V1", ]
    }
    reference <- stats::aov(value ~ group + covariate, data = sample_info)
    expect_equal(
      unname(stats::coef(result$raw_result$main_test$V1)),
      unname(stats::coef(reference))
    )
    expect_gt(nrow(result$tidy_result$post_hoc_test), 0)
  }
})

test_that("Kruskal-Wallis tests accept raw signed values", {
  skip_if_not_installed("FSA")
  exp <- make_log_transform_differential_exp()
  default <- gly_kruskal(exp, add_info = FALSE)
  logged <- gly_kruskal(exp, add_info = FALSE, log_transform = TRUE)
  expect_equal(default$tidy_result, logged$tidy_result)
  expression <- SummarizedExperiment::assay(exp)
  SummarizedExperiment::assay(exp) <- expression - max(expression) - 1
  result <- gly_kruskal(exp, add_info = FALSE, log_transform = FALSE)
  reference <- stats::kruskal.test(
    SummarizedExperiment::assay(exp)["V1", ],
    SummarizedExperiment::colData(exp)$group
  )
  expect_equal(result$tidy_result$main_test$p_val[1], reference$p.value)
  expect_equal(
    unname(result$tidy_result$main_test$statistic[1]),
    unname(reference$statistic)
  )
})

test_that("limma uses the selected scale for tests and reported contrasts", {
  skip_if_not_installed("limma")
  exp <- make_log_transform_differential_exp()
  expression <- SummarizedExperiment::assay(exp)
  sample_info <- as.data.frame(SummarizedExperiment::colData(exp))
  default <- gly_limma(exp, add_info = FALSE)
  logged <- gly_limma(exp, add_info = FALSE, log_transform = TRUE)
  expect_equal(default$tidy_result, logged$tidy_result)
  for (flag in c(TRUE, FALSE)) {
    result <- gly_limma(
      exp,
      covariate_cols = "covariate",
      contrasts = "A_vs_B",
      add_info = FALSE,
      log_transform = flag
    )
    analysis <- if (flag) log2(expression + 1e-6) else expression
    design <- stats::model.matrix(~ 0 + group + covariate, sample_info)
    colnames(design)[seq_len(3)] <- levels(sample_info$group)
    contrast <- limma::makeContrasts(B - A, levels = design)
    reference <- limma::lmFit(analysis, design) |>
      limma::contrasts.fit(contrast) |>
      limma::eBayes(trend = TRUE)
    reference_table <- limma::topTable(
      reference,
      coef = 1,
      number = Inf,
      sort.by = "none"
    )
    if (flag) {
      expect_equal(result$tidy_result$log2fc, reference_table$logFC)
    } else {
      expect_equal(result$tidy_result$estimate, reference_table$logFC)
      expect_true(all(is.na(result$tidy_result$log2fc)))
    }
    expect_equal(result$tidy_result$p_val, reference_table$P.Value)
    expect_equal(result$tidy_result$AveExpr, reference_table$AveExpr)
  }
})

test_that("formula models support supplied scales and named contrasts", {
  skip_if_not_installed("limma")
  exp <- make_log_transform_differential_exp(2)
  expression <- SummarizedExperiment::assay(exp)
  sample_info <- as.data.frame(SummarizedExperiment::colData(exp))
  default <- gly_linear_model(exp, ~ group + covariate, add_info = FALSE)
  logged <- gly_linear_model(
    exp,
    ~ group + covariate,
    add_info = FALSE,
    log_transform = TRUE
  )
  expect_equal(default$tidy_result, logged$tidy_result)
  for (flag in c(TRUE, FALSE)) {
    result <- gly_linear_model(
      exp,
      ~ group + covariate,
      contrasts = c(B_vs_A = "groupB"),
      add_info = FALSE,
      log_transform = flag
    )
    analysis <- if (flag) log2(expression + 1e-6) else expression
    design <- stats::model.matrix(~ group + covariate, sample_info)
    reference <- limma::lmFit(analysis, design) |>
      limma::eBayes(trend = TRUE)
    reference_table <- limma::topTable(
      reference,
      coef = "groupB",
      number = Inf,
      sort.by = "none"
    )
    expect_equal(unname(result$tidy_result$estimate), reference_table$logFC)
    expect_equal(result$tidy_result$p_val, reference_table$P.Value)
    expect_equal(result$tidy_result$AveExpr, reference_table$AveExpr)
  }
})
