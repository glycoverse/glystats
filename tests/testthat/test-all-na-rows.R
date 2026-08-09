all_na_row_experiment <- function(all_rows = FALSE) {
  set.seed(20260809)
  n_samples <- 50L
  n_variables <- 8L
  groups <- factor(rep(c("A", "B"), each = n_samples / 2))

  expr_mat <- matrix(
    stats::rnorm(n_variables * n_samples, mean = 5, sd = 0.25),
    nrow = n_variables,
    dimnames = list(
      paste0("V", seq_len(n_variables)),
      paste0("S", seq_len(n_samples))
    )
  )
  expr_mat[1:4, groups == "B"] <- expr_mat[1:4, groups == "B"] + 3
  expr_mat[5:8, groups == "A"] <- expr_mat[5:8, groups == "A"] + 3
  if (all_rows) {
    expr_mat[,] <- NA_real_
  } else {
    expr_mat["V1", ] <- NA_real_
  }

  sample_info <- S4Vectors::DataFrame(
    group = groups,
    age = seq_len(n_samples),
    time = stats::rexp(n_samples),
    event = rep(c(0, 1), length.out = n_samples),
    row.names = colnames(expr_mat)
  )
  variable_info <- S4Vectors::DataFrame(
    type = rep("feature", n_variables),
    row.names = rownames(expr_mat)
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr_mat),
    rowData = variable_info,
    colData = sample_info
  )
}

expect_na_variable_result <- function(tbl, columns, variable = "V1") {
  variable_result <- tbl[tbl$variable == variable, , drop = FALSE]
  expect_gt(nrow(variable_result), 0)
  expect_equal(
    vapply(variable_result[columns], \(x) all(is.na(x)), logical(1)),
    stats::setNames(rep(TRUE, length(columns)), columns)
  )
}

test_that("row-wise analyses retain all-NA variables as NA results", {
  exp <- all_na_row_experiment()

  ttest_result <- expect_no_error(
    suppressMessages(suppressWarnings(gly_ttest(exp)))
  )
  expect_na_variable_result(
    ttest_result$tidy_result,
    c("estimate", "p_val", "p_adj", "log2fc", "effect_size")
  )
  expect_identical(is.na(ttest_result$raw_result[["V1"]]), TRUE)

  wilcox_result <- expect_no_error(
    suppressMessages(suppressWarnings(gly_wilcox(exp)))
  )
  expect_na_variable_result(
    wilcox_result$tidy_result,
    c("statistic", "p_val", "p_adj", "log2fc", "effect_size")
  )
  expect_identical(is.na(wilcox_result$raw_result[["V1"]]), TRUE)

  anova_result <- expect_no_error(
    suppressMessages(suppressWarnings(gly_anova(exp)))
  )
  expect_na_variable_result(
    anova_result$tidy_result$main_test,
    c("statistic", "p_val", "p_adj", "effect_size", "post_hoc")
  )

  ancova_result <- expect_no_error(
    suppressMessages(suppressWarnings(gly_ancova(exp, "group", "age")))
  )
  expect_na_variable_result(
    ancova_result$tidy_result$main_test,
    c("statistic", "p_val", "p_adj", "post_hoc")
  )

  kruskal_result <- expect_no_error(
    suppressMessages(suppressWarnings(gly_kruskal(exp)))
  )
  expect_na_variable_result(
    kruskal_result$tidy_result$main_test,
    c("statistic", "p_val", "p_adj", "effect_size", "post_hoc")
  )

  fold_change_result <- expect_no_error(
    suppressMessages(gly_fold_change(exp))
  )
  expect_na_variable_result(fold_change_result$tidy_result, "log2fc")

  cox_result <- expect_no_error(
    suppressMessages(suppressWarnings(gly_cox(exp)))
  )
  expect_na_variable_result(
    cox_result$tidy_result,
    c("coefficient", "p_val", "p_adj", "hr")
  )
  expect_identical(is.na(cox_result$raw_result[["V1"]]), TRUE)
})

test_that("model-based analyses retain all-NA variables as NA results", {
  skip_if_not_installed("limma")
  exp <- all_na_row_experiment()

  limma_result <- expect_no_error(suppressMessages(gly_limma(exp)))
  expect_na_variable_result(
    limma_result$tidy_result,
    c("log2fc", "AveExpr", "t", "p_val", "p_adj", "b")
  )

  linear_model_result <- expect_no_error(
    gly_linear_model(exp, ~ group + age)
  )
  expect_na_variable_result(
    linear_model_result$tidy_result,
    c("estimate", "AveExpr", "t", "p_val", "p_adj", "b")
  )
})

test_that("correlation and ROC analyses retain all-NA variables", {
  skip_if_not_installed("Hmisc")
  skip_if_not_installed("pROC")
  exp <- all_na_row_experiment()

  variable_cor <- expect_no_error(
    suppressMessages(suppressWarnings(gly_cor(exp, on = "variable")))
  )
  na_cor <- variable_cor$tidy_result[
    variable_cor$tidy_result$variable1 == "V1" |
      variable_cor$tidy_result$variable2 == "V1",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(na_cor), 7)
  expect_equal(na_cor$cor, rep(NA_real_, 7))
  expect_equal(na_cor$p_val, rep(NA_real_, 7))

  expect_no_error(
    suppressMessages(suppressWarnings(gly_cor(exp, on = "sample")))
  )

  roc_result <- expect_no_error(
    suppressMessages(suppressWarnings(gly_roc(exp)))
  )
  expect_na_variable_result(
    roc_result$tidy_result$auc,
    c("auc", "auc_ci_low", "auc_ci_high")
  )
  expect_identical(is.na(roc_result$raw_result[["V1"]]), TRUE)
})

test_that("unsupervised analyses omit all-NA variables only while fitting", {
  exp <- all_na_row_experiment()

  pca_result <- expect_no_error(gly_pca(exp))
  expect_na_variable_result(pca_result$tidy_result$variables, "value")

  variable_kmeans <- expect_no_error(
    gly_kmeans(exp, on = "variable", centers = 2)
  )
  expect_na_variable_result(variable_kmeans$tidy_result, "cluster")
  expect_no_error(gly_kmeans(exp, on = "sample", centers = 2))

  variable_hclust <- expect_no_error(
    gly_hclust(exp, on = "variable", k_values = 2)
  )
  expect_na_variable_result(
    variable_hclust$tidy_result$clusters,
    "cluster_k2"
  )
  expect_no_error(gly_hclust(exp, on = "sample", k_values = 2))
})

test_that("variable k-means handles one retained variable with scaling", {
  exp <- all_na_row_experiment(all_rows = TRUE)
  SummarizedExperiment::assay(exp)["V1", ] <- seq_len(ncol(exp))

  result <- expect_no_error(
    gly_kmeans(exp, on = "variable", centers = 1)
  )

  expect_s3_class(result$raw_result, "kmeans")
  expect_identical(result$tidy_result$cluster[[1]], 1L)
  expect_equal(unname(result$tidy_result$cluster[-1]), rep(NA_integer_, 7))
})

test_that("embeddings ignore all-NA variables while fitting", {
  exp <- all_na_row_experiment()

  skip_if_not_installed("Rtsne")
  tsne_result <- expect_no_error(
    gly_tsne(exp, perplexity = 3, max_iter = 250)
  )
  expect_equal(nrow(tsne_result$tidy_result), ncol(exp))

  skip_if_not_installed("uwot")
  umap_result <- expect_no_error(gly_umap(exp, n_neighbors = 3))
  expect_equal(nrow(umap_result$tidy_result), ncol(exp))
})

test_that("PLS-DA and OPLS-DA restore all-NA variable outputs", {
  skip_if_not_installed("ropls")
  exp <- all_na_row_experiment()

  plsda_result <- expect_no_error(
    suppressMessages(gly_plsda(exp, ncomp = 2, permI = 0))
  )
  expect_na_variable_result(
    plsda_result$tidy_result$variables,
    setdiff(
      colnames(plsda_result$tidy_result$variables),
      c("variable", "type")
    )
  )
  expect_na_variable_result(plsda_result$tidy_result$vip, "vip")

  oplsda_result <- expect_no_error(
    suppressMessages(gly_oplsda(exp, ortho_i = 1, permI = 0))
  )
  expect_na_variable_result(
    oplsda_result$tidy_result$variables,
    setdiff(
      colnames(oplsda_result$tidy_result$variables),
      c("variable", "type")
    )
  )
  expect_na_variable_result(oplsda_result$tidy_result$vip, "vip")
})

test_that("entirely all-NA experiments return placeholders instead of errors", {
  exp <- all_na_row_experiment(all_rows = TRUE)
  calls <- list(
    ttest = \() gly_ttest(exp),
    wilcox = \() gly_wilcox(exp),
    anova = \() gly_anova(exp),
    ancova = \() gly_ancova(exp, covariate_cols = "age"),
    kruskal = \() gly_kruskal(exp),
    fold_change = \() gly_fold_change(exp),
    cox = \() gly_cox(exp),
    pca = \() gly_pca(exp),
    hclust = \() gly_hclust(exp, k_values = 2),
    kmeans = \() gly_kmeans(exp, centers = 2)
  )

  for (call in calls) {
    expect_no_error(suppressMessages(suppressWarnings(call())))
  }
})

test_that("optional analyses tolerate entirely all-NA experiments", {
  exp <- all_na_row_experiment(all_rows = TRUE)

  skip_if_not_installed("Hmisc")
  expect_no_error(gly_cor(exp, on = "variable"))
  expect_no_error(gly_cor(exp, on = "sample"))

  skip_if_not_installed("limma")
  expect_null(suppressMessages(gly_limma(exp))$raw_result)
  expect_null(gly_linear_model(exp, ~ group + age)$raw_result$fit)

  skip_if_not_installed("pROC")
  expect_no_error(suppressMessages(suppressWarnings(gly_roc(exp))))

  skip_if_not_installed("Rtsne")
  expect_null(gly_tsne(exp, perplexity = 3)$raw_result)

  skip_if_not_installed("uwot")
  expect_null(gly_umap(exp, n_neighbors = 3)$raw_result)

  skip_if_not_installed("ropls")
  expect_null(suppressMessages(gly_plsda(exp))$raw_result)
  expect_null(suppressMessages(gly_oplsda(exp))$raw_result)
})
