contrast_experiment <- function() {
  groups <- c("H", "M - dose", "C", "Y")
  subject <- rep(seq_len(10), length(groups))
  group <- factor(rep(groups, each = 10), levels = groups)
  withr::local_seed(19)
  values <- matrix(
    2^(rep(as.integer(group), each = 4) + rnorm(160, sd = 0.16)),
    nrow = 4,
    dimnames = list(paste0("v", 1:4), paste0("s", 1:40))
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = values),
    colData = S4Vectors::DataFrame(
      group = group,
      subject = subject,
      age = rnorm(40)
    )
  )
}

test_that("multigroup tests report the requested comparisons in either direction", {
  exp <- contrast_experiment()
  check_comparisons <- function(analyze, ...) {
    full <- analyze(exp, ...)
    expect_no_warning(
      selected <- analyze(exp, ..., contrasts = c("C_vs_M - dose", "H_vs_Y"))
    )
    main <- selected$tidy_result$main_test
    expect_identical(main$p_adj, full$tidy_result$main_test$p_adj)
    all_post <- full$tidy_result$post_hoc_test
    selected_post <- selected$tidy_result$post_hoc_test
    expect_equal(nrow(selected_post), 2L * nrow(exp))
    reversed <- subset(
      selected_post,
      ref_group == "C" & test_group == "M - dose"
    )
    original <- subset(all_post, ref_group == "M - dose" & test_group == "C")
    expect_equal(nrow(reversed), nrow(exp))
    expect_true(all(is.finite(reversed$log2fc)))
    expect_equal(reversed$log2fc, -original$log2fc)
    expect_equal(reversed$p_adj, original$p_adj)
    expect_identical(
      main$post_hoc,
      ifelse(
        original$p_adj[match(main$.variable, original$.variable)] < 0.05,
        "C_vs_M - dose;H_vs_Y",
        "H_vs_Y"
      )
    )

    expect_identical(
      selected$raw_result$post_hoc_test,
      full$raw_result$post_hoc_test
    )
  }
  if (requireNamespace("emmeans", quietly = TRUE)) {
    check_comparisons(gly_ancova, covariate_cols = "age")
    coded <- gly_ancova(
      exp,
      covariate_cols = "age",
      contrasts = list(group = "contr.sum")
    )
    expect_identical(
      coded$raw_result$main_test[[1]]$contrasts$group,
      "contr.sum"
    )
  }
  check_comparisons(gly_kruskal, subject_col = "subject")
  if (requireNamespace("FSA", quietly = TRUE)) check_comparisons(gly_kruskal)
})

test_that("parametric post-hoc tests use the model fitted with the supplied subset", {
  exp <- contrast_experiment()
  keep <- exp$subject != 1
  result <- gly_anova(exp, subset = keep)
  fit <- result$raw_result$main_test[[1]]
  expect_equal(
    result$raw_result$post_hoc_test[[1]]$group,
    stats::TukeyHSD(fit)$group
  )

  skip_if_not_installed("emmeans")
  result <- gly_ancova(exp, covariate_cols = "age", subset = keep)
  fit <- result$raw_result$main_test[[1]]
  expected <- as.data.frame(emmeans::contrast(
    emmeans::emmeans(fit, "group"),
    "pairwise"
  ))
  expect_equal(result$raw_result$post_hoc_test[[1]]$estimate, expected$estimate)
})

test_that("analyses without group contrasts reject them before fitting", {
  exp <- contrast_experiment()
  two <- exp[, exp$group %in% c("H", "C")]
  for (analyze in list(gly_ttest, gly_wilcox)) {
    expect_error(analyze(two, contrasts = "C_vs_H"), "ref_group")
  }
  for (analyze in list(
    gly_cox,
    gly_cor,
    gly_pca,
    gly_tsne,
    gly_umap,
    gly_hclust,
    gly_kmeans,
    gly_plsda,
    gly_oplsda
  )) {
    expect_error(analyze(exp, contrasts = "H_vs_C"), "contrasts")
  }
  result <- gly_anova(exp)
  expect_error(filter_sig_vars(exp, result, contrasts = "H_vs_C"), "contrasts")
})

test_that("post-hoc probabilities and weighted intervals match the fitted model", {
  skip_if_not_installed("emmeans")
  exp <- contrast_experiment()[1, ]
  w <- rep(c(0.5, 1, 3, 1), length.out = ncol(exp))
  offset <- rep(c(0, 0.25, -0.5), length.out = ncol(exp))
  for (weights in list(NULL, w, replace(w, exp$group == "M - dose", 0))) {
    result <- gly_anova(exp, weights = weights, offset = offset)
    fit <- result$raw_result$main_test[[1]]
    means <- emmeans::emmeans(fit, "group")
    if (!is.null(weights)) {
      means <- means[levels(exp$group) %in% exp$group[weights > 0]]
    }
    pairs <- emmeans::contrast(means, "pairwise")
    raw <- as.data.frame(summary(pairs, adjust = "none"))
    adjusted <- as.data.frame(summary(
      pairs,
      adjust = "tukey",
      infer = c(TRUE, TRUE)
    ))
    post <- result$tidy_result$post_hoc_test
    expect_equal(post$p_val, raw$p.value)
    expect_equal(post$p_adj, adjusted$p.value)
    table <- result$raw_result$post_hoc_test[[1]]$group
    expect_equal(unname(table[, "lwr"]), -adjusted$upper.CL)
    expect_equal(unname(table[, "upr"]), -adjusted$lower.CL)
  }
  result <- gly_ancova(exp, covariate_cols = "age")
  fit <- result$raw_result$main_test[[1]]
  pairs <- emmeans::contrast(emmeans::emmeans(fit, "group"), "pairwise")
  expect_equal(
    result$tidy_result$post_hoc_test$p_val,
    summary(pairs, adjust = "none")$p.value
  )
  expect_equal(
    result$tidy_result$post_hoc_test$p_adj,
    summary(pairs, adjust = "tukey")$p.value
  )
})

test_that("limma compares the original group levels when model names collide", {
  skip_if_not_installed("limma")
  exp <- contrast_experiment()
  levels(exp$group) <- c("M-dose", "M.dose", "C", "Y")
  selected <- gly_limma(exp, contrasts = "M-dose_vs_M.dose")
  post <- selected$tidy_result
  values <- log2(SummarizedExperiment::assay(exp) + 1e-6)
  expected <- rowMeans(values[, exp$group == "M.dose"]) -
    rowMeans(values[, exp$group == "M-dose"])
  expect_equal(post$log2fc, unname(expected[post$.variable]))
  expect_error(gly_limma(exp, contrasts = "C_vs_C"), "different groups")
})
