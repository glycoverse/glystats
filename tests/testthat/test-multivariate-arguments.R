test_that("PCA and UMAP preserve backend output options", {
  exp <- model_argument_experiment()
  fit <- gly_pca(
    exp,
    retx = FALSE,
    rank. = 2,
    scale. = FALSE,
    log_transform = FALSE
  )
  expected <- stats::prcomp(t(SummarizedExperiment::assay(exp)), rank. = 2)
  expect_null(fit$raw_result$x)
  expect_equal(fit$tidy_result$samples$value, as.vector(t(expected$x)))
  expect_equal(fit$tidy_result$eigenvalues$PC, 1:2)
  expect_equal(
    fit$tidy_result$eigenvalues$percent,
    unname(summary(expected)$importance[2, 1:2])
  )

  skip_if_not_installed("uwot")
  fit <- gly_umap(
    exp,
    ret_model = TRUE,
    n_neighbors = 4,
    n_epochs = 10,
    n_threads = 1,
    init = "random"
  )
  expect_equal(
    as.numeric(as.matrix(fit$tidy_result[c("umap1", "umap2")])),
    as.numeric(fit$raw_result$embedding)
  )
})

test_that("hierarchical clustering routes method and aligns cluster sizes", {
  exp <- model_argument_experiment()
  SummarizedExperiment::assay(exp)[1, ] <- NA_real_
  w <- seq_len(nrow(exp))
  mat <- SummarizedExperiment::assay(exp)[-1, ]
  got <- gly_hclust(
    exp,
    method = "average",
    members = w,
    scale = FALSE,
    log_transform = FALSE
  )$raw_result
  expected <- stats::hclust(
    stats::dist(mat),
    method = "average",
    members = w[-1]
  )
  expect_equal(got$height, expected$height)
  got <- gly_hclust(
    exp,
    method = "manhattan",
    scale = FALSE,
    log_transform = FALSE
  )$raw_result
  expect_equal(
    got$height,
    stats::hclust(stats::dist(mat, method = "manhattan"))$height
  )
})

test_that("model backends reject misspelled arguments for empty and populated data", {
  exp <- model_argument_experiment()
  for (empty in c(FALSE, TRUE)) {
    if (empty) {
      SummarizedExperiment::assay(exp)[] <- NA_real_
    }
    expect_error(gly_kmeans(exp, nstrat = 2), "nstrat")
    expect_error(gly_pca(exp, retxx = FALSE), "retxx")
    if (requireNamespace("Rtsne", quietly = TRUE)) {
      expect_error(gly_tsne(exp, max_itre = 5), "max_itre")
    }
    if (requireNamespace("uwot", quietly = TRUE)) {
      expect_error(gly_umap(exp, n_epohcs = 5), "n_epohcs")
    }
  }
})
