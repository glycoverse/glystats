multivariate_log_fixture <- function() {
  index <- seq_len(30)
  mat <- rbind(
    V1 = (index %% 7 + 1)^2 + index / 10,
    V2 = (index %% 5 + 2)^3 / 3,
    V3 = exp(index / 8) + index %% 3,
    V4 = (index * 3 %% 11) + 5,
    V5 = (index * 7 %% 13)^2 + 2,
    V6 = sqrt(index + 1) + index %% 4
  )
  colnames(mat) <- paste0("S", index)
  sample_info <- S4Vectors::DataFrame(
    group = factor(rep(c("A", "B"), each = 15)),
    time = (index * 7 %% 31) + 1,
    event = as.numeric(index %% 4 != 0),
    row.names = colnames(mat)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = mat),
    colData = sample_info
  )
}

test_that("multivariate log flags preserve defaults and positional arguments", {
  functions <- list(
    pca = list(gly_pca, .analyze_pca),
    cor = list(gly_cor, .analyze_cor),
    cox = list(gly_cox, .analyze_cox),
    roc = list(gly_roc, .analyze_roc),
    hclust = list(gly_hclust, .analyze_hclust),
    kmeans = list(gly_kmeans, .analyze_kmeans),
    oplsda = list(gly_oplsda, .analyze_oplsda),
    plsda = list(gly_plsda, .analyze_plsda),
    tsne = list(gly_tsne, .analyze_tsne),
    umap = list(gly_umap, .analyze_umap)
  )
  invalid_flags <- list(NA, NULL, logical(), c(TRUE, FALSE), 1, "TRUE")

  for (name in names(functions)) {
    for (fn in functions[[name]]) {
      args <- formals(fn)
      expect_identical(tail(names(args), 1), "log_transform")
      expect_identical(args$log_transform, !name %in% c("cox", "roc"))
      if (name != "roc") {
        expect_identical(tail(names(args), 2), c("...", "log_transform"))
      }
      for (flag in invalid_flags) {
        expect_error(do.call(fn, list(log_transform = flag)), "log_transform")
      }
    }
  }
})

test_that("PCA keeps natural-log defaults and can use signed input unchanged", {
  exp <- multivariate_log_fixture()
  mat <- SummarizedExperiment::assay(exp)
  default <- gly_pca(exp, scale = FALSE, add_info = FALSE)
  explicit <- gly_pca(exp, scale = FALSE, log_transform = TRUE, add_info = FALSE)
  expected <- stats::prcomp(t(log(mat + 1)), center = TRUE, scale. = FALSE)
  expect_equal(default$raw_result, expected)
  expect_equal(default$tidy_result, explicit$tidy_result)

  SummarizedExperiment::assay(exp) <- log(mat + 1) - 5
  raw_mat <- SummarizedExperiment::assay(exp)
  for (scale in c(FALSE, TRUE)) {
    expect_no_warning(result <- gly_pca(
      exp, scale = scale, add_info = FALSE, log_transform = FALSE
    ))
    expected <- stats::prcomp(t(raw_mat), center = TRUE, scale. = scale)
    expect_equal(result$raw_result, expected)
  }
})

test_that("correlation respects both log modes and matrix orientations", {
  skip_if_not_installed("Hmisc")
  exp <- multivariate_log_fixture()
  mat <- SummarizedExperiment::assay(exp)

  for (on in c("variable", "sample")) {
    default <- gly_cor(exp, on = on)
    explicit <- gly_cor(exp, on = on, log_transform = TRUE)
    expected_mat <- if (on == "variable") t(mat) else mat
    expected <- Hmisc::rcorr(log2(expected_mat + 1e-6), type = "pearson")
    expect_equal(default$raw_result, expected)
    expect_equal(default$tidy_result, explicit$tidy_result)

    signed_exp <- exp
    SummarizedExperiment::assay(signed_exp) <- mat - 100
    expect_no_warning(result <- gly_cor(
      signed_exp, on = on, log_transform = FALSE
    ))
    expected <- Hmisc::rcorr(expected_mat - 100, type = "pearson")
    expect_equal(result$raw_result, expected)
  }
})

test_that("hierarchical clustering transforms only when requested", {
  exp <- multivariate_log_fixture()
  mat <- SummarizedExperiment::assay(exp)

  for (on in c("variable", "sample")) {
    expected_mat <- if (on == "sample") t(mat) else mat
    default <- suppressMessages(gly_hclust(
      exp, on = on, scale = FALSE, add_info = FALSE
    ))
    explicit <- suppressMessages(gly_hclust(
      exp, on = on, scale = FALSE, add_info = FALSE, log_transform = TRUE
    ))
    expected <- stats::hclust(stats::dist(log2(expected_mat + 1e-6)))
    expect_equal(default$raw_result$merge, expected$merge)
    expect_equal(default$raw_result$height, expected$height)
    expect_equal(default$tidy_result, explicit$tidy_result)

    signed_exp <- exp
    SummarizedExperiment::assay(signed_exp) <- mat - 100
    result <- suppressMessages(gly_hclust(
      signed_exp, on = on, scale = FALSE, add_info = FALSE,
      dist.method = "manhattan", hclust.method = "average",
      log_transform = FALSE
    ))
    expected <- stats::hclust(
      stats::dist(expected_mat - 100, method = "manhattan"), method = "average"
    )
    expect_equal(result$raw_result$merge, expected$merge)
    expect_equal(result$raw_result$height, expected$height)
  }
})

test_that("k-means respects log modes without forwarding the flag", {
  exp <- multivariate_log_fixture()
  mat <- SummarizedExperiment::assay(exp)

  for (on in c("variable", "sample")) {
    expected_mat <- if (on == "sample") t(mat) else mat
    default <- gly_kmeans(
      exp, on = on, centers = 1, scale = FALSE, add_info = FALSE
    )
    explicit <- gly_kmeans(
      exp, on = on, centers = 1, scale = FALSE,
      add_info = FALSE, log_transform = TRUE
    )
    expected <- stats::kmeans(log2(expected_mat + 1e-6), centers = 1)
    expect_equal(default$raw_result$centers, expected$centers)
    expect_equal(default$raw_result$totss, expected$totss)
    expect_equal(default$tidy_result, explicit$tidy_result)

    signed_exp <- exp
    SummarizedExperiment::assay(signed_exp) <- mat - 100
    expect_no_warning(result <- gly_kmeans(
      signed_exp, on = on, centers = 1, scale = FALSE, add_info = FALSE,
      iter.max = 20, log_transform = FALSE
    ))
    expected <- stats::kmeans(expected_mat - 100, centers = 1, iter.max = 20)
    expect_equal(result$raw_result$centers, expected$centers)
    expect_equal(result$raw_result$totss, expected$totss)
  }
})

test_that("Cox models retain raw defaults and offer opt-in base-2 logs", {
  exp <- multivariate_log_fixture()
  set.seed(72)
  SummarizedExperiment::colData(exp)$time <- rexp(ncol(exp))
  SummarizedExperiment::colData(exp)$event <- rep(c(1, 1, 0), length.out = ncol(exp))
  mat <- SummarizedExperiment::assay(exp)
  default <- gly_cox(exp, add_info = FALSE, x = TRUE)
  raw <- gly_cox(exp, add_info = FALSE, x = TRUE, log_transform = FALSE)
  logged <- gly_cox(exp, add_info = FALSE, x = TRUE, log_transform = TRUE)
  expect_equal(default$tidy_result, raw$tidy_result)

  for (i in seq_len(nrow(mat))) {
    expect_equal(unname(default$raw_result[[i]]$x[, 1]), unname(mat[i, ]))
    expect_equal(
      unname(logged$raw_result[[i]]$x[, 1]), unname(log2(mat[i, ] + 1e-6))
    )
  }
})

test_that("ROC retains raw defaults and transforms predictors on request", {
  skip_if_not_installed("pROC")
  exp <- multivariate_log_fixture()
  mat <- SummarizedExperiment::assay(exp)
  default <- suppressMessages(gly_roc(exp, pos_class = "B", add_info = FALSE))
  raw <- suppressMessages(gly_roc(
    exp, pos_class = "B", add_info = FALSE, log_transform = FALSE
  ))
  logged <- suppressMessages(gly_roc(
    exp, pos_class = "B", add_info = FALSE, log_transform = TRUE
  ))
  expect_equal(default$tidy_result, raw$tidy_result)
  expect_equal(default$tidy_result$auc, logged$tidy_result$auc)

  for (i in seq_len(nrow(mat))) {
    expect_equal(unname(default$raw_result[[i]]$predictor), unname(mat[i, ]))
    expect_equal(
      unname(logged$raw_result[[i]]$predictor), unname(log2(mat[i, ] + 1e-6))
    )
  }
})

test_that("t-SNE applies the historical natural log only when enabled", {
  skip_if_not_installed("Rtsne")
  exp <- multivariate_log_fixture()
  mat <- SummarizedExperiment::assay(exp)
  set.seed(2025)
  default <- gly_tsne(
    exp, perplexity = 3, add_info = FALSE, pca = FALSE, max_iter = 250
  )
  set.seed(2025)
  explicit <- gly_tsne(
    exp, perplexity = 3, add_info = FALSE, pca = FALSE, max_iter = 250,
    log_transform = TRUE
  )
  expect_equal(default$tidy_result, explicit$tidy_result)
  set.seed(2025)
  expected <- Rtsne::Rtsne(
    t(log(mat + 1)), dims = 2, perplexity = 3, pca = FALSE,
    max_iter = 250, check_duplicates = FALSE
  )
  expect_equal(default$raw_result$Y, expected$Y)

  SummarizedExperiment::assay(exp) <- log(mat + 1) - 5
  set.seed(2025)
  raw <- gly_tsne(
    exp, perplexity = 3, add_info = FALSE, pca = FALSE, max_iter = 250,
    log_transform = FALSE
  )
  set.seed(2025)
  expected <- Rtsne::Rtsne(
    t(SummarizedExperiment::assay(exp)), dims = 2, perplexity = 3,
    pca = FALSE, max_iter = 250, check_duplicates = FALSE
  )
  expect_equal(raw$raw_result$Y, expected$Y)
})

test_that("UMAP applies the historical natural log only when enabled", {
  skip_if_not_installed("uwot")
  exp <- multivariate_log_fixture()
  mat <- SummarizedExperiment::assay(exp)
  set.seed(2025)
  default <- gly_umap(
    exp, n_neighbors = 4, add_info = FALSE, n_epochs = 10,
    init = "random", n_threads = 1
  )
  set.seed(2025)
  explicit <- gly_umap(
    exp, n_neighbors = 4, add_info = FALSE, n_epochs = 10,
    init = "random", n_threads = 1, log_transform = TRUE
  )
  expect_equal(default$tidy_result, explicit$tidy_result)
  set.seed(2025)
  expected <- uwot::umap(
    t(log(mat + 1)), n_neighbors = 4, n_components = 2, n_epochs = 10,
    init = "random", n_threads = 1, verbose = FALSE
  )
  expect_equal(default$raw_result, expected)

  SummarizedExperiment::assay(exp) <- log(mat + 1) - 5
  set.seed(2025)
  raw <- gly_umap(
    exp, n_neighbors = 4, add_info = FALSE, n_epochs = 10,
    init = "random", n_threads = 1, log_transform = FALSE
  )
  set.seed(2025)
  expected <- uwot::umap(
    t(SummarizedExperiment::assay(exp)), n_neighbors = 4, n_components = 2,
    n_epochs = 10, init = "random", n_threads = 1, verbose = FALSE
  )
  expect_equal(raw$raw_result, expected)
})

test_that("PLS-DA and OPLS-DA can consume prelogged data exactly once", {
  skip_if_not_installed("ropls")
  exp <- exp_topliss_valid()
  logged_exp <- exp
  SummarizedExperiment::assay(logged_exp) <- log(
    SummarizedExperiment::assay(exp) + 1
  )
  cases <- list(
    list(fn = gly_plsda, args = list(ncomp = 1)),
    list(fn = gly_oplsda, args = list(pred_i = 1, ortho_i = 1))
  )
  for (case in cases) {
    args <- c(case$args, list(
      scale = FALSE, add_info = FALSE, permI = 0, crossvalI = 3
    ))
    set.seed(2025)
    default <- suppressMessages(do.call(case$fn, c(list(exp = exp), args)))
    set.seed(2025)
    explicit <- suppressMessages(do.call(case$fn, c(
      list(exp = exp, log_transform = TRUE), args
    )))
    set.seed(2025)
    raw <- suppressMessages(do.call(case$fn, c(
      list(exp = logged_exp, log_transform = FALSE), args
    )))
    expect_equal(default$tidy_result, explicit$tidy_result)
    expect_equal(default$tidy_result, raw$tidy_result)
    expect_equal(default$raw_result@scoreMN, raw$raw_result@scoreMN)
    expect_equal(default$raw_result@loadingMN, raw$raw_result@loadingMN)
  }
})

test_that("disabling logs preserves all-missing variable handling", {
  exp <- multivariate_log_fixture()
  SummarizedExperiment::assay(exp)[1, ] <- NA_real_
  pca <- gly_pca(exp, add_info = FALSE, log_transform = FALSE)
  expect_true(all(is.na(
    pca$tidy_result$variables$value[pca$tidy_result$variables$variable == "V1"]
  )))
  kmeans <- gly_kmeans(exp, centers = 1, add_info = FALSE, log_transform = FALSE)
  expect_true(is.na(kmeans$tidy_result$cluster[1]))
  hclust <- suppressMessages(gly_hclust(
    exp, add_info = FALSE, log_transform = FALSE
  ))
  expect_true(is.na(hclust$tidy_result$clusters$cluster_k2[1]))
})
