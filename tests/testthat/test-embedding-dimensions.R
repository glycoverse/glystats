embedding_dimension_fixture <- function(all_missing = FALSE) {
  mat <- matrix(seq_len(32), nrow = 4)
  dimnames(mat) <- list(paste0("V", 1:4), paste0("S", 1:8))
  if (all_missing) mat[] <- NA_real_
  SummarizedExperiment::SummarizedExperiment(assays = list(abundance = mat))
}

test_that("t-SNE tidy output contains every fitted dimension", {
  skip_if_not_installed("Rtsne")
  local_mocked_bindings(
    Rtsne = function(X, dims, ...) {
      list(Y = matrix(seq_len(nrow(X) * dims), nrow = nrow(X)))
    },
    .package = "Rtsne"
  )
  for (dims in c(1L, 2L, 3L)) {
    result <- gly_tsne(embedding_dimension_fixture(), dims = dims, perplexity = 1, log_transform = FALSE, add_info = FALSE)
    expect_identical(names(result$tidy_result), c("sample", paste0("tsne", seq_len(dims))))
    expect_equal(unname(as.matrix(result$tidy_result[-1])), unname(result$raw_result$Y))
  }
  unnamed <- embedding_dimension_fixture()
  colnames(unnamed) <- NULL
  result <- gly_tsne(unnamed, dims = 2, perplexity = 1, log_transform = FALSE, add_info = FALSE)
  expect_equal(dim(result$tidy_result), c(8L, 3L))
  expect_identical(names(result$tidy_result), c("sample", "tsne1", "tsne2"))
})

test_that("UMAP tidy output contains every fitted dimension", {
  skip_if_not_installed("uwot")
  local_mocked_bindings(
    umap = function(X, n_components, ...) {
      matrix(seq_len(nrow(X) * n_components), nrow = nrow(X))
    },
    .package = "uwot"
  )
  for (dims in c(1L, 2L, 3L)) {
    result <- gly_umap(embedding_dimension_fixture(), n_components = dims, n_neighbors = 3, log_transform = FALSE, add_info = FALSE)
    expect_identical(names(result$tidy_result), c("sample", paste0("umap", seq_len(dims))))
    expect_equal(unname(as.matrix(result$tidy_result[-1])), unname(result$raw_result))
  }
  unnamed <- embedding_dimension_fixture()
  colnames(unnamed) <- NULL
  result <- gly_umap(unnamed, n_components = 2, n_neighbors = 3, log_transform = FALSE, add_info = FALSE)
  expect_equal(dim(result$tidy_result), c(8L, 3L))
  expect_identical(names(result$tidy_result), c("sample", "umap1", "umap2"))
})

test_that("all-missing embeddings respect their requested dimensions", {
  skip_if_not_installed("Rtsne")
  skip_if_not_installed("uwot")
  for (dims in c(1L, 3L)) {
    tsne <- gly_tsne(embedding_dimension_fixture(TRUE), dims = dims, add_info = FALSE)
    umap <- gly_umap(embedding_dimension_fixture(TRUE), n_components = dims, add_info = FALSE)
    expect_identical(names(tsne$tidy_result), c("sample", paste0("tsne", seq_len(dims))))
    expect_identical(names(umap$tidy_result), c("sample", paste0("umap", seq_len(dims))))
    expect_true(all(is.na(as.matrix(tsne$tidy_result[-1]))))
    expect_true(all(is.na(as.matrix(umap$tidy_result[-1]))))
    expect_equal(nrow(tsne$tidy_result), 8)
    expect_null(tsne$raw_result)
    expect_null(umap$raw_result)
  }
  unnamed <- embedding_dimension_fixture(TRUE)
  colnames(unnamed) <- NULL
  expect_equal(dim(gly_tsne(unnamed, dims = 3, add_info = FALSE)$tidy_result), c(8L, 4L))
  expect_equal(dim(gly_umap(unnamed, n_components = 3, add_info = FALSE)$tidy_result), c(8L, 4L))
})

test_that("actual one-dimensional backends return usable tidy results", {
  skip_if_not_installed("Rtsne")
  skip_if_not_installed("uwot")
  exp <- embedding_dimension_fixture()
  tsne <- withr::with_seed(1, gly_tsne(exp, dims = 1, perplexity = 1, pca = FALSE, max_iter = 50, log_transform = FALSE, add_info = FALSE))
  umap <- withr::with_seed(1, gly_umap(exp, n_components = 1, n_neighbors = 3, init = "random", n_epochs = 10, n_threads = 1, log_transform = FALSE, add_info = FALSE))
  expect_equal(dim(tsne$tidy_result), c(8L, 2L))
  expect_equal(dim(umap$tidy_result), c(8L, 2L))
  expect_true(all(is.finite(tsne$tidy_result$tsne1)))
  expect_true(all(is.finite(umap$tidy_result$umap1)))
})
