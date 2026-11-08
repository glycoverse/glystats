hclust_argument_fixture <- function() {
  mat <- rbind(
    V1 = c(1, 4, 2, 5),
    V2 = c(2, 8, 3, 7),
    V3 = c(9, 3, 8, 2),
    V4 = c(4, 6, 5, 9),
    V5 = c(7, 1, 6, 3)
  )
  colnames(mat) <- paste0("S", seq_len(ncol(mat)))
  SummarizedExperiment::SummarizedExperiment(assays = list(abundance = mat))
}

test_that("hierarchical clustering routes function-specific arguments", {
  exp <- hclust_argument_fixture()
  mat <- SummarizedExperiment::assay(exp)
  members <- c(2, 3, 1, 4, 2)
  expected <- stats::hclust(
    stats::dist(mat, method = "minkowski", p = 1, diag = TRUE, upper = TRUE),
    method = "average",
    members = members
  )
  result <- gly_hclust(
    exp,
    scale = FALSE,
    log_transform = FALSE,
    add_info = FALSE,
    dist.method = "minkowski",
    hclust.method = "average",
    p = 1,
    diag = TRUE,
    upper = TRUE,
    members = members
  )

  for (field in c("merge", "height", "order", "labels", "method", "dist.method")) {
    expect_equal(result$raw_result[[field]], expected[[field]])
  }
  expect_equal(
    unname(result$tidy_result$clusters$cluster_k2),
    unname(stats::cutree(expected, 2))
  )
})

test_that("hierarchical clustering reports clusters remaining after each merge", {
  exp <- hclust_argument_fixture()
  for (on in c("variable", "sample")) {
    result <- gly_hclust(
      exp,
      on = on,
      scale = FALSE,
      log_transform = FALSE,
      add_info = FALSE
    )
    n <- if (on == "variable") nrow(exp) else ncol(exp)
    expect_equal(result$tidy_result$heights$n_clusters, seq.int(n - 1L, 1L))
    expect_equal(tail(result$tidy_result$heights$n_clusters, 1), 1)
  }
})

test_that("hierarchical clustering does not silently discard unknown arguments", {
  expect_error(
    gly_hclust(
      hclust_argument_fixture(),
      scale = FALSE,
      log_transform = FALSE,
      misspelled_argument = TRUE
    ),
    "unused argument"
  )
})


test_that("merge counts do not depend on optional sample or feature labels", {
  exp <- hclust_argument_fixture()
  dimnames(exp) <- list(NULL, NULL)
  for (on in c("variable", "sample")) {
    result <- gly_hclust(
      exp,
      on = on,
      k_values = NULL,
      scale = FALSE,
      log_transform = FALSE,
      add_info = FALSE
    )
    n <- if (on == "variable") nrow(exp) else ncol(exp)
    expect_identical(result$raw_result$labels, result$identifiers[[on]]$data$id)
    expect_length(result$raw_result$order, n)
    expect_equal(result$tidy_result$heights$n_clusters, seq.int(n - 1L, 1L))
    expect_equal(tail(result$tidy_result$heights$n_clusters, 1), 1)
  }
})
