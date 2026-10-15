make_log_transform_set_exp <- function() {
  set.seed(42)
  mat <- matrix(exp(rnorm(72)), nrow = 3,
    dimnames = list(c("A", "B", "C"), paste0("S", 1:24)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = mat),
    colData = S4Vectors::DataFrame(group = factor(rep(c("a", "b"), each = 12)))
  )
}

test_that("set tests use the selected scale for construction and inference", {
  exp <- make_log_transform_set_exp()
  transformed <- exp
  SummarizedExperiment::assay(transformed) <- log2(
    SummarizedExperiment::assay(exp) + 1e-6
  )
  for (sets in list(NULL, list(pair = c("A", "B")))) {
    default <- gly_set_test(exp, sets = sets)
    explicit <- gly_set_test(exp, sets = sets, log_transform = TRUE)
    supplied <- gly_set_test(transformed, sets = sets, log_transform = FALSE)
    expect_equal(default, explicit)
    expect_equal(default, supplied)
  }
  raw <- gly_set_test(exp, sets = list(single = "A"), log_transform = FALSE)
  mat <- SummarizedExperiment::assay(exp)
  groups <- SummarizedExperiment::colData(exp)$group
  expected <- mean(mat["A", groups == levels(groups)[2]]) -
    mean(mat["A", groups == levels(groups)[1]])
  expect_equal(raw$tidy_result$members$marginal_estimate, expected)
})

test_that("fold change preserves arithmetic means and supports log-scale means", {
  mat <- rbind(A = c(1, 4, 2, 16, 8, 32), B = c(0.1, 0.2, 0.4, 0.8, 0.3, 0.6))
  colnames(mat) <- paste0("S", seq_len(ncol(mat)))
  groups <- factor(rep(c("a", "b", "c"), each = 2))
  for (n in c(4, 6)) {
    x <- mat[, seq_len(n), drop = FALSE]
    g <- droplevels(groups[seq_len(n)])
    exp <- SummarizedExperiment::SummarizedExperiment(
      assays = list(expression = x),
      colData = S4Vectors::DataFrame(group = g)
    )
    default <- gly_fold_change(exp, add_info = FALSE)
    expect_equal(default, gly_fold_change(exp, add_info = FALSE, log_transform = FALSE))
    logged <- gly_fold_change(exp, add_info = FALSE, log_transform = TRUE)$tidy_result
    pairs <- combn(levels(g), 2, simplify = FALSE)
    expected_raw <- unlist(lapply(pairs, function(pair) {
      log2(rowMeans(x[, g == pair[2], drop = FALSE]) /
        rowMeans(x[, g == pair[1], drop = FALSE]))
    }), use.names = FALSE)
    z <- log2(x + 1e-6)
    expected_log <- unlist(lapply(pairs, function(pair) {
      rowMeans(z[, g == pair[2], drop = FALSE]) -
        rowMeans(z[, g == pair[1], drop = FALSE])
    }), use.names = FALSE)
    expect_equal(unname(default$tidy_result$log2fc), expected_raw)
    expect_equal(unname(logged$log2fc), expected_log)
    expect_true(all(is.finite(logged$log2fc)))
  }
})

test_that("set and fold-change flags require one non-missing logical", {
  exp <- make_log_transform_set_exp()
  for (invalid in list(NA, NULL, logical(), c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(gly_set_test(exp, log_transform = invalid), "log_transform")
    expect_error(gly_fold_change(exp, log_transform = invalid), "log_transform")
  }
})


test_that("fold-change pseudocount handles zeros without altering raw defaults", {
  mat <- matrix(c(0, 0, 2, 8), nrow = 1,
    dimnames = list("zero", paste0("S", 1:4)))
  groups <- factor(c("a", "a", "b", "b"))
  expect_equal(unname(.analyze_fold_change(mat, groups)$log2fc), Inf)
  expected <- mean(log2(c(2, 8) + 1e-6)) - log2(1e-6)
  expect_equal(unname(.analyze_fold_change(mat, groups, TRUE)$log2fc), expected)
})

test_that("paired set testing skips an already-applied transform", {
  exp <- make_log_transform_set_exp()
  SummarizedExperiment::colData(exp)$subject <- rep(seq_len(12), 2)
  transformed <- exp
  SummarizedExperiment::assay(transformed) <- log2(
    SummarizedExperiment::assay(exp) + 1e-6
  )
  for (sets in list(NULL, list(pair = c("A", "B")))) {
    expect_equal(
      gly_set_test(exp, sets = sets, subject_col = "subject"),
      gly_set_test(transformed, sets = sets, subject_col = "subject", log_transform = FALSE)
    )
  }
})
