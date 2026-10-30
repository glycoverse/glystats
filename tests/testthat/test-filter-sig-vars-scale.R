raw_limma_filter_experiment <- function(n_groups = 2L) {
  groups <- factor(rep(LETTERS[seq_len(n_groups)], each = 6))
  values <- outer(c(100, 100, 200, 300), rep(1, length(groups))) +
    outer(c(1, 0.8, 1.5, 2.5), rep(c(-2, -1, 0, 0, 1, 2), n_groups))
  values[1, groups == "B"] <- values[1, groups == "B"] + 10
  dimnames(values) <- list(
    c("signal", "control1", "control2", "control3"),
    paste0("s", seq_along(groups))
  )
  glyexp::GlycoproteomicSE(
    values,
    rowData = S4Vectors::DataFrame(
      protein = paste0("P", seq_len(nrow(values))),
      protein_site = seq_len(nrow(values)),
      glycan_composition = glyrepr::as_glycan_composition(rep(
        "H5N2",
        nrow(values)
      )),
      row.names = rownames(values)
    ),
    colData = S4Vectors::DataFrame(
      group = groups,
      row.names = colnames(values)
    ),
    metadata = list(
      glycan_type = "N",
      exp_type = "glycoproteomics",
      scale = "raw"
    )
  )
}

test_that("raw-scale limma filtering keeps significant glycoproteomics variables", {
  skip_if_not_installed("limma")
  glyco <- raw_limma_filter_experiment()
  for (exp in list(glyco, methods::as(glyco, "SummarizedExperiment"))) {
    result <- suppressMessages(gly_limma(
      exp,
      log_transform = FALSE,
      add_info = FALSE
    ))
    original <- exp
    expect_identical(result$preprocessing$analysis_scale, "raw")
    expect_true(all(is.na(result$tidy_result$log2fc)))
    expect_identical(
      result$tidy_result$variable[result$tidy_result$p_adj < 0.05],
      "signal"
    )
    expect_identical(rownames(filter_sig_vars(exp, result)), "signal")
    expect_identical(
      rownames(filter_sig_vars(exp, result, fc_cutoff = NULL)),
      "signal"
    )
    expect_identical(
      rownames(filter_sig_vars(
        exp,
        result,
        p_adj_cutoff = NULL,
        p_val_cutoff = 0.05
      )),
      "signal"
    )
    expect_error(
      filter_sig_vars(exp, result, fc_cutoff = 2),
      "raw-scale",
      class = "rlang_error"
    )
    expect_error(
      filter_sig_vars(exp, result, p_adj_cutoff = NULL, fc_cutoff = 2),
      "raw-scale",
      class = "rlang_error"
    )
    expect_identical(exp, original)
  }
})

test_that("raw-scale limma filtering respects comparison and p-value selection", {
  skip_if_not_installed("limma")
  exp <- raw_limma_filter_experiment(3L)
  result <- suppressMessages(gly_limma(
    exp,
    log_transform = FALSE,
    add_info = FALSE
  ))
  expect_identical(rownames(filter_sig_vars(exp, result)), "signal")
  expect_identical(
    rownames(filter_sig_vars(exp, result, comparison = "A_vs_B")),
    "signal"
  )
  expect_identical(
    rownames(filter_sig_vars(exp, result, comparison = "B_vs_C")),
    "signal"
  )
  expect_length(
    rownames(filter_sig_vars(exp, result, comparison = "A_vs_C")),
    0
  )
  unadjusted <- suppressMessages(gly_limma(
    exp,
    log_transform = FALSE,
    p_adj_method = NULL,
    add_info = FALSE
  ))
  expect_identical(
    rownames(filter_sig_vars(
      exp,
      unadjusted,
      p_adj_cutoff = NULL,
      p_val_cutoff = 0.05
    )),
    "signal"
  )
})

test_that("logged and legacy limma results retain fold-change filtering", {
  skip_if_not_installed("limma")
  exp <- raw_limma_filter_experiment()
  result <- suppressMessages(gly_limma(exp, add_info = FALSE))
  expect_identical(result$preprocessing$analysis_scale, "log2")
  expect_length(rownames(filter_sig_vars(exp, result)), 0)
  expect_identical(
    rownames(filter_sig_vars(exp, result, fc_cutoff = 1)),
    "signal"
  )
  result$preprocessing <- NULL
  expect_length(rownames(filter_sig_vars(exp, result)), 0)
  expect_identical(
    rownames(filter_sig_vars(exp, result, fc_cutoff = 1)),
    "signal"
  )
})
