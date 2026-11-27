test_that("add_info parameter works correctly for functions returning tibbles with variable column", {
  # Use test_gp_exp and filter to 2 groups
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5) |>
    as_test_se()

  # Test gly_fold_change
  result_no_info <- suppressMessages(gly_fold_change(
    exp_2group,
    add_info = FALSE
  ))
  result_with_info <- suppressMessages(gly_fold_change(
    exp_2group,
    add_info = TRUE
  ))

  expect_equal(
    colnames(result_no_info$tidy_result),
    c(".variable", "log2fc")
  )

  expect_equal(
    colnames(result_with_info$tidy_result),
    c(
      ".variable",
      "protein",
      "glycan_composition",
      "protein_site",
      "gene",
      "log2fc"
    )
  )
})

test_that("add_info parameter works correctly for functions returning tibbles with sample column", {
  # Use test_gp_exp
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10) |>
    as_test_se()

  # Test gly_pca
  result_pca_no_info <- suppressMessages(gly_pca(exp_subset, add_info = FALSE))
  result_pca_with_info <- suppressMessages(gly_pca(exp_subset, add_info = TRUE))

  tbl_no_info <- get_tidy_result(result_pca_no_info, "samples")
  tbl_with_info <- get_tidy_result(result_pca_with_info, "samples")

  expect_equal(
    colnames(tbl_no_info),
    c(".sample", "PC", "value")
  )

  expect_equal(
    colnames(tbl_with_info),
    c(".sample", "group", "PC", "value")
  )
})

test_that("limma keeps native estimates alongside colliding feature annotations", {
  skip_if_not_installed("limma")
  exp <- model_argument_experiment()
  SummarizedExperiment::rowData(exp)$estimate <- paste0(
    "label",
    seq_len(nrow(exp))
  )
  for (scale in c("raw", "log", "log10", "logx")) {
    S4Vectors::metadata(exp)$scale <- scale
    expected <- gly_limma(
      exp,
      log_transform = FALSE,
      add_info = FALSE
    )$tidy_result
    got <- gly_limma(exp, log_transform = FALSE)$tidy_result
    expect_equal(
      got$estimate,
      expected$estimate[match(got$.variable, expected$.variable)]
    )
    expect_identical(
      got$estimate_metadata,
      SummarizedExperiment::rowData(exp)$estimate
    )
  }
})
