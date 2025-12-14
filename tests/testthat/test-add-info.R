test_that("add_info parameter works correctly for functions returning tibbles with variable column", {
  # Use test_gp_exp and filter to 2 groups
  exp_2group <- test_gp_exp |>
    glyexp::filter_obs(group %in% c("C", "H")) |>
    glyexp::slice_sample_var(n = 5)

  # Test gly_fold_change
  result_no_info <- suppressMessages(gly_fold_change(exp_2group, add_info = FALSE))
  result_with_info <- suppressMessages(gly_fold_change(exp_2group, add_info = TRUE))

  expect_equal(
    colnames(result_no_info),
    c("variable", "log2fc")
  )

  expect_equal(
    colnames(result_with_info),
    c("variable", "protein", "gene", "glycan_composition", "glycan_structure", "protein_site", "log2fc")
  )
})

test_that("add_info parameter works correctly for functions returning tibbles with sample column", {
  # Use test_gp_exp
  exp_subset <- test_gp_exp |>
    glyexp::slice_sample_var(n = 10)

  # Test gly_pca
  result_pca_no_info <- suppressMessages(gly_pca(exp_subset, add_info = FALSE))
  result_pca_with_info <- suppressMessages(gly_pca(exp_subset, add_info = TRUE))

  tbl_no_info <- get_tidy_result(result_pca_no_info, "samples")
  tbl_with_info <- get_tidy_result(result_pca_with_info, "samples")

  expect_equal(
    colnames(tbl_no_info),
    c("sample", "PC", "value")
  )

  expect_equal(
    colnames(tbl_with_info),
    c("sample", "group", "PC", "value")
  )
})
