test_that("gly_wgcna works correctly", {
  skip_if_not_installed("WGCNA")

  result <- suppressWarnings(suppressMessages(gly_wgcna(test_gp_exp)))

  expect_s3_class(result, c("glystats_wgcna_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("modules", "eigenvalues"))
  expect_true(tibble::is_tibble(result$modules))
  expect_true(tibble::is_tibble(result$eigenvalues))
})

test_that("gly_wgcna_ works correctly", {
  skip_if_not_installed("WGCNA")

  expr_mat <- glyexp::get_expr_mat(test_gp_exp)
  result <- suppressWarnings(suppressMessages(gly_wgcna_(expr_mat)))

  expect_s3_class(result, c("glystats_wgcna_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("modules", "eigenvalues"))
  expect_true(tibble::is_tibble(result$modules))
  expect_true(tibble::is_tibble(result$eigenvalues))
})