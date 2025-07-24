test_that("gly_wgcna works correctly", {
  skip_if_not_installed("WGCNA")

  result <- suppressWarnings(suppressMessages(gly_wgcna(test_gp_exp)))

  expect_s3_class(result, c("glystats_wgcna_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))

  # Test tidy_result structure
  expect_type(result$tidy_result, "list")
  expect_setequal(names(result$tidy_result), c("modules", "eigenvalues"))
  expect_true(tibble::is_tibble(result$tidy_result$modules))
  expect_true(tibble::is_tibble(result$tidy_result$eigenvalues))

  # Test raw_result
  expect_type(result$raw_result, "list")
  expect_true("colors" %in% names(result$raw_result))
  expect_true("MEs" %in% names(result$raw_result))
})

test_that("gly_wgcna_ works correctly", {
  skip_if_not_installed("WGCNA")

  expr_mat <- glyexp::get_expr_mat(test_gp_exp)
  result <- suppressWarnings(suppressMessages(gly_wgcna_(expr_mat)))

  expect_s3_class(result, c("glystats_wgcna_res", "glystats_res"))
  expect_type(result, "list")
  expect_setequal(names(result), c("tidy_result", "raw_result"))

  # Test tidy_result structure
  expect_type(result$tidy_result, "list")
  expect_setequal(names(result$tidy_result), c("modules", "eigenvalues"))
  expect_true(tibble::is_tibble(result$tidy_result$modules))
  expect_true(tibble::is_tibble(result$tidy_result$eigenvalues))

  # Test raw_result
  expect_type(result$raw_result, "list")
  expect_true("colors" %in% names(result$raw_result))
  expect_true("MEs" %in% names(result$raw_result))
})