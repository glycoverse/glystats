test_that("analysis functions use experiment inputs", {
  analysis_functions <- grep(
    "^gly_",
    getNamespaceExports("glystats"),
    value = TRUE
  )

  expect_false(any(grepl("_$", analysis_functions)))
  purrr::walk(analysis_functions, function(function_name) {
    analysis_function <- getExportedValue("glystats", function_name)
    expect_identical(names(formals(analysis_function))[1], "exp")
  })
})

test_that("analysis functions reject matrix inputs", {
  expr_mat <- matrix(1:4, nrow = 2)

  expect_error(gly_fold_change(expr_mat), "glyexp_experiment")
  expect_error(gly_pca(expr_mat), "glyexp_experiment")
  expect_error(gly_cox(expr_mat), "glyexp_experiment")
})
