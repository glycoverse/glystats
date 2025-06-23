test_that("gly_pca works", {
  # Note: this integration test only makes sure the function runs,
  # it doesn't promise the result is correct.
  pca_res <- gly_pca(test_gp_exp)
  expect_s3_class(pca_res, c("glystats_pca_res", "glystats_res"))
  expect_type(pca_res, "list")
  expect_setequal(names(pca_res), c("samples", "variables", "eigenvalues"))
})