#' Principal Component Analysis (PCA)
#'
#' Perform principal component analysis on the expression data.
#' The function uses `prcomp()` to perform PCA and `broom::tidy()` to tidy the results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param center A logical indicating whether to center the data. Default is TRUE.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param ... Additional arguments passed to `prcomp()`.
#'
#' @return A list containing three tibbles:
#'  - `samples`: PCA scores for each sample
#'  - `variables`: PCA loadings for each variable
#'  - `eigenvalues`: PCA eigenvalues
#' 
#' @examples
#' \dontrun{
#' pca_res <- gly_pca(exp)
#' }
#'
#' @export
gly_pca <- function(exp, center = TRUE, scale = TRUE, ...) {
  mat <- log(t(exp$expr_mat) + 1)
  prcomp_res <- stats::prcomp(mat, center = center, scale = scale, ...)
  res <- list(
    "samples" = broom::tidy(prcomp_res, matrix = "samples"),
    "variables" = broom::tidy(prcomp_res, matrix = "variables"),
    "eigenvalues" = broom::tidy(prcomp_res, matrix = "eigenvalues")
  )
  structure(res, class = "gly_pca")
}
