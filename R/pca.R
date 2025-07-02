#' Principal Component Analysis (PCA)
#'
#' Perform principal component analysis on the expression data.
#' The function uses `prcomp()` to perform PCA and `broom::tidy()` to tidy the results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param center A logical indicating whether to center the data. Default is TRUE.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw prcomp object.
#' @param ... Additional arguments passed to `prcomp()`.
#'
#' @section Required packages:
#' This function only uses base R packages and does not require additional dependencies.
#'
#' @return A list containing three tibbles (when return_raw = FALSE):
#'  - `samples`: PCA scores for each sample
#'  - `variables`: PCA loadings for each variable
#'  - `eigenvalues`: PCA eigenvalues
#' When return_raw = TRUE, returns the raw prcomp object.
#'
#' @export
gly_pca <- function(exp, center = TRUE, scale = TRUE, return_raw = FALSE, ...) {
  checkmate::check_logical(return_raw, len = 1)
  
  mat <- log(t(exp$expr_mat) + 1)
  prcomp_res <- stats::prcomp(mat, center = center, scale = scale, ...)
  
  # Return raw results if requested
  if (return_raw) {
    return(prcomp_res)
  }
  
  res <- list(
    "samples" = broom::tidy(prcomp_res, matrix = "samples"),
    "variables" = broom::tidy(prcomp_res, matrix = "variables"),
    "eigenvalues" = broom::tidy(prcomp_res, matrix = "eigenvalues")
  )
  structure(res, class = c("glystats_pca_res", "glystats_res"))
}
