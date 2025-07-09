#' Principal Component Analysis (PCA)
#'
#' Perform principal component analysis on the expression data.
#' The function uses `prcomp()` to perform PCA and `broom::tidy()` to tidy the results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param center A logical indicating whether to center the data. Default is TRUE.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the PCA results are returned.
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
#' @seealso [stats::prcomp()]
#' @export
gly_pca <- function(exp, center = TRUE, scale = TRUE, add_info = TRUE, return_raw = FALSE, ...) {
  checkmate::check_logical(add_info, len = 1)
  checkmate::check_logical(return_raw, len = 1)

  mat <- log(t(exp$expr_mat) + 1)
  prcomp_res <- stats::prcomp(mat, center = center, scale = scale, ...)

  # Return raw results if requested
  if (return_raw) {
    return(prcomp_res)
  }

  # Get tidy results and rename columns to be consistent
  samples_tbl <- broom::tidy(prcomp_res, matrix = "samples")
  variables_tbl <- broom::tidy(prcomp_res, matrix = "variables")
  eigenvalues_tbl <- broom::tidy(prcomp_res, matrix = "eigenvalues")

  # Rename columns to be consistent with other glystats functions
  if ("row" %in% colnames(samples_tbl)) {
    samples_tbl <- dplyr::rename(samples_tbl, sample = .data$row)
  }
  if ("column" %in% colnames(variables_tbl)) {
    variables_tbl <- dplyr::rename(variables_tbl, variable = .data$column)
  }

  res <- list(
    "samples" = samples_tbl,
    "variables" = variables_tbl,
    "eigenvalues" = eigenvalues_tbl
  )

  # Process results with add_info logic
  res <- .process_results_add_info(res, exp, add_info)

  structure(res, class = c("glystats_pca_res", "glystats_res"))
}
