#' Principal Component Analysis (PCA)
#'
#' Perform principal component analysis on the expression data.
#' The function uses `prcomp()` to perform PCA and `broom::tidy()` to tidy the results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param center A logical indicating whether to center the data. Default is TRUE.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the PCA results are returned.
#'  Only applicable to `gly_pca()`.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw prcomp object.
#' @param ... Additional arguments passed to `prcomp()`.
#'
#' @section Required packages:
#' This function only uses base R packages and does not require additional dependencies.
#'
#' @details
#' The function performs log transformation on the expression data (log(x + 1)) before
#' PCA analysis.
#'
#' `gly_pca()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_pca_()` is the underlying API that works with matrices directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @return A list containing three tibbles (when return_raw = FALSE):
#'  - `samples`: PCA scores for each sample
#'  - `variables`: PCA loadings for each variable
#'  - `eigenvalues`: PCA eigenvalues
#' When return_raw = TRUE, returns the raw prcomp object.
#' @seealso [stats::prcomp()]
#' @export
gly_pca <- function(exp, center = TRUE, scale = TRUE, add_info = TRUE, return_raw = FALSE, ...) {
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)

  # Call the underlying API
  result <- gly_pca_(expr_mat, center, scale, return_raw, ...)

  # If raw results requested, return directly (no add_info processing needed)
  if (return_raw) {
    return(result)
  }

  # Process results with add_info logic
  .process_results_add_info(result, exp, add_info)
}

#' @rdname gly_pca
#' @export
gly_pca_ <- function(expr_mat, center = TRUE, scale = TRUE, return_raw = FALSE, ...) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_logical(center, len = 1)
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

  # Prepare data for PCA (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)
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
    samples_tbl <- dplyr::rename(samples_tbl, c(sample = "row"))
  }
  if ("column" %in% colnames(variables_tbl)) {
    variables_tbl <- dplyr::rename(variables_tbl, c(variable = "column"))
  }

  # Create result list
  result <- list(
    "samples" = samples_tbl,
    "variables" = variables_tbl,
    "eigenvalues" = eigenvalues_tbl
  )

  # Add S3 class
  structure(result, class = c("glystats_pca_res", "glystats_res"))
}
