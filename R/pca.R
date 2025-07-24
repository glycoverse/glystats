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
#' @return A list containing:
#'  - `tidy_result`: A list of tibbles with PCA results:
#'    - `samples`: PCA scores for each sample containing the following columns:
#'      - `sample`: Sample name
#'      - `PC`: Principal component name (PC1, PC2, etc.)
#'      - `value`: Score value for the principal component
#'    - `variables`: PCA loadings for each variable containing the following columns:
#'      - `variable`: Variable name
#'      - `PC`: Principal component name (PC1, PC2, etc.)
#'      - `value`: Loading value for the principal component
#'    - `eigenvalues`: PCA eigenvalues containing the following columns:
#'      - `PC`: Principal component name (PC1, PC2, etc.)
#'      - `std.dev`: Standard deviation
#'      - `percent`: Percentage of variance explained
#'      - `cumulative`: Cumulative percentage of variance explained
#'  - `raw_result`: The raw prcomp object from `stats::prcomp()`
#' @seealso [stats::prcomp()]
#' @export
gly_pca <- function(exp, center = TRUE, scale = TRUE, add_info = TRUE, ...) {
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)

  # Call the underlying API
  result <- gly_pca_(expr_mat, center, scale, ...)
  result$tidy_result <- .process_results_add_info(result$tidy_result, exp, add_info)
  result
}

#' @rdname gly_pca
#' @export
gly_pca_ <- function(expr_mat, center = TRUE, scale = TRUE, ...) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_logical(center, len = 1)
  checkmate::assert_logical(scale, len = 1)

  # Prepare data for PCA (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)
  prcomp_res <- stats::prcomp(mat, center = center, scale = scale, ...)

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

  # Create tidy result list
  tidy_result <- list(
    "samples" = samples_tbl,
    "variables" = variables_tbl,
    "eigenvalues" = eigenvalues_tbl
  )

  # Add S3 class to tidy_result
  tidy_result <- structure(tidy_result, class = c("glystats_pca_res", "glystats_res"))

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = prcomp_res
    ),
    class = c("glystats_pca_res", "glystats_res")
  )
}
