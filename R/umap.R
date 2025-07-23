#' Uniform Manifold Approximation and Projection (UMAP)
#'
#' Perform UMAP dimensionality reduction on the expression data.
#' The function uses `uwot::umap()` to perform UMAP analysis.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param n_neighbors Number of neighbors to consider for each point. Default is 15.
#' @param n_components Number of output dimensions. Default is 2.
#' @param min_dist Minimum distance between embedded points. Default is 0.1.
#' @param spread Controls how tightly the embedding is packed. Default is 1.0.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'  will be added to the result tibble. If FALSE, only the UMAP coordinates are returned.
#'  Only applicable to `gly_umap()`.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw umap result matrix.
#' @param ... Additional arguments passed to `uwot::umap()`.
#'
#' @section Required packages:
#' This function requires the `uwot` package to be installed for UMAP analysis.
#'
#' @details
#' `gly_umap()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_umap_()` is the underlying API that works with matrices directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @return A tibble with UMAP coordinates (umap1, umap2) when return_raw = FALSE,
#'   or raw umap result matrix when return_raw = TRUE.
#' @seealso [uwot::umap()]
#' @export
gly_umap <- function(
  exp,
  n_neighbors = 15,
  n_components = 2,
  min_dist = 0.1,
  spread = 1,
  add_info = TRUE,
  return_raw = FALSE,
  ...
) {
  .check_pkg_available("uwot")
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_flag(add_info)

  expr_mat <- glyexp::get_expr_mat(exp)
  result <- gly_umap_(expr_mat, n_components, n_neighbors, min_dist, spread, return_raw, ...)

  .process_results_add_info(result, exp, add_info)
}

#' @rdname gly_umap
#' @export
gly_umap_ <- function(
  expr_mat,
  n_components = 2,
  n_neighbors = 15,
  min_dist = 0.1,
  spread = 1,
  return_raw = FALSE,
  ...
) {
  .check_pkg_available("uwot")

  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_int(n_components, lower = 1)
  checkmate::assert_int(n_neighbors, lower = 1)
  checkmate::assert_number(min_dist, lower = 0, upper = 1)
  checkmate::assert_number(spread, lower = 0)
  checkmate::assert_flag(return_raw)

  # Prepare data (samples as rows, variables as columns)
  mat <- t(expr_mat)

  # Apply log transformation
  mat <- log(mat + 1)

  # Perform UMAP
  umap_res <- uwot::umap(
    X = mat,
    n_components = n_components,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    spread = spread,
    verbose = FALSE,
    ...
  )

  # Return raw results if requested
  if (return_raw) {
    return(umap_res)
  }

  # Create result tibble
  result <- tibble::tibble(
    sample = rownames(mat),
    umap1 = umap_res[, 1],
    umap2 = umap_res[, 2]
  )

  # Handle more than 2 components
  if (n_components > 2) {
    for (i in 3:n_components) {
      result[[paste0("umap", i)]] <- umap_res[, i]
    }
  }

  # Add S3 class
  structure(result, class = c("glystats_umap_res", "glystats_res", class(result)))
}