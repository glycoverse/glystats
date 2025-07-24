#' Uniform Manifold Approximation and Projection (UMAP)
#'
#' Perform UMAP dimensionality reduction on the expression data.
#' The function uses `uwot::umap()` to perform UMAP analysis.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param n_neighbors Number of neighbors to consider for each point. Default is 15.
#' @param n_components Number of output dimensions. Default is 2.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'  will be added to the result tibble. If FALSE, only the UMAP coordinates are returned.
#'  Only applicable to `gly_umap()`.
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
#' @return A list with two elements:
#' - `tidy_result`: A tibble with UMAP coordinates (umap1, umap2, etc.) and sample names
#' - `raw_result`: The raw UMAP result matrix
#' The list has classes `glystats_umap_res` and `glystats_res`.
#' @seealso [uwot::umap()]
#' @export
gly_umap <- function(
  exp,
  n_neighbors = 15,
  n_components = 2,
  add_info = TRUE,
  ...
) {
  .check_pkg_available("uwot")
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_flag(add_info)

  expr_mat <- glyexp::get_expr_mat(exp)
  result <- gly_umap_(expr_mat, n_neighbors, n_components, ...)

  result$tidy_result <- .process_results_add_info(result$tidy_result, exp, add_info)
  result
}

#' @rdname gly_umap
#' @export
gly_umap_ <- function(
  expr_mat,
  n_neighbors = 15,
  n_components = 2,
  ...
) {
  .check_pkg_available("uwot")

  checkmate::assert_matrix(expr_mat, mode = "numeric")

  # Prepare data (samples as rows, variables as columns)
  mat <- t(expr_mat)

  # Apply log transformation
  mat <- log(mat + 1)

  # Perform UMAP
  umap_res <- uwot::umap(
    X = mat,
    n_components = n_components,
    n_neighbors = n_neighbors,
    verbose = FALSE,
    ...
  )

  # Create tidy result tibble
  tidy_result <- tibble::tibble(
    sample = rownames(mat),
    umap1 = umap_res[, 1],
    umap2 = umap_res[, 2]
  )

  # Handle more than 2 components
  if (n_components > 2) {
    for (i in 3:n_components) {
      tidy_result[[paste0("umap", i)]] <- umap_res[, i]
    }
  }

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = umap_res
    ),
    class = c("glystats_umap_res", "glystats_res")
  )
}