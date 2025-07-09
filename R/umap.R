#' Uniform Manifold Approximation and Projection (UMAP)
#'
#' Perform UMAP dimensionality reduction on the expression data.
#' The function uses `uwot::umap()` to perform UMAP analysis.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param n_neighbors Number of neighbors to consider for each point. Default is 15.
#' @param n_components Number of output dimensions. Default is 2.
#' @param min_dist Minimum distance between embedded points. Default is 0.1.
#' @param spread Controls how tightly the embedding is packed. Default is 1.0.
#' @param metric Distance metric to use. Default is "euclidean".
#' @param n_epochs Number of training epochs. Default is 200.
#' @param learning_rate Learning rate for the optimization. Default is 1.0.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'  will be added to the result tibble. If FALSE, only the UMAP coordinates are returned.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw umap result matrix.
#' @param ... Additional arguments passed to `uwot::umap()`.
#'
#' @section Required packages:
#' This function requires the `uwot` package to be installed for UMAP analysis.
#'
#' @return A tibble with UMAP coordinates (umap1, umap2) when return_raw = FALSE,
#'   or raw umap result matrix when return_raw = TRUE.
#' @seealso [uwot::umap()]
#' @export
gly_umap <- function(exp,
                     n_neighbors = 15,
                     n_components = 2,
                     min_dist = 0.1,
                     spread = 1.0,
                     metric = "euclidean",
                     n_epochs = 200,
                     learning_rate = 1.0,
                     add_info = TRUE,
                     return_raw = FALSE,
                     ...) {

  .check_pkg_available("uwot")

  checkmate::check_logical(add_info, len = 1)
  checkmate::check_logical(return_raw, len = 1)

  # Extract expression matrix and sample info
  mat <- t(exp$expr_mat)  # Samples as rows, variables as columns
  sample_info <- exp$sample_info

  # Check if n_neighbors is appropriate for the number of samples
  if (n_neighbors >= nrow(mat)) {
    cli::cli_warn("n_neighbors should be smaller than the number of samples.")
    n_neighbors <- max(1, nrow(mat) - 1)
    cli::cli_inform("Setting n_neighbors to {n_neighbors}.")
  }

  # Apply log transformation (common for expression data)
  mat <- log(mat + 1)

  # Capture output to avoid console messages
  umap_res <- suppressMessages(
    uwot::umap(
      X = mat,
      n_neighbors = n_neighbors,
      n_components = n_components,
      min_dist = min_dist,
      spread = spread,
      metric = metric,
      n_epochs = n_epochs,
      learning_rate = learning_rate,
      verbose = FALSE,
      ...
    )
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

  # Process results with add_info logic
  result <- .process_results_add_info(result, exp, add_info)

  # Set S3 class
  structure(result, class = c("glystats_umap_res", "glystats_res", class(result)))
} 