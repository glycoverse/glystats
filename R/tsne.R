#' t-Distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' Perform t-SNE dimensionality reduction on the expression data.
#' The function uses `Rtsne::Rtsne()` to perform t-SNE analysis.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param dims Number of output dimensions. Default is 2.
#' @param perplexity Perplexity parameter for t-SNE. Default is 30.
#' @param theta Speed/accuracy trade-off parameter. Default is 0.5.
#' @param max_iter Maximum number of iterations. Default is 1000.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw Rtsne object.
#' @param ... Additional arguments passed to `Rtsne::Rtsne()`.
#'
#' @section Required packages:
#' This function requires the `Rtsne` package to be installed for t-SNE analysis.
#'
#' @return A tibble with t-SNE coordinates (tsne1, tsne2) when return_raw = FALSE,
#'   or raw Rtsne object when return_raw = TRUE.
#' 
#' @examples
#' \dontrun{
#' tsne_res <- gly_tsne(exp)
#' }
#'
#' @export
gly_tsne <- function(exp, dims = 2, perplexity = 30, theta = 0.5, max_iter = 1000, return_raw = FALSE, ...) {
  
  .check_pkg_available("Rtsne")
  
  checkmate::check_logical(return_raw, len = 1)
  
  # Extract expression matrix and sample info
  mat <- t(exp$expr_mat)  # Samples as rows, variables as columns
  sample_info <- exp$sample_info
  
  # Check if perplexity is appropriate for the number of samples
  # Perplexity should be smaller than the number of samples
  max_perplexity <- floor((nrow(mat) - 1) / 3)
  if (perplexity >= nrow(mat) || perplexity > max_perplexity) {
    cli::cli_warn("Perplexity should be smaller than the number of samples.")
    perplexity <- max(1, max_perplexity)
    cli::cli_inform("Setting perplexity to {perplexity}.")
  }
  
  # Apply log transformation (common for expression data)
  mat <- log(mat + 1)
  
  # Perform t-SNE
  tsne_res <- Rtsne::Rtsne(
    X = mat,
    dims = dims,
    perplexity = perplexity,
    theta = theta,
    max_iter = max_iter,
    check_duplicates = FALSE,  # Set to FALSE for better performance
    ...
  )
  
  # Return raw results if requested
  if (return_raw) {
    return(tsne_res)
  }

  # Create result tibble
  result <- tibble::tibble(
    sample = rownames(mat),
    tsne1 = tsne_res$Y[, 1],
    tsne2 = tsne_res$Y[, 2]
  )
  
  # Set S3 class
  structure(result, class = c("glystats_tsne_res", "glystats_res", class(result)))
}
