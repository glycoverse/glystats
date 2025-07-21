#' t-Distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' Perform t-SNE dimensionality reduction on the expression data.
#' The function uses `Rtsne::Rtsne()` to perform t-SNE analysis.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param dims Number of output dimensions. Default is 2.
#' @param perplexity Perplexity parameter for t-SNE. Default is 30.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'  will be added to the result tibble. If FALSE, only the t-SNE coordinates are returned.
#'  Only applicable to `gly_tsne()`.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw Rtsne object.
#' @param ... Additional arguments passed to `Rtsne::Rtsne()`.
#'
#' @section Required packages:
#' This function requires the `Rtsne` package to be installed for t-SNE analysis.
#'
#' @details
#' `gly_tsne()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_tsne_()` is the underlying API that works with matrices directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @return A tibble with t-SNE coordinates (tsne1, tsne2) when return_raw = FALSE,
#'   or raw Rtsne object when return_raw = TRUE.
#' @seealso [Rtsne::Rtsne()]
#' @export
gly_tsne <- function(exp, dims = 2, perplexity = 30, add_info = TRUE, return_raw = FALSE, ...) {

  .check_pkg_available("Rtsne")

  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

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

  # Process results with add_info logic
  result <- .process_results_add_info(result, exp, add_info)

  # Set S3 class
  structure(result, class = c("glystats_tsne_res", "glystats_res", class(result)))
}

#' @rdname gly_tsne
#' @export
gly_tsne_ <- function(expr_mat, dims = 2, perplexity = 30, return_raw = FALSE, ...) {
  .check_pkg_available("Rtsne")

  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_logical(return_raw, len = 1)

  # Prepare data (samples as rows, variables as columns)
  mat <- t(expr_mat)

  # Check if perplexity is appropriate for the number of samples
  max_perplexity <- floor((nrow(mat) - 1) / 3)
  if (perplexity >= nrow(mat) || perplexity > max_perplexity) {
    cli::cli_warn("Perplexity should be smaller than the number of samples.")
    perplexity <- max(1, max_perplexity)
    cli::cli_inform("Setting perplexity to {perplexity}.")
  }

  # Apply log transformation
  mat <- log(mat + 1)

  # Perform t-SNE
  tsne_res <- Rtsne::Rtsne(
    X = mat,
    dims = dims,
    perplexity = perplexity,
    ...
  )

  # Return raw results if requested
  if (return_raw) {
    return(tsne_res)
  }

  # Create result tibble
  tsne_coords <- tsne_res$Y
  colnames(tsne_coords) <- paste0("tsne", seq_len(ncol(tsne_coords)))

  result_tbl <- tibble::tibble(
    sample = rownames(mat)
  ) %>%
    dplyr::bind_cols(tibble::as_tibble(tsne_coords))

  # Add S3 class
  structure(result_tbl, class = c("glystats_tsne_res", "glystats_res", class(result_tbl)))
}
