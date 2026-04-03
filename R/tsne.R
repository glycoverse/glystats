#' t-Distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' Perform t-SNE dimensionality reduction on the expression data.
#' The function uses `Rtsne::Rtsne()` to perform t-SNE analysis.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat (Only for [gly_tsne_()]) A numeric matrix with variables as rows and samples as columns.
#' @param dims Number of output dimensions. Default is 2.
#' @param perplexity Perplexity parameter for t-SNE. Default is 30.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'  will be added to the result tibble. If FALSE, only the t-SNE coordinates are returned.
#'  Only applicable to `gly_tsne()`.
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
#' @return A list with three elements:
#' - `tidy_result`: A tibble with t-SNE coordinates containing the following columns:
#'   - `sample`: Sample name
#'   - `tsne1`: First t-SNE dimension
#'   - `tsne2`: Second t-SNE dimension
#' - `raw_result`: The raw Rtsne object
#' - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_tsne_res` and `glystats_res`.
#' @seealso [Rtsne::Rtsne()]
#' @export
gly_tsne <- function(exp, dims = 2, perplexity = 30, add_info = TRUE, ...) {
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_logical(add_info, len = 1)

  expr_mat <- glyexp::get_expr_mat(exp)
  result <- gly_tsne_(expr_mat, dims, perplexity, ...)

  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' @rdname gly_tsne
#' @export
gly_tsne_ <- function(expr_mat, dims = 2, perplexity = 30, ...) {
  rlang::check_installed("Rtsne")

  checkmate::assert_matrix(expr_mat, mode = "numeric")

  mat <- t(expr_mat) # Samples as rows, variables as columns

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
  dots <- rlang::list2(...)
  if ("X" %in% names(dots)) {
    cli::cli_abort(
      "{.field X} should not be supplied through `...`; data comes from the function inputs."
    )
  }
  call_args <- c(
    list(
      X = mat,
      dims = dims,
      perplexity = perplexity
    ),
    dots
  )
  if (!"check_duplicates" %in% names(call_args)) {
    call_args$check_duplicates <- FALSE
  }
  tsne_res <- do.call(Rtsne::Rtsne, call_args)

  # Create tidy result tibble
  tidy_result <- tibble::tibble(
    sample = rownames(mat),
    tsne1 = tsne_res$Y[, 1],
    tsne2 = tsne_res$Y[, 2]
  )

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = tsne_res
    ),
    class = c("glystats_tsne_res", "glystats_res")
  )
}
