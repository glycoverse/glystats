#' t-Distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' Perform t-SNE dimensionality reduction on the expression data.
#' The function uses `Rtsne::Rtsne()` to perform t-SNE analysis.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing an expression matrix and
#'   sample information.
#' @param dims Number of output dimensions. Default is 2.
#' @param perplexity Perplexity parameter for t-SNE. Default is 30.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'  will be added to the result tibble. If FALSE, only the t-SNE coordinates are returned.
#' @param ... Additional arguments passed to `Rtsne::Rtsne()`.
#'
#' @section Required packages:
#' This function requires the `Rtsne` package to be installed for t-SNE analysis.
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
  .assert_data_container(exp)
  checkmate::assert_logical(add_info, len = 1)

  expr_mat <- .get_expr_mat(exp)
  result <- .analyze_tsne(expr_mat, dims, perplexity, ...)

  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_tsne()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_tsne <- function(expr_mat, dims = 2, perplexity = 30, ...) {
  rlang::check_installed("Rtsne")

  checkmate::assert_matrix(expr_mat, mode = "numeric")

  sample_order <- colnames(expr_mat)
  expr_mat <- .remove_all_na_variables(expr_mat)
  if (nrow(expr_mat) == 0) {
    return(structure(
      list(
        tidy_result = tibble::tibble(
          sample = sample_order,
          tsne1 = NA_real_,
          tsne2 = NA_real_
        ),
        raw_result = NULL
      ),
      class = c("glystats_tsne_res", "glystats_res")
    ))
  }

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
