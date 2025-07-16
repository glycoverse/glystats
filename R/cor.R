#' Correlation Analysis for Glycomics and Glycoproteomics Data
#'
#' Perform pairwise correlation analysis on variables or samples in the expression data.
#' The function calculates correlation coefficients and p-values for all pairs,
#' with optional multiple testing correction.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param on A character string specifying what to correlate. Either "variable" (default) to correlate
#'   variables/features, or "sample" to correlate samples/observations.
#' @param method A character string indicating which correlation coefficient is to be computed.
#'   One of "pearson" (default) or "spearman". Note: "kendall" is not supported by Hmisc::rcorr.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'   See `p.adjust.methods` for available methods. Default is "BH".
#'   If NULL, no adjustment is performed.

#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw rcorr object.
#' @param ... Additional arguments passed to `Hmisc::rcorr()`.
#'
#' @section Required packages:
#' This function requires the `Hmisc` package for efficient correlation calculation.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' correlation analysis. When `on = "variable"` (default), correlations are calculated between
#' variables across samples. When `on = "sample"`, correlations are calculated between
#' samples across variables.
#'
#' **Correlation Calculation:**
#' Correlation coefficients and p-values are calculated using `Hmisc::rcorr()` which is
#' more efficient than pairwise `stats::cor.test()` calls.
#'
#' **Multiple Testing Correction:**
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @return A tibble with correlation results (when return_raw = FALSE):
#'  - First column: First element of the pair (variable1/sample1)
#'  - Second column: Second element of the pair (variable2/sample2)
#'  - `cor`: Correlation coefficient
#'  - `p_value`: Raw p-value
#'  - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#' When return_raw = TRUE, returns the raw rcorr object.
#'
#' @seealso [Hmisc::rcorr()], [stats::cor()]
#' @export
gly_cor <- function(
  exp,
  on = "variable",
  method = "pearson",
  p_adj_method = "BH",
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_choice(on, c("variable", "sample"))
  checkmate::assert_choice(method, c("pearson", "spearman"))
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_logical(return_raw, len = 1)

  # Check if Hmisc is available
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop("Package 'Hmisc' is required for gly_cor. Please install it with: install.packages('Hmisc')")
  }

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)

  # Prepare data for correlation based on 'on' parameter
  # expr_mat is variables x samples (variables as rows, samples as columns)
  if (on == "sample") {
    # Correlate samples: need samples as columns for rcorr
    mat <- expr_mat
    cor_type <- "sample"
  } else {
    # Correlate variables: need variables as columns for rcorr
    mat <- t(expr_mat)
    cor_type <- "variable"
  }

  # Apply log transformation
  mat <- log2(mat + 1)

  # Calculate correlation using Hmisc::rcorr
  rcorr_result <- Hmisc::rcorr(mat, type = method, ...)

  # Return raw results if requested
  if (return_raw) {
    return(rcorr_result)
  }

  # Extract correlation and p-value matrices from rcorr result
  cor_matrix <- rcorr_result$r
  p_matrix <- rcorr_result$P

  # Convert matrices to long format tibble, excluding diagonal and duplicates
  # Create a tibble with all combinations, then filter for upper triangle
  result_tbl <- tibble::tibble(
    item1 = rep(rownames(cor_matrix), each = ncol(cor_matrix)),
    item2 = rep(colnames(cor_matrix), times = nrow(cor_matrix)),
    cor = as.vector(cor_matrix),
    p_value = as.vector(p_matrix)
  ) %>%
    # Add row indices to filter upper triangle
    dplyr::mutate(
      row_idx = rep(seq_len(nrow(cor_matrix)), each = ncol(cor_matrix)),
      col_idx = rep(seq_len(ncol(cor_matrix)), times = nrow(cor_matrix))
    ) %>%
    # Keep only upper triangle (excluding diagonal)
    dplyr::filter(.data$row_idx < .data$col_idx) %>%
    # Remove helper columns
    dplyr::select(-c("row_idx", "col_idx"))

  # Set appropriate column names based on correlation type
  if (cor_type == "sample") {
    colnames(result_tbl)[1:2] <- c("sample1", "sample2")
  } else {
    colnames(result_tbl)[1:2] <- c("variable1", "variable2")
  }

  # Adjust p-values if requested
  if (!is.null(p_adj_method)) {
    result_tbl$p_adj <- stats::p.adjust(result_tbl$p_value, method = p_adj_method)
  }

  structure(result_tbl, class = c("glystats_cor_res", "glystats_res", class(result_tbl)))
}