#' Principal Component Analysis (PCA)
#'
#' Perform principal component analysis on the expression data.
#' The function uses `prcomp()` to perform PCA and `broom::tidy()` to tidy the results.
#' If `scale = TRUE`, constant variables (zero variance) will be removed before PCA.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing an expression matrix and
#'   sample information.
#' @param center A logical indicating whether to center the data. Default is TRUE.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the PCA results are returned.
#' @param ... Additional arguments passed to `prcomp()`.
#'
#' @section Required packages:
#' This function only uses base R packages and does not require additional dependencies.
#'
#' @details
#' The function performs log transformation on the expression data (log(x + 1)) before
#' PCA analysis.
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
#'  - `meta_data`: A list containing metadata from the input experiment
#' @seealso [stats::prcomp()]
#' @export
gly_pca <- function(exp, center = TRUE, scale = TRUE, add_info = TRUE, ...) {
  .assert_data_container(exp)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- .get_expr_mat(exp)

  # Run the internal computation
  result <- .analyze_pca(expr_mat, center, scale, ...)
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_pca()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_pca <- function(expr_mat, center = TRUE, scale = TRUE, ...) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_logical(center, len = 1)
  checkmate::assert_logical(scale, len = 1)

  variable_order <- rownames(expr_mat)
  missing_variables <- .all_na_variables(expr_mat)
  expr_mat <- .remove_all_na_variables(expr_mat)

  if (nrow(expr_mat) == 0) {
    return(structure(
      list(
        tidy_result = list(
          samples = tibble::tibble(
            sample = colnames(expr_mat),
            PC = NA_integer_,
            value = NA_real_
          ),
          variables = tibble::tibble(
            variable = variable_order,
            PC = NA_integer_,
            value = NA_real_
          ),
          eigenvalues = tibble::tibble(
            PC = integer(),
            std.dev = double(),
            percent = double(),
            cumulative = double()
          )
        ),
        raw_result = NULL
      ),
      class = c("glystats_pca_res", "glystats_res")
    ))
  }

  # Prepare data for PCA (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)

  # If scaling, check for constant columns (zero variance)
  # These will cause prcomp to fail with "cannot rescale a constant/zero column to unit variance"
  if (scale) {
    col_vars <- apply(mat, 2, stats::var, na.rm = TRUE)
    constant_cols <- which(col_vars == 0 | is.na(col_vars))

    if (length(constant_cols) > 0) {
      constant_names <- colnames(mat)[constant_cols]
      cli::cli_warn(c(
        "Removed {length(constant_cols)} constant column{?s} before PCA (zero variance after log transformation).",
        "i" = "Affected variable{?s}: {.val {constant_names}}"
      ))
      mat <- mat[, -constant_cols, drop = FALSE]

      # Check if we have any columns left
      if (ncol(mat) == 0) {
        cli::cli_abort(
          "No columns with non-zero variance remain after removing constant columns."
        )
      }
    }
  }

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
  variables_tbl <- .restore_all_na_variable_components(
    variables_tbl,
    missing_variables,
    variable_order,
    "PC"
  )

  # Create tidy result list
  tidy_result <- list(
    "samples" = samples_tbl,
    "variables" = variables_tbl,
    "eigenvalues" = eigenvalues_tbl
  )

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = prcomp_res
    ),
    class = c("glystats_pca_res", "glystats_res")
  )
}
