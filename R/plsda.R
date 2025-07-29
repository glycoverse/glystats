#' Partial Least Squares Discriminant Analysis (PLS-DA)
#'
#' Perform partial least squares discriminant analysis on the expression data.
#' The function uses `mixOmics::plsda()` to perform PLS-DA and returns tidy results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor or character vector specifying group membership for each sample.
#'   Character vectors will be automatically converted to factors.
#' @param group_col A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param ncomp An integer indicating the number of components to include. Default is 2.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the PLS-DA results are returned.
#'  Only applicable to `gly_plsda()`.
#' @param ... Additional arguments passed to `mixOmics::plsda()`.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `mixOmics` for PLS-DA analysis
#'
#' @details
#' `gly_plsda()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_plsda_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @section Sample size requirements:
#' According to the Topliss ratio principle, the ratio of samples to variables (n/p)
#' should be at least 5 to avoid overfitting and ensure reliable results. This function
#' will throw an error if n/p < 5. For datasets with high dimensionality relative to
#' sample size, consider:
#' - Feature selection before analysis
#' - Collecting more samples
#'
#' @return A list containing:
#'  - `tidy_result`: A list of tibbles with PLS-DA results:
#'    - `samples`: PLS-DA scores for each sample containing the following columns:
#'      - `sample`: Sample name
#'      - `group`: Group assignment
#'      - `comp1`, `comp2`, etc.: PLS-DA component scores
#'    - `variables`: PLS-DA loadings for each variable containing the following columns:
#'      - `variable`: Variable name
#'      - `comp1`, `comp2`, etc.: PLS-DA component loadings
#'    - `variance`: PLS-DA explained variance containing the following columns:
#'      - `comp`: Component name (comp1, comp2, etc.)
#'      - `explained_variance`: Percentage of variance explained by each component
#'    - `vip`: Variable Importance in Projection scores containing the following columns:
#'      - `variable`: Variable name
#'      - `vip`: VIP score
#'  - `raw_result`: The raw mixOmics plsda object from `mixOmics::plsda()`
#' @seealso [mixOmics::plsda()]
#' @export
gly_plsda <- function(exp, group_col = "group", ncomp = 2, scale = TRUE, add_info = TRUE, ...) {
  # Check package availability
  .check_pkg_available("mixOmics")

  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_int(ncomp, lower = 1)
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "PLS-DA"
  )
  groups <- group_info$groups

  # Validate sample-to-variable ratio (Topliss ratio)
  n_samples <- length(groups)
  n_variables <- nrow(expr_mat)
  topliss_ratio <- n_samples / n_variables

  if (topliss_ratio < 5) {
    cli::cli_abort(c(
      "Insufficient sample-to-variable ratio for reliable PLS-DA analysis.",
      "x" = "Current ratio: {n_samples}/{n_variables} = {round(topliss_ratio, 2)}",
      "!" = "According to the Topliss ratio principle, n/p should be >= 5 to avoid overfitting.",
      "i" = "Consider:",
      "*" = "Collecting more samples (need >= {ceiling(n_variables * 5)} samples)",
      "*" = "Reducing variables through feature selection"
    ))
  }

  # Call the underlying API
  result <- gly_plsda_(expr_mat, groups, ncomp, scale, ...)
  result$tidy_result <- .process_results_add_info(result$tidy_result, exp, add_info)
  result
}

#' @rdname gly_plsda
#' @export
gly_plsda_ <- function(expr_mat, groups, ncomp = 2, scale = TRUE, ...) {
  .check_pkg_available("mixOmics")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  groups <- .convert_groups_to_factor(groups)
  checkmate::assert_factor(groups, len = ncol(expr_mat))

  # Prepare data (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)

  # Perform PLS-DA
  plsda_res <- mixOmics::plsda(X = mat, Y = groups, ncomp = ncomp, scale = scale, ...)

  # Extract and format results
  tidy_result <- .format_plsda_results(plsda_res, groups, NULL, FALSE)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = plsda_res
    ),
    class = c("glystats_plsda_res", "glystats_res")
  )
}

# Helper function to format PLS-DA results
.format_plsda_results <- function(plsda_res, groups, sample_info, add_info = TRUE) {
  # Extract sample scores
  samples_tbl <- tibble::as_tibble(plsda_res$variates$X, .name_repair = "minimal")
  colnames(samples_tbl) <- paste0("comp", seq_len(ncol(samples_tbl)))
  samples_tbl$sample <- rownames(plsda_res$variates$X)

  # Note: We don't add group or sample_info here because it will be handled by .process_results_add_info
  # This ensures consistent behavior across all functions

  # Extract variable loadings
  variables_tbl <- tibble::as_tibble(plsda_res$loadings$X, .name_repair = "minimal")
  colnames(variables_tbl) <- paste0("comp", seq_len(ncol(variables_tbl)))
  variables_tbl$variable <- rownames(plsda_res$loadings$X)

  # Extract explained variance information
  # For PLS-DA, we use prop_expl_var$X which contains the proportion of variance explained
  n_comp <- length(plsda_res$prop_expl_var$X)

  variance_tbl <- tibble::tibble(
    component = paste0("comp", seq_len(n_comp)),
    prop_var_explained = plsda_res$prop_expl_var$X,
    cumulative_prop_var = cumsum(plsda_res$prop_expl_var$X)
  )

  # Calculate VIP (Variable Importance in Projection) scores
  vip_matrix <- mixOmics::vip(plsda_res)
  # Calculate overall VIP score as the square root of the sum of squared VIP scores across components
  overall_vip <- sqrt(rowSums(vip_matrix^2))

  vip_tbl <- tibble::tibble(
    variable = rownames(vip_matrix),
    VIP = overall_vip
  )

  list(
    samples = samples_tbl,
    variables = variables_tbl,
    variance = variance_tbl,
    vip = vip_tbl
  )
}
