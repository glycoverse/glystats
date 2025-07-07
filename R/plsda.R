#' Partial Least Squares Discriminant Analysis (PLS-DA)
#'
#' Perform partial least squares discriminant analysis on the expression data.
#' The function uses `mixOmics::plsda()` to perform PLS-DA and returns tidy results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param group_col A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param ncomp An integer indicating the number of components to include. Default is 2.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw mixOmics plsda object.
#' @param ... Additional arguments passed to `mixOmics::plsda()`.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `mixOmics` for PLS-DA analysis
#'
#' @section Sample size requirements:
#' According to the Topliss ratio principle, the ratio of samples to variables (n/p)
#' should be at least 5 to avoid overfitting and ensure reliable results. This function
#' will throw an error if n/p < 5. For datasets with high dimensionality relative to
#' sample size, consider:
#' - Feature selection before analysis
#' - Collecting more samples
#'
#' @return A list containing four tibbles (when return_raw = FALSE):
#'  - `samples`: PLS-DA scores for each sample with group information
#'  - `variables`: PLS-DA loadings for each variable
#'  - `variance`: PLS-DA explained variance information
#'  - `vip`: Variable Importance in Projection (VIP) scores for each variable
#' When return_raw = TRUE, returns the raw mixOmics plsda object.
#' @seealso [mixOmics::plsda()]
#' @export
gly_plsda <- function(exp, group_col = "group", ncomp = 2, scale = TRUE, return_raw = FALSE, ...) {
  # Check package availability
  .check_pkg_available("mixOmics")

  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_int(ncomp, lower = 1)
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

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
      "*" = "Collecting more samples (need ≥ {ceiling(n_variables >= 5)} samples)",
      "*" = "Reducing variables through feature selection",
      "*" = "Using dimensionality reduction techniques first"
    ))
  }

  # Prepare data matrix (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)

  # Perform PLS-DA
  plsda_res <- mixOmics::plsda(X = mat, Y = groups, ncomp = ncomp, scale = scale, ...)

  # Return raw results if requested
  if (return_raw) {
    return(plsda_res)
  }

  # Extract and format results
  res <- .format_plsda_results(plsda_res, groups, sample_info)
  structure(res, class = c("glystats_plsda_res", "glystats_res"))
}

# Helper function to format PLS-DA results
.format_plsda_results <- function(plsda_res, groups, sample_info) {
  # Extract sample scores
  samples_tbl <- tibble::as_tibble(plsda_res$variates$X, .name_repair = "minimal")
  colnames(samples_tbl) <- paste0("comp", seq_len(ncol(samples_tbl)))
  samples_tbl$sample <- rownames(plsda_res$variates$X)
  samples_tbl$group <- as.character(groups)

  # Add sample information if available (excluding group column to avoid duplication)
  if (!is.null(sample_info) && "sample" %in% colnames(sample_info)) {
    sample_info_subset <- sample_info[, !colnames(sample_info) %in% "group", drop = FALSE]
    if (ncol(sample_info_subset) > 1) {  # Only join if there are columns other than sample
      samples_tbl <- dplyr::left_join(samples_tbl, sample_info_subset, by = "sample")
    }
  }

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
