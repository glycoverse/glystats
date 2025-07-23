#' Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#'
#' Perform orthogonal partial least squares discriminant analysis on the expression data.
#' The function uses `ropls::opls()` to perform OPLS-DA and returns tidy results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor vector specifying group membership for each sample.
#' @param group_col A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param pred_i An integer indicating the number of predictive components to include. Default is 1.
#' @param ortho_i An integer indicating the number of orthogonal components to include. Default is NA (automatic).
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the OPLS-DA results are returned.
#'  Only applicable to `gly_oplsda()`.
#' @param ... Additional arguments passed to `ropls::opls()`.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `ropls` for OPLS-DA analysis
#'
#' @details
#' `gly_oplsda()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_oplsda_()` is the underlying API that works with matrices and factor vectors directly,
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
#'  - `tidy_result`: A list of tibbles with OPLS-DA results:
#'    - `samples`: OPLS-DA scores for each sample with group information
#'    - `variables`: OPLS-DA loadings for each variable
#'    - `variance`: OPLS-DA explained variance information
#'    - `vip`: Variable Importance in Projection (VIP) scores for each variable
#'  - `raw_result`: The raw ropls opls object from `ropls::opls()`
#' @seealso [ropls::opls()]
#' @export
gly_oplsda <- function(exp, group_col = "group", pred_i = 1, ortho_i = NA, scale = TRUE, add_info = TRUE, ...) {
  # Check package availability
  .check_pkg_available("ropls")

  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_int(pred_i, lower = 1)
  if (!is.na(ortho_i)) {
    checkmate::assert_int(ortho_i, lower = 0)
  }
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups (OPLS-DA only supports binary classification)
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "OPLS-DA"
  )
  groups <- group_info$groups

  # Validate sample-to-variable ratio (Topliss ratio)
  n_samples <- length(groups)
  n_variables <- nrow(expr_mat)
  topliss_ratio <- n_samples / n_variables

  if (topliss_ratio < 5) {
    cli::cli_abort(c(
      "Insufficient sample-to-variable ratio for reliable OPLS-DA analysis.",
      "x" = "Current ratio: {n_samples}/{n_variables} = {round(topliss_ratio, 2)}",
      "!" = "According to the Topliss ratio principle, n/p should be >= 5 to avoid overfitting.",
      "i" = "Consider:",
      "*" = "Collecting more samples (need >= {ceiling(n_variables * 5)} samples)",
      "*" = "Reducing variables through feature selection"
    ))
  }

  # Call the underlying API
  result <- gly_oplsda_(expr_mat, groups, pred_i, ortho_i, scale, ...)
  result$tidy_result <- .process_results_add_info(result$tidy_result, exp, add_info)
  result
}

#' @rdname gly_oplsda
#' @export
gly_oplsda_ <- function(expr_mat, groups, pred_i = 1, ortho_i = NA, scale = TRUE, ...) {
  .check_pkg_available("ropls")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))

  # Prepare data matrix (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)

  # Perform OPLS-DA
  # Set appropriate cross-validation folds based on sample size
  n_samples <- nrow(mat)
  crossval_i <- min(7, n_samples - 1)  # Default is 7, but must be less than sample size

  # Suppress plotting to prevent Rplots.pdf generation
  # Open a null device to capture any plotting output
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # For small datasets, reduce permutation tests to allow model building
  perm_i <- if (n_samples < 10) 0 else 20

  # Handle ortho_i parameter - if NA, let ropls decide automatically
  if (is.na(ortho_i)) {
    oplsda_res <- ropls::opls(x = mat, y = groups, predI = pred_i,
                              scaleC = if (scale) "standard" else "none",
                              crossvalI = crossval_i, permI = perm_i,
                              fig.pdfC = "none", info.txtC = "none", ...)
  } else {
    oplsda_res <- ropls::opls(x = mat, y = groups, predI = pred_i, orthoI = ortho_i,
                              scaleC = if (scale) "standard" else "none",
                              crossvalI = crossval_i, permI = perm_i,
                              fig.pdfC = "none", info.txtC = "none", ...)
  }

  # Extract and format results
  tidy_result <- .format_oplsda_results(oplsda_res, groups, NULL, FALSE)

  # Add S3 class to tidy_result
  tidy_result <- structure(tidy_result, class = c("glystats_oplsda_res", "glystats_res"))

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = oplsda_res
    ),
    class = c("glystats_oplsda_res", "glystats_res")
  )
}

# Helper function to format OPLS-DA results
.format_oplsda_results <- function(oplsda_res, groups, sample_info, add_info = TRUE) {
  # Check if model was successfully built
  if (length(oplsda_res@scoreMN) == 0) {
    # Model building failed - provide informative error
    cli::cli_abort(c(
      "OPLS-DA model building failed.",
      "i" = "This usually happens when:",
      "*" = "The data doesn't have sufficient discriminatory power between groups",
      "*" = "The first predictive component is not statistically significant",
      "*" = "Sample size is too small relative to the number of variables",
      "!" = "Consider using PLS-DA instead, or check your data quality."
    ))
  }

  # Extract sample scores (predictive + orthogonal components)
  # ropls stores scores in @scoreMN for predictive and @orthoScoreMN for orthogonal
  pred_scores <- oplsda_res@scoreMN
  ortho_scores <- oplsda_res@orthoScoreMN

  # Combine predictive and orthogonal scores
  if (!is.null(ortho_scores) && ncol(ortho_scores) > 0 && nrow(ortho_scores) == nrow(pred_scores)) {
    all_scores <- cbind(pred_scores, ortho_scores)
    # Create column names: p1, p2, ... for predictive, o1, o2, ... for orthogonal
    pred_names <- paste0("p", seq_len(ncol(pred_scores)))
    ortho_names <- paste0("o", seq_len(ncol(ortho_scores)))
    colnames(all_scores) <- c(pred_names, ortho_names)
  } else {
    all_scores <- pred_scores
    colnames(all_scores) <- paste0("p", seq_len(ncol(pred_scores)))
  }

  samples_tbl <- tibble::as_tibble(all_scores, .name_repair = "minimal")
  samples_tbl$sample <- rownames(all_scores)

  # Note: We don't add group here because it will be handled by .process_results_add_info
  # This ensures consistent behavior across all functions

  # Extract variable loadings (predictive components)
  pred_loadings <- oplsda_res@loadingMN
  variables_tbl <- tibble::as_tibble(pred_loadings, .name_repair = "minimal")
  colnames(variables_tbl) <- paste0("p", seq_len(ncol(pred_loadings)))
  variables_tbl$variable <- rownames(pred_loadings)

  # Extract explained variance information from modelDF
  # modelDF contains individual R2X values for each component
  model_df <- oplsda_res@modelDF

  # Remove the 'sum' row if present
  if ("sum" %in% rownames(model_df)) {
    model_df <- model_df[rownames(model_df) != "sum", , drop = FALSE]
  }

  # Get component names and R2X values
  component_names <- rownames(model_df)
  r2x_values <- model_df[["R2X"]]

  # Remove NA values
  valid_idx <- !is.na(r2x_values)
  component_names <- component_names[valid_idx]
  r2x_values <- r2x_values[valid_idx]

  variance_tbl <- tibble::tibble(
    component = component_names,
    prop_var_explained = r2x_values,
    cumulative_prop_var = cumsum(r2x_values)
  )

  # Calculate VIP (Variable Importance in Projection) scores
  # ropls provides VIP scores in @vipVn
  vip_scores <- oplsda_res@vipVn

  vip_tbl <- tibble::tibble(
    variable = names(vip_scores),
    VIP = as.numeric(vip_scores)
  )

  list(
    samples = samples_tbl,
    variables = variables_tbl,
    variance = variance_tbl,
    vip = vip_tbl
  )
}
