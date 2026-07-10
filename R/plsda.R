#' Partial Least Squares Discriminant Analysis (PLS-DA)
#'
#' Perform partial least squares discriminant analysis on the expression data.
#' The function uses `ropls::opls()` to perform PLS-DA and returns tidy results.
#'
#' @param exp A `glyexp::experiment()` or `SummarizedExperiment` object containing
#'   an expression matrix and sample information.
#' @param group_col A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param ncomp An integer indicating the number of components to include. Default is 2.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the PLS-DA results are returned.
#' @param ... Additional arguments passed to `ropls::opls()`.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `ropls` for PLS-DA analysis
#'
#' @return A list containing:
#'  - `tidy_result`: A list of tibbles with PLS-DA results:
#'    - `samples`: PLS-DA scores for each sample containing the following columns:
#'      - `sample`: Sample name
#'      - `group`: Group assignment
#'      - `p1`, `p2`, etc.: PLS-DA component scores
#'    - `variables`: PLS-DA loadings for each variable containing the following columns:
#'      - `variable`: Variable name
#'      - `p1`, `p2`, etc.: PLS-DA component loadings
#'      - `pcorr1`, `pcorr2`, etc.: Correlation between each variable and component
#'    - `variance`: PLS-DA explained variance containing the following columns:
#'      - `component`: Component name (p1, p2, etc.)
#'      - `prop_var_explained`: Proportion of variance explained by each component
#'      - `cumulative_prop_var`: Cumulative proportion of variance explained
#'    - `vip`: Variable Importance in Projection scores containing the following columns:
#'      - `variable`: Variable name
#'      - `vip`: VIP score
#'    - `perm_test`: Permutation test results containing the following columns:
#'      - `model`: Model type ("Original" for the original model, "Permutation" for permuted models)
#'      - `perm_id`: Permutation ID (0 for original model, 1+ for permutations)
#'      - Additional columns from the permutation test matrix (e.g., R2X, R2Y, Q2, etc.)
#'  - `raw_result`: The raw ropls opls object from `ropls::opls()`
#'  - `meta_data`: A list containing metadata from the input experiment
#' @seealso [ropls::opls()]
#' @export
gly_plsda <- function(
  exp,
  group_col = "group",
  ncomp = 2,
  scale = TRUE,
  add_info = TRUE,
  ...
) {
  # Check package availability
  rlang::check_installed("ropls")

  # Validate inputs
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_int(ncomp, lower = 1)
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "PLS-DA"
  )
  groups <- group_info$groups

  # Run the internal computation
  result <- .analyze_plsda(expr_mat, groups, ncomp, scale, ...)
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_plsda()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_plsda <- function(expr_mat, groups, ncomp = 2, scale = TRUE, ...) {
  rlang::check_installed("ropls")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))

  # Prepare data matrix (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)

  # Perform PLS-DA using ropls::opls with orthoI = 0
  # Set appropriate cross-validation folds based on sample size
  n_samples <- nrow(mat)
  crossval_i <- min(7, n_samples - 1) # Default is 7, but must be less than sample size

  # Suppress plotting to prevent Rplots.pdf generation
  # Open a null device to capture any plotting output
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # For small datasets, reduce permutation tests to allow model building
  perm_i <- if (n_samples < 10) 0 else 20

  dots <- rlang::list2(...)
  disallowed_args <- intersect(names(dots), c("x", "y"))
  if (length(disallowed_args) > 0) {
    cli::cli_abort(
      "Arguments {cli::format_inline(disallowed_args)} should not be supplied through `...`; data comes from the function inputs."
    )
  }

  call_args <- c(
    list(
      x = mat,
      y = groups,
      predI = ncomp
    ),
    dots
  )
  if (!"orthoI" %in% names(call_args)) {
    call_args$orthoI <- 0
  }
  if (!"scaleC" %in% names(call_args)) {
    call_args$scaleC <- if (scale) "standard" else "none"
  }
  if (!"crossvalI" %in% names(call_args)) {
    call_args$crossvalI <- crossval_i
  }
  if (!"permI" %in% names(call_args)) {
    call_args$permI <- perm_i
  }
  if (!"fig.pdfC" %in% names(call_args)) {
    call_args$fig.pdfC <- "none"
  }
  if (!"info.txtC" %in% names(call_args)) {
    call_args$info.txtC <- "none"
  }

  # Perform PLS-DA by setting orthoI = 0 (no orthogonal components)
  plsda_res <- do.call(ropls::opls, call_args)

  # Extract and format results
  tidy_result <- .format_oplsda_results(plsda_res, groups, NULL, FALSE)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = plsda_res
    ),
    class = c("glystats_plsda_res", "glystats_res")
  )
}
