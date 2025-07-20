#' Linear Models for Differential Expression Analysis using limma
#'
#' Perform differential expression analysis using linear models with empirical Bayes
#' moderation from the limma package.
#'
#' @param exp A `glyexp_experiment` object containing expression data and sample information.
#' @param group_col A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param p_adj_method A character string specifying the method for multiple testing correction.
#'   Must be one of the methods supported by `stats::p.adjust()`. Default is "BH" (Benjamini-Hochberg).
#'   Set to NULL to skip p-value adjustment.
#' @param ref_group A character string specifying the reference group.
#'  If NULL (default), the first level of the group factor is used as the reference.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw limma fit objects as a list.
#' @param ... Additional arguments passed to `limma::lmFit()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. The analysis uses linear models with empirical Bayes moderation
#' to improve statistical power, especially for small sample sizes. Exactly 2 groups are
#' required in the grouping variable.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `limma` for linear model fitting and empirical Bayes moderation
#'
#' @returns A tibble with limma results including log2 fold change (log2fc),
#'   or a list of limma fit objects if `return_raw` is TRUE.
#' @seealso [limma::lmFit()], [limma::eBayes()]
#' @export
gly_limma <- function(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  ref_group = NULL,
  add_info = TRUE,
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

  # Check package availability
  .check_pkg_available("limma")

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "limma"
  )
  groups <- group_info$groups
  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Perform limma analysis
  result <- .gly_limma_2groups(expr_mat, groups, p_adj_method, ref_group, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Process results with add_info logic
  result <- .process_results_add_info(result, exp, add_info)

  # Add S3 class
  structure(result, class = c("glystats_dea_res_limma", "glystats_dea_res", "glystats_res", class(result)))
}

# Internal helper functions for limma analysis

# Limma-specific 2-group analysis function
.gly_limma_2groups <- function(expr_mat, groups, p_adj_method, ref_group, return_raw = FALSE, ...) {
  # Reorder groups if ref_group is specified
  if (!is.null(ref_group)) {
    groups <- .reorder_groups_for_ref(groups, ref_group)
  }
  cli::cli_alert_info("Reference group: {.val {levels(groups)[1]}}")

  # Prepare data for limma
  log_expr_mat <- log2(expr_mat + 1)

  # Create design matrix
  design <- stats::model.matrix(~ groups)
  colnames(design) <- c("Intercept", paste0(levels(groups)[2], "_vs_", levels(groups)[1]))

  # Fit linear model
  fit <- limma::lmFit(log_expr_mat, design, ...)

  # Apply empirical Bayes moderation
  fit <- limma::eBayes(fit)

  if (return_raw) {
    return(fit)
  }

  # Extract results and convert to tibble
  .gly_limma_tibblify(fit, p_adj_method, expr_mat, groups)
}

# Convert limma fit object to tibble
.gly_limma_tibblify <- function(fit, p_adj_method, expr_mat, groups) {
  # Use topTable to extract results - this is the standard limma approach
  # The coefficient of interest is the second column (group comparison)
  coef_idx <- 2

  # Extract all results using topTable
  top_table <- limma::topTable(
    fit,
    coef = coef_idx,
    number = Inf,  # Get all genes
    adjust.method = if (is.null(p_adj_method)) "none" else p_adj_method,
    sort.by = "none"  # Keep original order
  )

  # Convert to tibble and rename columns to match glystats conventions
  result_tbl <- tibble::as_tibble(top_table, rownames = "variable") %>%
    dplyr::rename(
      log2fc = "logFC",
      p = "P.Value",
      p_adj = "adj.P.Val",
      t = "t",
      b = "B"
    )

  # Remove p_adj column if p_adj_method was NULL
  if (is.null(p_adj_method)) {
    result_tbl <- dplyr::select(result_tbl, -"p_adj")
  }

  result_tbl
}

# Helper function to reorder groups so that ref_group becomes the first level
.reorder_groups_for_ref <- function(groups, ref_group) {
  current_levels <- levels(groups)

  # Move ref_group to the first position
  new_levels <- c(ref_group, setdiff(current_levels, ref_group))

  # Reorder the factor levels
  factor(groups, levels = new_levels)
}
