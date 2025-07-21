#' Calculate fold change
#'
#' Calculate fold change for a given expression matrix and group information.
#' It could only be used for 2-group analysis.
#' When you run this function, you will see message about "Group 1" and "Group 2".
#' "Group 1" is the reference group, and "Group 2" is the test group.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor vector specifying group membership for each sample.
#'   Must have exactly 2 levels.
#' @param group_col The column name of the group information in the sample information.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the fold change results are returned.
#'  Only applicable to `gly_fold_change()`.
#'
#' @details
#' `gly_fold_change()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_fold_change_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @return A tibble with two columns: `variable` and `log2fc`.
#' @export
gly_fold_change <- function(exp, group_col = "group", add_info = TRUE) {
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_logical(add_info, len = 1)

  # Extract and validate groups using helper function
  group_info <- .extract_and_validate_groups(
    sample_info = exp$sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "fold_change",
    show_info = TRUE
  )
  groups <- group_info$groups

  expr_mat <- exp$expr_mat

  # Call the underlying API
  result <- gly_fold_change_(expr_mat, groups)

  # Process results with add_info logic
  result <- .process_results_add_info(result, exp, add_info)

  structure(result, class = c("glystats_fc_res", class(result)))
}

#' @rdname gly_fold_change
#' @export
gly_fold_change_ <- function(expr_mat, groups) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))

  # Validate group count
  if (length(levels(groups)) != 2) {
    cli::cli_abort("groups must have exactly 2 levels for fold change calculation")
  }

  # Calculate mean expression for each group
  group1_mean <- rowMeans(expr_mat[, groups == levels(groups)[1], drop = FALSE])
  group2_mean <- rowMeans(expr_mat[, groups == levels(groups)[2], drop = FALSE])

  # Calculate log2 fold change
  log2fc <- log2(group2_mean / group1_mean)

  # Create result tibble
  result <- tibble::tibble(variable = rownames(expr_mat), log2fc = log2fc)

  # Add S3 class
  structure(result, class = c("glystats_fc_res", class(result)))
}
