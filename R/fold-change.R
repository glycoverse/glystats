#' Calculate fold change
#'
#' Calculate fold change for a given expression matrix and group information.
#' It could only be used for 2-group analysis.
#' When you run this function, you will see message about "Group 1" and "Group 2".
#' "Group 1" is the reference group, and "Group 2" is the test group.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param group_col The column name of the group information in the sample information.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the fold change results are returned.
#'
#' @return A tibble with two columns: `variable` and `log2fc`.
#' @export
gly_fold_change <- function(exp, group_col = "group", add_info = TRUE) {
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_string(group_col)
  checkmate::check_logical(add_info, len = 1)

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

  # Calculate mean expression for each group
  group1_mean <- rowMeans(expr_mat[, groups == levels(groups)[1], drop = FALSE])
  group2_mean <- rowMeans(expr_mat[, groups == levels(groups)[2], drop = FALSE])

  # Calculate log2 fold change
  log2fc <- log2(group2_mean / group1_mean)

  res <- tibble::tibble(variable = rownames(expr_mat), log2fc = log2fc)

  # Process results with add_info logic
  res <- .process_results_add_info(res, exp, add_info)

  structure(res, class = c("glystats_fc_res", class(res)))
}
