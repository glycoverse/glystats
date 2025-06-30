#' Calculate fold change
#'
#' Calculate fold change for a given expression matrix and group information.
#' It could only be used for 2-group analysis.
#' When you run this function, you will see message about "Group 1" and "Group 2".
#' "Group 1" is the reference group, and "Group 2" is the test group.
#' 
#' @param exp A `glystats_exp` object.
#' @param group_col The column name of the group information in the `exp$var_info` data frame.
#'
#' @return A tibble with two columns: `variable` and `log2fc`.
#' @export
gly_fold_change <- function(exp, group_col = "group") {
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_string(group_col)

  if (!group_col %in% colnames(exp$sample_info)) {
    cli::cli_abort("Column {.field {group_col}} not found in sample information")
  }

  groups <- exp$sample_info[[group_col]]
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }

  if (length(levels(groups)) != 2) {
    cli::cli_abort("Group column {.field {group_col}} must have exactly 2 levels")
  }

  cli::cli_alert_info("Group 1: {.val {levels(groups)[1]}}")
  cli::cli_alert_info("Group 2: {.val {levels(groups)[2]}}")

  expr_mat <- exp$expr_mat
  
  # Calculate mean expression for each group
  group1_mean <- rowMeans(expr_mat[, groups == levels(groups)[1], drop = FALSE])
  group2_mean <- rowMeans(expr_mat[, groups == levels(groups)[2], drop = FALSE])
  
  # Calculate log2 fold change
  log2fc <- log2(group2_mean / group1_mean)
  
  res <- tibble::tibble(variable = rownames(expr_mat), log2fc = log2fc)
  structure(res, class = c("glystats_fc_res", class(res)))
}
