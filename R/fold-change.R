#' Calculate fold change
#'
#' Calculate fold change for a given expression matrix and group information.
#' It could only be used for 2-group analysis.
#' When you run this function, you will see message about "Ref Group" and "Test Group".
#' "Ref Group" is the reference group, and "Test Group" is the test/treatment/case group.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param expr_mat (Only for [gly_fold_change_()]) A numeric matrix with variables as rows and samples as columns.
#' @param groups (Only for [gly_fold_change_()]) A factor or character vector specifying group membership for each sample.
#'   Character vectors will be automatically converted to factors.
#'   If two groups, the first level is the reference group.
#'   If more than two groups, pairwise comparisons will be performed,
#'   with levels coming first as reference groups.
#' @param group_col (Only for [gly_fold_change()]) The column name of the group information in the sample information.
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
#' @returns A tibble with fold change results containing the following columns:
#'   - `variable`: Variable name
#'   - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean))
#'   If more than two groups, two additional columns `ref_group` and `test_group` will be added.
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
    method = "fold_change",
    show_info = TRUE
  )
  groups <- group_info$groups

  expr_mat <- exp$expr_mat

  # Call the underlying API
  result <- gly_fold_change_(expr_mat, groups)

  # Process results with add_info logic
  .process_results_add_info(result, exp, add_info)
}

#' @rdname gly_fold_change
#' @export
gly_fold_change_ <- function(expr_mat, groups) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  groups <- .convert_groups_to_factor(groups)
  checkmate::assert_factor(groups, len = ncol(expr_mat))

  # Calculate fold change
  if (length(levels(groups)) == 2) {
    result <- .fc_2groups(expr_mat, groups)
  } else {
    result <- .fc_multi_groups(expr_mat, groups)
  }
  structure(result, class = c("glystats_fc_res", "glystats_res", class(result)))
}

#' Calculate fold change for 2-group analysis
#'
#' @param expr_mat (Only for `.fc_2groups()`) A numeric matrix with variables as rows and samples as columns.
#' @param groups (Only for `.fc_2groups()`) A factor vector specifying group membership for each sample.
#'   Must have exactly 2 levels. The first level is the reference group.
#'
#' @return A tibble with fold change results containing the following columns:
#'   - `variable`: Variable name
#'   - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean))
#' @noRd
.fc_2groups <- function(expr_mat, groups) {
  group1_mean <- rowMeans(expr_mat[, groups == levels(groups)[1], drop = FALSE])
  group2_mean <- rowMeans(expr_mat[, groups == levels(groups)[2], drop = FALSE])
  log2fc <- log2(group2_mean / group1_mean)
  tibble::tibble(variable = rownames(expr_mat), log2fc = log2fc)
}

#' Calculate fold change for multi-group analysis
#'
#' @param expr_mat (Only for `.fc_multi_groups()`) A numeric matrix with variables as rows and samples as columns.
#' @param groups (Only for `.fc_multi_groups()`) A factor vector specifying group membership for each sample.
#'   Must have more than 2 levels.
#'
#' @return A tibble with fold change results containing the following columns:
#'   - `variable`: Variable name
#'   - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean))
#'   - `ref_group`: Reference group
#'   - `test_group`: Test/treatment/case group
#' @noRd
.fc_multi_groups <- function(expr_mat, groups) {
  comparisons <- .make_comparisons(levels(groups))
  result_list <- purrr::map(comparisons, function(c) {
    sample_mask <- groups %in% c
    sub_groups <- factor(groups[sample_mask], levels = c)
    .fc_2groups(expr_mat[, sample_mask, drop = FALSE], sub_groups)
  })
  comparison_str <- purrr::map_chr(comparisons, ~ stringr::str_c(.x[1], "_vs_", .x[2]))
  names(result_list) <- comparison_str
  dplyr::bind_rows(result_list, .id = "comparison") |>
    tidyr::separate(col = "comparison", into = c("ref_group", "test_group"), sep = "_vs_")
}
