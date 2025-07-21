#' Two-sample t-test for Differential Expression Analysis
#'
#' Perform two-sample t-test for glycomics or glycoproteomics data.
#' The function supports Student's t-test for comparing two groups.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor vector specifying group membership for each sample.
#'   Must have exactly 2 levels.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param ref_group A character string specifying the reference group.
#'  If NULL (default), the first level of the group factor is used as the reference.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#'  Only applicable to `gly_ttest()`.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to `stats::t.test()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. Exactly 2 groups are required in the grouping variable.
#'
#' `gly_ttest()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_ttest_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @returns A tibble with t-test results including log2 fold change (log2fc),
#'   or a list of `t.test` models if `return_raw` is TRUE.
#' @seealso [stats::t.test()]
#' @export
gly_ttest <- function(
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

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "t-test"
  )
  groups <- group_info$groups
  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Call the underlying API
  result <- gly_ttest_(expr_mat, groups, p_adj_method, ref_group, return_raw, ...)

  # If raw results requested, return directly (no add_info processing needed)
  if (return_raw) {
    return(result)
  }

  # Process results with add_info logic
  result <- .process_results_add_info(result, exp, add_info)

  # Add S3 class
  structure(result, class = c("glystats_dea_res_ttest", "glystats_dea_res", "glystats_res", class(result)))
}

#' @rdname gly_ttest
#' @export
gly_ttest_ <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  ref_group = NULL,
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_logical(return_raw, len = 1)

  # Validate group count
  if (length(levels(groups)) != 2) {
    cli::cli_abort("groups must have exactly 2 levels for t-test")
  }

  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Perform t-test
  result <- .gly_dea_2groups(expr_mat, groups, stats::t.test, p_adj_method, ref_group, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class
  structure(result, class = c("glystats_dea_res_ttest", "glystats_dea_res", "glystats_res", class(result)))
}

#' Wilcoxon rank-sum test for Differential Expression Analysis
#'
#' Perform Wilcoxon rank-sum test (Mann-Whitney U test) for glycomics or glycoproteomics data.
#' The function supports non-parametric comparison of two groups.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor vector specifying group membership for each sample.
#'   Must have exactly 2 levels.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param ref_group A character string specifying the reference group.
#'  If NULL (default), the first level of the group factor is used as the reference.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#'  Only applicable to `gly_wilcox()`.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to `stats::wilcox.test()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. Exactly 2 groups are required in the grouping variable.
#'
#' `gly_wilcox()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_wilcox_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @returns
#' A tibble with Wilcoxon test results including log2 fold change (log2fc),
#' or a list of `wilcox.test` models if `return_raw` is TRUE.
#'
#' @seealso [stats::wilcox.test()]
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#'
#' @export
gly_wilcox <- function(
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

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "wilcoxon"
  )
  groups <- group_info$groups
  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Call the underlying API
  result <- gly_wilcox_(expr_mat, groups, p_adj_method, ref_group, return_raw, ...)

  # If raw results requested, return directly (no add_info processing needed)
  if (return_raw) {
    return(result)
  }

  # Process results with add_info logic
  result <- .process_results_add_info(result, exp, add_info)

  # Add S3 class
  structure(result, class = c("glystats_dea_res_wilcoxon", "glystats_dea_res", "glystats_res", class(result)))
}

#' @rdname gly_wilcox
#' @export
gly_wilcox_ <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  ref_group = NULL,
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_logical(return_raw, len = 1)

  # Validate group count
  if (length(levels(groups)) != 2) {
    cli::cli_abort("groups must have exactly 2 levels for Wilcoxon test")
  }

  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Perform Wilcoxon test
  result <- .gly_dea_2groups(expr_mat, groups, stats::wilcox.test, p_adj_method, ref_group, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class
  structure(result, class = c("glystats_dea_res_wilcoxon", "glystats_dea_res", "glystats_res", class(result)))
}

# Internal helper functions for t-test and Wilcoxon test

.gly_dea_2groups <- function(expr_mat, groups, .f, p_adj_method, ref_group, return_raw = FALSE, ...) {
  # Reorder groups if ref_group is specified
  if (!is.null(ref_group)) {
    groups <- .reorder_groups_for_ref(groups, ref_group)
  }
  cli::cli_alert_info("Reference group: {.val {levels(groups)[1]}}")

  mod_list <- .gly_dea_2groups_raw(expr_mat, groups, .f, ...)
  if (return_raw) {
    return(mod_list)
  }
  .gly_dea_2groups_tibblify(mod_list, .f, p_adj_method, expr_mat, groups)
}

# Generate raw model list for 2-group analysis
.gly_dea_2groups_raw <- function(expr_mat, groups, .f, ...) {
  data <- expr_mat %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(group = groups) %>%
    tidyr::pivot_longer(cols = -all_of(c("sample", "group")), names_to = "variable", values_to = "value") %>%
    dplyr::mutate(log_value = log2(.data$value + 1))

  # Perform statistical tests and store raw results
  nested_data <- data %>%
    dplyr::nest_by(.data$variable) %>%
    dplyr::mutate(test_result = list(.f(log_value ~ group, data = .data$data)))

  # Return named list of raw results
  raw_results <- nested_data$test_result
  names(raw_results) <- nested_data$variable
  raw_results
}

# Convert raw model list to tibble for 2-group analysis
.gly_dea_2groups_tibblify <- function(mod_list, .f, p_adj_method, expr_mat, groups) {
  # Create a tibble from the model list
  var_names <- names(mod_list)

  result_tbl <- tibble::tibble(
    variable = var_names,
    test_result = mod_list
  ) %>%
    dplyr::mutate(
      params = purrr::map(.data$test_result, ~ parameters::model_parameters(.x)),
      params = purrr::map(.data$params, ~ parameters::standardize_names(.x)),
    ) %>%
    dplyr::select(all_of(c("variable", "params"))) %>%
    tidyr::unnest(all_of("params")) %>%
    dplyr::ungroup() %>%
    janitor::clean_names()

  if (!is.null(p_adj_method)) {
    result_tbl <- dplyr::mutate(result_tbl, p_adj = stats::p.adjust(.data$p, method = p_adj_method))
  }

  # Calculate log2 fold change
  result_tbl <- .add_log2fc_to_result(result_tbl, expr_mat, groups)

  result_tbl
}

# Helper function to add log2 fold change to DEA results
.add_log2fc_to_result <- function(result_tbl, expr_mat, groups) {
  # Calculate mean expression for each group for each variable
  group_levels <- levels(groups)
  group1_indices <- which(groups == group_levels[1])
  group2_indices <- which(groups == group_levels[2])

  # Calculate log2 fold change for each variable
  log2fc_values <- purrr::map_dbl(result_tbl$variable, function(var_name) {
    # Get expression values for this variable
    var_expr <- expr_mat[var_name, ]

    # Calculate mean expression for each group
    group1_mean <- mean(var_expr[group1_indices], na.rm = TRUE)
    group2_mean <- mean(var_expr[group2_indices], na.rm = TRUE)

    # Calculate log2 fold change (group2 vs group1)
    log2(group2_mean / group1_mean)
  })

  # Add log2fc column to result
  result_tbl$log2fc <- log2fc_values

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