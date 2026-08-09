#' Two-sample t-test for Differential Expression Analysis
#'
#' Perform two-sample t-test for glycomics or glycoproteomics data.
#' The function supports Student's t-test for comparing two groups.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing an expression matrix and
#'   sample information.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param ref_group A character string specifying the reference group.
#'  If NULL (default), the first level of the group factor is used as the reference.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#' @param ... Additional arguments passed to `stats::t.test()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1e-6)) before
#' statistical testing. Exactly 2 groups are required in the grouping variable.
#'
#' @returns A list with three elements:
#' - `tidy_result`: A tibble with t-test results containing the following columns:
#'   - `variable`: Variable name
#'   - `estimate`: Difference in group means (group2 - group1)
#'   - `estimate1`: Mean of group 1
#'   - `estimate2`: Mean of group 2
#'   - `statistic`: t-statistic
#'   - `p_val`: Raw p-value from t-test
#'   - `parameter`: Degrees of freedom
#'   - `conf_low`: Lower bound of 95% confidence interval
#'   - `conf_high`: Upper bound of 95% confidence interval
#'   - `effect_size`: Cohen's d
#'   - `method`: Statistical method used
#'   - `alternative`: Alternative hypothesis
#'   - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'   - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean))
#' - `raw_result`: A list of `t.test` model objects
#' - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_ttest_res` and `glystats_res`.
#' @seealso [stats::t.test()]
#' @export
gly_ttest <- function(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  ref_group = NULL,
  add_info = TRUE,
  ...
) {
  # Validate inputs
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)

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

  # Run the internal computation
  result <- .analyze_ttest(expr_mat, groups, p_adj_method, ref_group, ...)

  # Process results with add_info logic
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_ttest()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_ttest <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  ref_group = NULL,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )

  # Validate group count
  if (length(levels(groups)) != 2) {
    cli::cli_abort("groups must have exactly 2 levels for t-test")
  }

  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Perform t-test
  result <- .gly_dea_2groups(
    expr_mat,
    groups,
    stats::t.test,
    p_adj_method,
    ref_group,
    ...
  )

  # Add S3 class
  structure(result, class = c("glystats_ttest_res", "glystats_res"))
}

#' Wilcoxon rank-sum test for Differential Expression Analysis
#'
#' Perform Wilcoxon rank-sum test (Mann-Whitney U test) for glycomics or glycoproteomics data.
#' The function supports non-parametric comparison of two groups.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing an expression matrix and
#'   sample information.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param ref_group A character string specifying the reference group.
#'  If NULL (default), the first level of the group factor is used as the reference.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#' @param ... Additional arguments passed to `stats::wilcox.test()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1e-6)) before
#' statistical testing. Exactly 2 groups are required in the grouping variable.
#'
#' @returns
#' A list with three elements:
#' - `tidy_result`: A tibble with Wilcoxon test results containing the following columns:
#'   - `variable`: Variable name
#'   - `statistic`: Wilcoxon test statistic
#'   - `p_val`: Raw p-value from Wilcoxon test
#'   - `effect_size`: Rank-biserial correlation
#'   - `method`: Statistical method used
#'   - `alternative`: Alternative hypothesis
#'   - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'   - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean))
#'   Additional columns from experiment metadata may be included if add_info = TRUE.
#' - `raw_result`: A list of `wilcox.test` model objects
#' - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_wilcox_res` and `glystats_res`.
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
  ...
) {
  # Validate inputs
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)

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

  # Run the internal computation
  result <- .analyze_wilcox(expr_mat, groups, p_adj_method, ref_group, ...)

  # Process results with add_info logic
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_wilcox()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_wilcox <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  ref_group = NULL,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )

  # Validate group count
  if (length(levels(groups)) != 2) {
    cli::cli_abort("groups must have exactly 2 levels for Wilcoxon test")
  }

  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Perform Wilcoxon test
  result <- .gly_dea_2groups(
    expr_mat,
    groups,
    stats::wilcox.test,
    p_adj_method,
    ref_group,
    ...
  )

  # Add S3 class
  structure(result, class = c("glystats_wilcox_res", "glystats_res"))
}

# Internal helper functions for t-test and Wilcoxon test

.gly_dea_2groups <- function(
  expr_mat,
  groups,
  .f,
  p_adj_method,
  ref_group,
  ...
) {
  # Reorder groups if ref_group is specified
  if (!is.null(ref_group)) {
    groups <- .reorder_groups_for_ref(groups, ref_group)
  }

  log_expr_mat <- .log_transform_expr_mat(expr_mat)

  mod_list <- .gly_dea_2groups_raw(log_expr_mat, groups, .f, ...)
  tidy_result <- .gly_dea_2groups_tibblify(
    mod_list,
    .f,
    p_adj_method,
    expr_mat,
    groups,
    log_expr_mat
  )

  # Return list with both tidy and raw results
  list(
    tidy_result = tidy_result,
    raw_result = mod_list
  )
}

# Generate raw model list for 2-group analysis
.gly_dea_2groups_raw <- function(log_expr_mat, groups, .f, ...) {
  data <- log_expr_mat %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(group = groups) %>%
    tidyr::pivot_longer(
      cols = -all_of(c("sample", "group")),
      names_to = "variable",
      values_to = "log_value"
    )

  # Perform statistical tests and store raw results
  dots <- rlang::list2(...)
  disallowed_args <- intersect(names(dots), c("formula", "data"))
  if (length(disallowed_args) > 0) {
    cli::cli_abort(
      "Arguments {cli::format_inline(disallowed_args)} should not be supplied through `...`; they are controlled internally."
    )
  }
  safe_f <- purrr::possibly(function(...) rlang::exec(.f, ...), otherwise = NA)
  group_levels <- levels(groups)
  nested_data <- data %>%
    dplyr::nest_by(.data$variable)

  # Run tests as x = group2, y = group1 so direction-sensitive outputs align
  # with the package convention used by log2fc.
  raw_results <- purrr::map(nested_data$data, function(df) {
    safe_f(
      x = df$log_value[df$group == group_levels[2]],
      y = df$log_value[df$group == group_levels[1]],
      !!!dots
    )
  })
  names(raw_results) <- nested_data$variable
  n_na <- sum(is.na(raw_results))
  if (n_na > 0) {
    cli::cli_warn("{.val {n_na}} variables failed to fit the model")
  }
  raw_results
}

# Convert raw model list to tibble for 2-group analysis
.gly_dea_2groups_tibblify <- function(
  mod_list,
  .f,
  p_adj_method,
  expr_mat,
  groups,
  log_expr_mat
) {
  # Create a tibble from the model list
  var_names <- names(mod_list)
  all_failed <- all(purrr::map_lgl(mod_list, rlang::is_na))

  if (all_failed && identical(.f, stats::t.test)) {
    result_tbl <- tibble::tibble(
      variable = var_names,
      estimate = NA_real_,
      estimate1 = NA_real_,
      estimate2 = NA_real_,
      statistic = NA_real_,
      p_val = NA_real_,
      parameter = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      method = NA_character_,
      alternative = NA_character_
    )
  } else if (all_failed && identical(.f, stats::wilcox.test)) {
    result_tbl <- tibble::tibble(
      variable = var_names,
      statistic = NA_real_,
      p_val = NA_real_,
      method = NA_character_,
      alternative = NA_character_
    )
  } else {
    result_tbl <- tibble::tibble(
      variable = var_names,
      test_result = mod_list
    ) %>%
      dplyr::mutate(
        params = purrr::map(
          .data$test_result,
          ~ {
            if (rlang::is_na(.x)) NULL else broom::tidy(.x)
          }
        ),
      ) %>%
      dplyr::select(all_of(c("variable", "params"))) %>%
      tidyr::unnest(all_of("params"), keep_empty = TRUE) %>%
      dplyr::ungroup() %>%
      janitor::clean_names()

    # Rename p_value to p_val for consistency
    if ("p_value" %in% colnames(result_tbl)) {
      result_tbl <- dplyr::rename(result_tbl, p_val = "p_value")
    }
  }

  if (!is.null(p_adj_method)) {
    result_tbl <- dplyr::mutate(
      result_tbl,
      p_adj = stats::p.adjust(.data$p_val, method = p_adj_method)
    )
  }

  # Calculate log2 fold change
  result_tbl <- .add_log2fc_to_result(result_tbl, expr_mat, groups)
  result_tbl <- .add_effect_size_to_2group_result(
    result_tbl,
    log_expr_mat,
    groups,
    .f
  )

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

    group1_values <- var_expr[group1_indices]
    group2_values <- var_expr[group2_indices]
    if (all(is.na(group1_values)) || all(is.na(group2_values))) {
      return(NA_real_)
    }

    # Calculate mean expression for each group
    group1_mean <- mean(group1_values, na.rm = TRUE)
    group2_mean <- mean(group2_values, na.rm = TRUE)

    # Calculate log2 fold change (group2 vs group1)
    log2(group2_mean / group1_mean)
  })

  # Add log2fc column to result
  result_tbl$log2fc <- log2fc_values

  result_tbl
}

#' Add effect sizes to two-group differential analysis results
#'
#' @param result_tbl A tibble containing the tidy test results.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor specifying group membership for each sample.
#' @param .f The statistical test function used to generate the results.
#'
#' @return A tibble with an `effect_size` column added.
#' @noRd
.add_effect_size_to_2group_result <- function(
  result_tbl,
  expr_mat,
  groups,
  .f
) {
  effect_size_values <- purrr::map_dbl(result_tbl$variable, function(var_name) {
    if (identical(.f, stats::t.test)) {
      .calculate_cohens_d(expr_mat, groups, var_name)
    } else if (identical(.f, stats::wilcox.test)) {
      .calculate_rank_biserial(expr_mat, groups, var_name)
    } else {
      NA_real_
    }
  })

  result_tbl$effect_size <- effect_size_values
  result_tbl
}

#' Calculate Cohen's d for a two-group comparison
#'
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor specifying group membership for each sample.
#' @param var_name A character string specifying the variable to evaluate.
#'
#' @return A numeric scalar containing Cohen's d.
#' @noRd
.calculate_cohens_d <- function(expr_mat, groups, var_name) {
  log_values <- expr_mat[var_name, ]
  group_levels <- levels(groups)
  ref_values <- log_values[groups == group_levels[1]]
  test_values <- log_values[groups == group_levels[2]]
  ref_values <- ref_values[!is.na(ref_values)]
  test_values <- test_values[!is.na(test_values)]

  if (length(ref_values) < 2 || length(test_values) < 2) {
    return(NA_real_)
  }

  pooled_sd <- sqrt(
    (((length(ref_values) - 1) * stats::sd(ref_values)^2) +
      ((length(test_values) - 1) * stats::sd(test_values)^2)) /
      (length(ref_values) + length(test_values) - 2)
  )

  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }

  (mean(test_values) - mean(ref_values)) / pooled_sd
}

#' Calculate rank-biserial correlation for a two-group comparison
#'
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor specifying group membership for each sample.
#' @param var_name A character string specifying the variable to evaluate.
#'
#' @return A numeric scalar containing the rank-biserial correlation.
#' @noRd
.calculate_rank_biserial <- function(expr_mat, groups, var_name) {
  log_values <- expr_mat[var_name, ]
  group_levels <- levels(groups)
  ref_values <- log_values[groups == group_levels[1]]
  test_values <- log_values[groups == group_levels[2]]
  ref_values <- ref_values[!is.na(ref_values)]
  test_values <- test_values[!is.na(test_values)]

  if (length(ref_values) == 0 || length(test_values) == 0) {
    return(NA_real_)
  }

  combined_values <- c(ref_values, test_values)
  test_ranks <- rank(combined_values)[
    (length(ref_values) + 1):length(combined_values)
  ]
  test_size <- length(test_values)
  ref_size <- length(ref_values)
  u_stat <- sum(test_ranks) - test_size * (test_size + 1) / 2

  (2 * u_stat / (ref_size * test_size)) - 1
}

#' Log-transform an expression matrix
#'
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#'
#' @return A numeric matrix with log2(x + 1e-6) values.
#' @noRd
.log_transform_expr_mat <- function(expr_mat) {
  log2(expr_mat + 1e-6)
}

# Helper function to reorder groups so that ref_group becomes the first level
.reorder_groups_for_ref <- function(groups, ref_group) {
  current_levels <- levels(groups)

  # Move ref_group to the first position
  new_levels <- c(ref_group, setdiff(current_levels, ref_group))

  # Reorder the factor levels
  factor(groups, levels = new_levels)
}
