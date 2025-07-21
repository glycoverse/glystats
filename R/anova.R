#' One-way ANOVA for Differential Expression Analysis
#'
#' Perform one-way ANOVA for glycomics or glycoproteomics data.
#' The function supports parametric comparison of multiple groups.
#' For significant results, Tukey's HSD post-hoc test is automatically performed.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor vector specifying group membership for each sample.
#'   Must have at least 2 levels.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#'  Only applicable to `gly_anova()`.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to `stats::aov()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. At least 2 groups are required in the grouping variable.
#'
#' `gly_anova()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_anova_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' **Post-hoc Test:**
#' Tukey's HSD test for pairwise comparisons (`stats::TukeyHSD()`) is performed
#' for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' A list containing two tibbles: `main_test` with ANOVA results and a post_hoc column
#' indicating significant group pairs, and `post_hoc_test` with detailed pairwise
#' comparison results in long format (variable, group1, group2, p_value, p_adj columns).
#' If `return_raw` is TRUE, returns a list of `aov` models and `TukeyHSD` objects.
#'
#' @seealso [stats::aov()], [stats::TukeyHSD()]
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#'
#' @export
gly_anova <- function(exp, group_col = "group", p_adj_method = "BH", add_info = TRUE, return_raw = FALSE, ...) {
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
    max_count = NULL,
    method = "anova"
  )
  groups <- group_info$groups

  # Call the underlying API
  result <- gly_anova_(expr_mat, groups, p_adj_method, return_raw, ...)

  # If raw results requested, return directly (no add_info processing needed)
  if (return_raw) {
    return(result)
  }

  # Process results with add_info logic for main_test
  result$main_test <- .process_results_add_info(result$main_test, exp, add_info)

  # Add S3 class to the entire list
  structure(result, class = c("glystats_dea_res_anova", "glystats_dea_res", "glystats_res", class(result)))
}

#' @rdname gly_anova
#' @export
gly_anova_ <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_logical(return_raw, len = 1)

  # Validate group count
  if (length(levels(groups)) < 2) {
    cli::cli_abort("groups must have at least 2 levels for ANOVA")
  }

  # Perform ANOVA
  result <- .gly_dea_multi_groups(expr_mat, groups, stats::aov, stats::TukeyHSD, p_adj_method, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class to the entire list
  structure(result, class = c("glystats_dea_res_anova", "glystats_dea_res", "glystats_res", class(result)))
}

#' Kruskal-Wallis test for Differential Expression Analysis
#'
#' Perform Kruskal-Wallis test for glycomics or glycoproteomics data.
#' The function supports non-parametric comparison of multiple groups.
#' For significant results, Dunn's post-hoc test is automatically performed.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor vector specifying group membership for each sample.
#'   Must have at least 2 levels.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#'  Only applicable to `gly_kruskal()`.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to `stats::kruskal.test()`.
#'
#' @section Required packages:
#' This function requires the `FSA` package for Dunn's post-hoc test.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. At least 2 groups are required in the grouping variable.
#'
#' `gly_kruskal()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_kruskal_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' **Post-hoc Test:**
#' Dunn's test with Holm correction for multiple comparisons (`FSA::dunnTest()`) is performed
#' for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' A list containing two tibbles: `main_test` with Kruskal-Wallis test results and a post_hoc column
#' indicating significant group pairs, and `post_hoc_test` with detailed pairwise
#' comparison results in long format (variable, group1, group2, p_value, p_adj columns).
#' If `return_raw` is TRUE, returns a list of `kruskal.test` objects and `dunnTest` objects.
#'
#' @seealso [stats::kruskal.test()], [FSA::dunnTest()]
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#'
#' @export
gly_kruskal <- function(exp, group_col = "group", p_adj_method = "BH", add_info = TRUE, return_raw = FALSE, ...) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

  # Check package availability
  .check_pkg_available("FSA")

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "kruskal"
  )
  groups <- group_info$groups

  # Call the underlying API
  result <- gly_kruskal_(expr_mat, groups, p_adj_method, return_raw, ...)

  # If raw results requested, return directly (no add_info processing needed)
  if (return_raw) {
    return(result)
  }

  # Process results with add_info logic for main_test
  result$main_test <- .process_results_add_info(result$main_test, exp, add_info)

  # Add S3 class to the entire list
  structure(result, class = c("glystats_dea_res_kruskal", "glystats_dea_res", "glystats_res", class(result)))
}

#' @rdname gly_kruskal
#' @export
gly_kruskal_ <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  return_raw = FALSE,
  ...
) {
  # Check package availability
  .check_pkg_available("FSA")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_logical(return_raw, len = 1)

  # Validate group count
  if (length(levels(groups)) < 2) {
    cli::cli_abort("groups must have at least 2 levels for Kruskal-Wallis test")
  }

  # Perform Kruskal-Wallis test
  result <- .gly_dea_multi_groups(expr_mat, groups, stats::kruskal.test, FSA::dunnTest, p_adj_method, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class to the entire list
  structure(result, class = c("glystats_dea_res_kruskal", "glystats_dea_res", "glystats_res", class(result)))
}

# Internal helper functions for multi-group analysis (ANOVA and Kruskal-Wallis)

# Multi-group analysis function
.gly_dea_multi_groups <- function(expr_mat, groups, .f, .ph, p_adj_method, return_raw = FALSE, ...) {
  if (length(levels(groups)) < 2) {
    cli::cli_abort("Multi-group analysis requires at least 2 groups")
  }
  mod_list <- .gly_dea_multi_groups_raw(expr_mat, groups, .f, p_adj_method, ...)
  if (return_raw) {
    return(mod_list)
  }
  .gly_dea_multi_groups_tibblify(mod_list, .f, p_adj_method)
}

# Helper function to perform raw post-hoc tests (returns raw objects)
.perform_raw_posthoc_test <- function(data_nested, .f) {
  if (identical(.f, stats::aov)) {
    # TukeyHSD for ANOVA - return raw TukeyHSD object
    aov_model <- stats::aov(log_value ~ group, data = data_nested)
    stats::TukeyHSD(aov_model)
  } else {
    # Dunn test for Kruskal-Wallis - return raw dunnTest object
    FSA::dunnTest(log_value ~ group, data = data_nested, method = "holm")
  }
}

# Helper function to generate raw main test results
.generate_raw_main_results <- function(data, .f, ...) {
  main_test_raw <- data %>%
    dplyr::nest_by(.data$variable) %>%
    dplyr::mutate(
      test_result = list(.f(log_value ~ group, data = .data$data))
    )

  main_test_list <- main_test_raw$test_result
  names(main_test_list) <- main_test_raw$variable
  main_test_list
}

# Helper function to generate raw post-hoc results
.generate_raw_posthoc_results <- function(main_test_raw, data, .f, p_adj_method) {
  # First, we need to determine which variables are significant
  # We'll need to extract p-values from the raw results
  significant_vars <- c()

  for (var_name in names(main_test_raw)) {
    raw_result <- main_test_raw[[var_name]]

    # Extract p-value from raw result
    if (identical(.f, stats::aov)) {
      # For ANOVA, extract p-value from summary
      p_val <- summary(raw_result)[[1]][["Pr(>F)"]][1]
    } else {
      # For Kruskal-Wallis, p-value is in $p.value
      p_val <- raw_result$p.value
    }

    # Apply p-adjustment if needed
    if (!is.null(p_val) && !is.na(p_val)) {
      if (p_val < 0.05) {  # Use unadjusted p-value for initial filtering
        significant_vars <- c(significant_vars, var_name)
      }
    }
  }

  # Apply p-adjustment to all p-values if requested
  if (!is.null(p_adj_method) && length(significant_vars) > 0) {
    all_p_vals <- sapply(names(main_test_raw), function(var_name) {
      raw_result <- main_test_raw[[var_name]]
      if (identical(.f, stats::aov)) {
        summary(raw_result)[[1]][["Pr(>F)"]][1]
      } else {
        raw_result$p.value
      }
    })

    adj_p_vals <- stats::p.adjust(all_p_vals, method = p_adj_method)
    significant_vars <- names(adj_p_vals)[adj_p_vals < 0.05 & !is.na(adj_p_vals)]
  }

  if (length(significant_vars) > 0) {
    # Perform post-hoc tests for significant variables and return raw objects
    posthoc_raw_results <- data %>%
      dplyr::filter(.data$variable %in% significant_vars) %>%
      dplyr::nest_by(.data$variable) %>%
      dplyr::mutate(posthoc_raw = list(.perform_raw_posthoc_test(.data$data, .f))) %>%
      dplyr::select(all_of(c("variable", "posthoc_raw")))

    # Convert to named list
    raw_posthoc_list <- posthoc_raw_results$posthoc_raw
    names(raw_posthoc_list) <- posthoc_raw_results$variable

    return(raw_posthoc_list)
  } else {
    # No significant results, return empty named list
    return(list())
  }
}

# Helper function to prepare data for analysis
.prepare_multi_group_data <- function(expr_mat, groups) {
  expr_mat %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(group = groups) %>%
    tidyr::pivot_longer(
      cols = -all_of(c("sample", "group")),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::mutate(log_value = log2(.data$value + 1))
}

# Generate raw model list for multi-group analysis
.gly_dea_multi_groups_raw <- function(expr_mat, groups, .f, p_adj_method, ...) {
  data <- .prepare_multi_group_data(expr_mat, groups)
  main_test_list <- .generate_raw_main_results(data, .f, ...)
  posthoc_raw <- .generate_raw_posthoc_results(main_test_list, data, .f, p_adj_method)

  list(
    main_test = main_test_list,
    post_hoc = posthoc_raw
  )
}

# Convert raw model list to tibble for multi-group analysis
.gly_dea_multi_groups_tibblify <- function(mod_list, .f, p_adj_method) {
  main_test_tbl <- .tibblify_main_test_results(mod_list$main_test, .f, p_adj_method)
  post_hoc_vec <- .format_posthoc_results(mod_list$post_hoc, .f, main_test_tbl$variable)
  main_test_tbl$post_hoc <- post_hoc_vec

  # Create post-hoc test tibble in long format
  post_hoc_test_tbl <- .tibblify_posthoc_results(mod_list$post_hoc, .f)

  # Return list with both tibbles
  list(
    main_test = main_test_tbl,
    post_hoc_test = post_hoc_test_tbl
  )
}

# Helper function to convert main test raw results to tibble
.tibblify_main_test_results <- function(main_test_raw, .f, p_adj_method) {
  var_names <- names(main_test_raw)

  result_tbl <- tibble::tibble(
    variable = var_names,
    test_result = main_test_raw
  ) %>%
    dplyr::mutate(
      params = purrr::map(.data$test_result, ~ parameters::model_parameters(.x)),
      params = purrr::map(.data$params, ~ parameters::standardize_names(.x)),
    ) %>%
    dplyr::select(all_of(c("variable", "params"))) %>%
    tidyr::unnest(all_of("params")) %>%
    dplyr::ungroup() %>%
    janitor::clean_names()

  # For ANOVA, filter to keep only group effects (not residuals)
  # For Kruskal-Wallis, keep all results (there's only one row per variable anyway)
  if (identical(.f, stats::aov)) {
    result_tbl <- result_tbl %>%
      dplyr::filter(.data$parameter == "group")
  }

  if (!is.null(p_adj_method)) {
    result_tbl <- dplyr::mutate(result_tbl, p_adj = stats::p.adjust(.data$p, method = p_adj_method))
  }

  result_tbl
}

# Helper function to convert post-hoc raw results to tibble
.format_posthoc_results <- function(posthoc_raw, .f, variables) {
  if (length(posthoc_raw) == 0) {
    return(rep(NA_character_, length(variables)))
  }
  posthoc_map <- purrr::imap(posthoc_raw, function(raw_result, var_name) {
    if (identical(.f, stats::aov)) {
      tukey_df <- as.data.frame(raw_result$group)
      sig_pairs <- rownames(tukey_df)[tukey_df$"p adj" < 0.05]
      sig_str <- if (length(sig_pairs) == 0) NA_character_ else paste(sig_pairs, collapse = ";")
    } else {
      dunn_df <- raw_result$res
      sig_pairs <- dunn_df$Comparison[dunn_df$P.adj < 0.05]
      sig_str <- if (length(sig_pairs) == 0) NA_character_ else paste(sig_pairs, collapse = ";")
    }
    sig_str
  })
  purrr::map_chr(variables, ~ posthoc_map[[.x]] %||% NA_character_)
}

# Helper function to convert post-hoc raw results to long format tibble
.tibblify_posthoc_results <- function(posthoc_raw, .f) {
  if (length(posthoc_raw) == 0) {
    # Return empty tibble with correct structure
    return(tibble::tibble(
      variable = character(0),
      group1 = character(0),
      group2 = character(0),
      p_value = numeric(0),
      p_adj = numeric(0)
    ))
  }

  # Convert each post-hoc result to tibble and combine
  result_list <- purrr::imap(posthoc_raw, function(raw_result, var_name) {
    if (identical(.f, stats::aov)) {
      # For TukeyHSD results
      tukey_df <- as.data.frame(raw_result$group)
      tukey_df$comparison <- rownames(tukey_df)

      # Parse group comparisons (format: "group2-group1")
      comparison_parts <- stringr::str_split(tukey_df$comparison, "-", simplify = TRUE)

      tibble::tibble(
        variable = var_name,
        group1 = comparison_parts[, 2],  # Second part is group1
        group2 = comparison_parts[, 1],  # First part is group2
        p_value = tukey_df$`p adj`,      # TukeyHSD already provides adjusted p-values
        p_adj = tukey_df$`p adj`
      )
    } else {
      # For Dunn test results
      dunn_df <- raw_result$res

      # Parse group comparisons (format: "group1 - group2")
      comparison_parts <- stringr::str_split(dunn_df$Comparison, " - ", simplify = TRUE)

      tibble::tibble(
        variable = var_name,
        group1 = comparison_parts[, 1],
        group2 = comparison_parts[, 2],
        p_value = dunn_df$P.unadj,       # Unadjusted p-value
        p_adj = dunn_df$P.adj            # Adjusted p-value
      )
    }
  })

  # Combine all results
  dplyr::bind_rows(result_list)
}
