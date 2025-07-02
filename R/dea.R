#' Two-sample t-test for Differential Expression Analysis
#'
#' Perform two-sample t-test for glycomics or glycoproteomics data.
#' The function supports Student's t-test for comparing two groups.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to `stats::t.test()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. Exactly 2 groups are required in the grouping variable.
#'
#' @returns
#' A tibble with t-test results.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#'
#' @export
gly_ttest <- function(exp, group_col = "group", p_adj_method = "BH", return_raw = FALSE, ...) {
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_string(group_col)
  checkmate::check_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::check_logical(return_raw, len = 1)

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

  # Perform t-test
  result <- .gly_dea_2groups(expr_mat, groups, stats::t.test, p_adj_method, return_raw, ...)

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
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to `stats::wilcox.test()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. Exactly 2 groups are required in the grouping variable.
#'
#' @returns
#' A tibble with Wilcoxon test results.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#' 
#' @export
gly_wilcox <- function(exp, group_col = "group", p_adj_method = "BH", return_raw = FALSE, ...) {
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_string(group_col)
  checkmate::check_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::check_logical(return_raw, len = 1)

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

  # Perform Wilcoxon test
  result <- .gly_dea_2groups(expr_mat, groups, stats::wilcox.test, p_adj_method, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class
  structure(result, class = c("glystats_dea_res_wilcoxon", "glystats_dea_res", "glystats_res", class(result)))
}

#' One-way ANOVA for Differential Expression Analysis
#'
#' Perform one-way ANOVA for glycomics or glycoproteomics data.
#' The function supports parametric comparison of multiple groups.
#' For significant results, Tukey's HSD post-hoc test is automatically performed.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to `stats::aov()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. At least 2 groups are required in the grouping variable.
#' 
#' **Post-hoc Test:**
#' Tukey's HSD test for pairwise comparisons (`stats::TukeyHSD()`) is performed
#' for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' A tibble with ANOVA results and a post_hoc column indicating significant group pairs.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#' 
#' @export
gly_anova <- function(exp, group_col = "group", p_adj_method = "BH", return_raw = FALSE, ...) {
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_string(group_col)
  checkmate::check_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::check_logical(return_raw, len = 1)

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

  # Perform ANOVA
  result <- .gly_dea_multi_groups(expr_mat, groups, stats::aov, stats::TukeyHSD, p_adj_method, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class
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
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
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
#' **Post-hoc Test:**
#' Dunn's test with Holm correction for multiple comparisons (`FSA::dunnTest()`) is performed
#' for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' A tibble with Kruskal-Wallis test results and a post_hoc column indicating significant group pairs.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#' 
#' @export
gly_kruskal <- function(exp, group_col = "group", p_adj_method = "BH", return_raw = FALSE, ...) {
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_string(group_col)
  checkmate::check_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::check_logical(return_raw, len = 1)

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

  # Perform Kruskal-Wallis test
  result <- .gly_dea_multi_groups(expr_mat, groups, stats::kruskal.test, FSA::dunnTest, p_adj_method, return_raw, ...)

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class
  structure(result, class = c("glystats_dea_res_kruskal", "glystats_dea_res", "glystats_res", class(result)))
}

.gly_dea_2groups <- function(expr_mat, groups, .f, p_adj_method, return_raw = FALSE, ...) {
  mod_list <- .gly_dea_2groups_raw(expr_mat, groups, .f, ...)
  if (return_raw) {
    return(mod_list)
  }
  .gly_dea_2groups_tibblify(mod_list, .f, p_adj_method)
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
.gly_dea_2groups_tibblify <- function(mod_list, .f, p_adj_method) {
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

  result_tbl
}

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
  main_test_tbl
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
