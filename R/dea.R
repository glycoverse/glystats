#' Differential Expression Analysis (DEA)
#'
#' Perform differential expression analysis for glycomics or glycoproteomics data
#' using various statistical methods. The function supports both two-group comparisons
#' (t-test, Wilcoxon test) and multi-group comparisons (ANOVA, Kruskal-Wallis test).
#' For multi-group methods, post-hoc tests are automatically performed for significant results.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param method A character string specifying the statistical method to use:
#'   - `"t-test"`: Student's t-test for two-group comparison (parametric)
#'   - `"wilcoxon"`: Wilcoxon rank-sum test for two-group comparison (non-parametric)
#'   - `"anova"`: One-way ANOVA for multi-group comparison (parametric)
#'   - `"kruskal"`: Kruskal-Wallis test for multi-group comparison (non-parametric)
#'  If not provided, default method is t-test for 2 groups, anova for more than 2 groups.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw statistical model objects as a list.
#' @param ... Additional arguments passed to the underlying statistical functions.
#'
#' @section Required packages:
#' Depending on the method used, this function may require additional packages:
#' - `FSA` package is required when using the "kruskal" method for Dunn's post-hoc test
#' - All other methods use base R statistical functions
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. For two-group methods (t-test, Wilcoxon), exactly 2 groups are required.
#' For multi-group methods (ANOVA, Kruskal-Wallis), at least 2 groups are required.
#'
#' **Underlying Statistical Functions:**
#' - `"t-test"`: `stats::t.test()`
#' - `"wilcoxon"`: `stats::wilcox.test()`
#' - `"anova"`: `stats::aov()`
#' - `"kruskal"`: `stats::kruskal.test()`
#'
#' **Post-hoc Tests:**
#' - For ANOVA: Tukey's HSD test for pairwise comparisons (`stats::TukeyHSD()`)
#' - For Kruskal-Wallis: Dunn's test with Holm correction for multiple comparisons (`FSA::dunnTest()`)
#'
#' Post-hoc tests are only performed for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' For two-group methods (`"t-test"`, `"wilcoxon"`):
#' A tibble with statistical test results.
#'
#' For multi-group methods (`"anova"`, `"kruskal"`):
#' A list containing two elements:
#' - `main_test`: A tibble with main statistical test results
#' - `post_hoc`: A tibble with post-hoc pairwise comparison results (empty if no significant results)
#'
#' @examples
#' \dontrun{
#' # Create example experiment object
#' sample_info <- tibble::tibble(
#'   sample = paste0("S", 1:6),
#'   group = rep(c("Control", "Treatment"), each = 3)
#' )
#' 
#' var_info <- tibble::tibble(
#'   variable = paste0("Glycan_", 1:100)
#' )
#' 
#' expr_mat <- matrix(rnorm(600), nrow = 100, ncol = 6)
#' rownames(expr_mat) <- var_info$variable
#' colnames(expr_mat) <- sample_info$sample
#' 
#' exp <- glyexp::experiment(
#'   expr_mat = expr_mat,
#'   sample_info = sample_info,
#'   var_info = var_info,
#'   exp_type = "glycomics",
#'   glycan_type = "N"
#' )
#' 
#' # Two-group comparison with t-test
#' result_ttest <- gly_dea(exp, method = "t-test")
#' 
#' # Two-group comparison with Wilcoxon test
#' result_wilcox <- gly_dea(exp, method = "wilcoxon")
#' 
#' # Multi-group comparison with ANOVA
#' sample_info$group <- rep(c("A", "B", "C"), each = 2)
#' exp_multi <- glyexp::experiment(expr_mat, sample_info, var_info, "glycomics", "N")
#' result_anova <- gly_dea(exp_multi, method = "anova")
#' 
#' # Access main test results
#' result_anova$main_test
#' 
#' # Access post-hoc results
#' result_anova$post_hoc
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#' 
#' @export
gly_dea <- function(exp, method = NULL, group_col = "group", p_adj_method = "BH", return_raw = FALSE, ...) {
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_choice(method, c("t-test", "wilcoxon", "anova", "kruskal"), null.ok = TRUE)
  checkmate::check_string(group_col)
  checkmate::check_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::check_logical(return_raw, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Check if group column exists
  if (!group_col %in% colnames(sample_info)) {
    cli::cli_abort("Column {.field {group_col}} not found in sample information")
  }

  # Get group information
  groups <- sample_info[[group_col]]
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }

  # Check method
  n_groups <- length(levels(groups))
  if (is.null(method)) {
    # Default method is t-test for 2 groups, anova for more than 2 groups
    method <- if (n_groups == 2) "t-test" else "anova"
  } else if (method %in% c("t-test", "wilcoxon")) {
    # Or check user-specified method
    if (n_groups != 2) {
      cli::cli_abort(c(
        "{.field {group_col}} must be a factor with exactly 2 levels for {.val {method}}",
        "i" = "Current levels: {.val {levels(groups)}}"
      ))
    }
    cli::cli_alert_info("Group 1: {.val {levels(groups)[1]}}")
    cli::cli_alert_info("Group 2: {.val {levels(groups)[2]}}")
  } else if (method %in% c("anova", "kruskal")) {
    if (n_groups < 2) {
      cli::cli_abort(c(
        "{.field {group_col}} must be a factor with at least 2 levels for {.val {method}}",
        "i" = "Current levels: {.val {levels(groups)}}"
      ))
    }
    cli::cli_alert_info("Number of groups: {.val {n_groups}}")
    cli::cli_alert_info("Groups: {.val {levels(groups)}}")
  }

  # Check package availability based on method
  if (method == "kruskal") {
    .check_pkg_available("FSA")
  }

  # Perform DEA
  result <- switch(method,
    "t-test" = .gly_dea_2groups(expr_mat, groups, stats::t.test, p_adj_method, return_raw, ...),
    "wilcoxon" = .gly_dea_2groups(expr_mat, groups, stats::wilcox.test, p_adj_method, return_raw, ...),
    "anova" = .gly_dea_multi_groups(expr_mat, groups, stats::aov, stats::TukeyHSD, p_adj_method, return_raw, ...),
    "kruskal" = .gly_dea_multi_groups(expr_mat, groups, stats::kruskal.test, FSA::dunnTest, p_adj_method, return_raw, ...),
    cli::cli_abort("Invalid method: {.val {method}}")
  )

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Add S3 class
  subclass <- switch(method,
    "t-test" = "glystats_dea_res_ttest",
    "wilcoxon" = "glystats_dea_res_wilcoxon",
    "anova" = "glystats_dea_res_anova",
    "kruskal" = "glystats_dea_res_kruskal"
  )
  structure(result, class = c(subclass, "glystats_dea_res", "glystats_res", class(result)))
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
  # Convert main test results to tibble
  main_test_tbl <- .tibblify_main_test_results(mod_list$main_test, .f, p_adj_method)
  
  # Convert post-hoc results to tibble
  posthoc_tbl <- .tibblify_posthoc_results(mod_list$post_hoc, .f)
  
  list(main_test = main_test_tbl, post_hoc = posthoc_tbl)
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
.tibblify_posthoc_results <- function(posthoc_raw, .f) {
  if (length(posthoc_raw) == 0) {
    return(tibble::tibble())
  }
  
  # Convert each raw post-hoc result to tibble format
  posthoc_list <- list()
  
  for (var_name in names(posthoc_raw)) {
    raw_result <- posthoc_raw[[var_name]]
    
    if (identical(.f, stats::aov)) {
      # TukeyHSD result processing
      tukey_df <- as.data.frame(raw_result$group)
      tukey_df$comparison <- rownames(tukey_df)
      tukey_df$variable <- var_name
      posthoc_list[[var_name]] <- tukey_df
    } else {
      # Dunn test result processing
      dunn_df <- raw_result$res
      dunn_df$variable <- var_name
      posthoc_list[[var_name]] <- dunn_df
    }
  }
  
  # Combine all results
  result_tbl <- dplyr::bind_rows(posthoc_list) %>%
    janitor::clean_names()
  
  result_tbl
}
