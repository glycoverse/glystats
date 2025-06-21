#' Differential Expression Analysis (DEA)
#'
#' Perform differential expression analysis for glycomics or glycoproteomics data
#' using various statistical methods. The function supports both two-group comparisons
#' (t-test, Wilcoxon test) and multi-group comparisons (ANOVA, Kruskal-Wallis test).
#' For multi-group methods, post-hoc tests are automatically performed for significant results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param method A character string specifying the statistical method to use:
#'   - `"t-test"`: Student's t-test for two-group comparison (parametric)
#'   - `"wilcoxon"`: Wilcoxon rank-sum test for two-group comparison (non-parametric)
#'   - `"anova"`: One-way ANOVA for multi-group comparison (parametric)
#'   - `"kruskal"`: Kruskal-Wallis test for multi-group comparison (non-parametric)
#' @param group_col A character string specifying the column name of the grouping variable
#'   in the sample information. Default is `"group"`.
#' @param ... Additional arguments passed to the underlying statistical functions.
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
#' Post-hoc tests are only performed for variables with significant main effects (p < 0.05).
#'
#' @returns
#' For two-group methods (`"t-test"`, `"wilcoxon"`):
#' A tibble with statistical test results. For t-test, an additional `log2fc` column
#' (log2 fold change) is included.
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
gly_dea <- function(exp, method, group_col = "group", ...) {
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_choice(method, c("t-test", "wilcoxon", "anova", "kruskal"))
  checkmate::check_string(group_col)

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
  
  # Check group requirements based on method
  n_groups <- length(levels(groups))
  if (method %in% c("t-test", "wilcoxon")) {
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

  # Perform DEA
  result <- switch(method,
    "t-test" = .gly_dea_2groups(expr_mat, groups, stats::t.test, ...),
    "wilcoxon" = .gly_dea_2groups(expr_mat, groups, stats::wilcox.test, ...),
    "anova" = .gly_dea_multi_groups(expr_mat, groups, stats::aov, stats::TukeyHSD, ...),
    "kruskal" = .gly_dea_multi_groups(expr_mat, groups, stats::kruskal.test, FSA::dunnTest, ...),
    cli::cli_abort("Invalid method: {.val {method}}")
  )

  # Add S3 class
  subclass <- switch(method,
    "t-test" = "glystats_dea_res_ttest",
    "wilcoxon" = "glystats_dea_res_wilcoxon",
    "anova" = "glystats_dea_res_anova",
    "kruskal" = "glystats_dea_res_kruskal"
  )
  class(result) <- c(subclass, "glystats_dea_res", class(result))

  result
}

.gly_dea_2groups <- function(expr_mat, groups, .f, ...) {
  data <- expr_mat %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(group = groups) %>%
    tidyr::pivot_longer(cols = -all_of(c("sample", "group")), names_to = "variable", values_to = "value") %>%
    dplyr::mutate(log_value = log2(.data$value + 1))

  ttest_res <- data %>%
    dplyr::nest_by(.data$variable) %>%
    dplyr::mutate(
      ttest = list(.f(log_value ~ group, data = .data$data)),
      params = list(parameters::model_parameters(.data$ttest)),
      params = list(parameters::standardize_names(.data$params)),
      ) %>%
    dplyr::select(all_of(c("variable", "params"))) %>%
    tidyr::unnest(all_of("params")) %>%
    dplyr::ungroup() %>%
    janitor::clean_names()

  if (identical(.f, stats::t.test)) {
    ttest_res <- dplyr::mutate(ttest_res, log2fc = .data$mean_group1 - .data$mean_group2)
  }

  ttest_res
}


.gly_dea_multi_groups <- function(expr_mat, groups, .f, .ph, ...) {
  # .f is aov or kruskal.test
  # .ph is post-hoc test, i.e., TukeyHSD for aov or FSA::dunnTest for kruskal.test

  # Helper function to perform post-hoc tests
  .perform_posthoc_test <- function(data_nested, .f) {
    if (identical(.f, stats::aov)) {
      # TukeyHSD for ANOVA
      aov_model <- stats::aov(log_value ~ group, data = data_nested)
      tukey_result <- stats::TukeyHSD(aov_model)
      # Extract pairwise comparisons
      tukey_df <- as.data.frame(tukey_result$group)
      tukey_df$comparison <- rownames(tukey_df)
      tukey_df
    } else {
      # Dunn test for Kruskal-Wallis
      dunn_result <- FSA::dunnTest(log_value ~ group, data = data_nested, method = "holm")
      dunn_result$res
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

  # Helper function to run main tests on all variables
  .run_main_tests <- function(data, .f) {
    result <- data %>%
      dplyr::nest_by(.data$variable) %>%
      dplyr::mutate(
        test_result = list(.f(log_value ~ group, data = .data$data)),
        params = list(parameters::model_parameters(.data$test_result)),
        params = list(parameters::standardize_names(.data$params)),
      ) %>%
      dplyr::select(all_of(c("variable", "params"))) %>%
      tidyr::unnest(all_of("params")) %>%
      dplyr::ungroup() %>%
      janitor::clean_names()

    # For ANOVA, filter to keep only group effects (not residuals)
    # For Kruskal-Wallis, keep all results (there's only one row per variable anyway)
    if (identical(.f, stats::aov)) {
      result <- result %>%
        dplyr::filter(.data$parameter == "group")
    }

    result
  }

  # Helper function to generate post-hoc results
  .generate_posthoc_results <- function(main_test_res, data, .f) {
    # Filter significant results (p < 0.05)
    significant_vars <- main_test_res %>%
      dplyr::filter(.data$p < 0.05) %>%
      dplyr::pull(.data$variable)
    
    if (length(significant_vars) > 0) {
      # Perform post-hoc tests for significant variables
      posthoc_results <- data %>%
        dplyr::filter(.data$variable %in% significant_vars) %>%
        dplyr::nest_by(.data$variable) %>%
        dplyr::mutate(posthoc = list(.perform_posthoc_test(.data$data, .f))) %>%
        dplyr::select(all_of(c("variable", "posthoc"))) %>%
        tidyr::unnest(all_of("posthoc")) %>%
        dplyr::ungroup() %>%
        janitor::clean_names()
      
      return(posthoc_results)
    } else {
      # No significant results, return empty tibble with proper structure
      return(tibble::tibble())
    }
  }

  if (length(levels(groups)) < 2) {
    cli::cli_abort("Multi-group analysis requires at least 2 groups")
  }
  data <- .prepare_multi_group_data(expr_mat, groups)
  main_test_res <- .run_main_tests(data, .f)
  posthoc_res <- .generate_posthoc_results(main_test_res, data, .f)

  list(main_test = main_test_res, post_hoc = posthoc_res)
}
