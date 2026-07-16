#' One-way ANOVA for Differential Expression Analysis
#'
#' Perform one-way ANOVA for glycomics or glycoproteomics data.
#' The function supports parametric comparison of multiple groups.
#' For significant results, Tukey's HSD post-hoc test is automatically performed.
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
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#' @param ... Additional arguments passed to `stats::aov()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1e-6)) before
#' statistical testing. At least 2 groups are required in the grouping variable.
#'
#' For any variable failed to fit a `stats::aov()` model,
#' NAs will be assigned to the results in both main test and post-hoc test.
#'
#' **Post-hoc Test:**
#' Tukey's HSD test for pairwise comparisons (`stats::TukeyHSD()`) is performed
#' for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' A list containing three elements:
#'   - `tidy_result`: A list containing:
#'     - `main_test`: A tibble with ANOVA results containing the following columns:
#'       - `variable`: Variable name
#'       - `term`: ANOVA term (usually "groups")
#'       - `df`: Degrees of freedom
#'       - `sumsq`: Sum of squares
#'       - `meansq`: Mean squares
#'       - `statistic`: F-statistic
#'       - `p_val`: Raw p-value from ANOVA
#'       - `effect_size`: Eta-squared
#'       - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'       - `post_hoc`: Significant group pairs from post-hoc test, in the format of "ref_vs_test".
#'     - `post_hoc_test`: A tibble with pairwise comparison results containing the following columns:
#'       - `variable`: Variable name
#'       - `ref_group`: Reference group
#'       - `test_group`: Test/treatment/case group
#'       - `p_val`: Raw p-value from Tukey's HSD test
#'       - `p_adj`: Adjusted p-value from Tukey's HSD test
#'   - `raw_result`: A list containing:
#'     - `main_test`: A list of raw `aov` model objects.
#'     - `post_hoc_test`: A list of raw `TukeyHSD` objects. Post-hoc
#'       comparison labels follow the package direction, i.e.
#'       `test_group - ref_group`.
#'   - `meta_data`: A list containing metadata from the input experiment.
#'
#' @seealso [stats::aov()], [stats::TukeyHSD()]
#' @export
gly_anova <- function(
  exp,
  group_col = "group",
  p_adj_method = "BH",
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
    max_count = NULL,
    method = "anova"
  )
  groups <- group_info$groups

  # Run the internal computation
  result <- .analyze_anova(expr_mat, groups, p_adj_method, ...)
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)
  result
}

#' Run the internal computation for `gly_anova()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_anova <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
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
  if (length(levels(groups)) < 2) {
    cli::cli_abort("groups must have at least 2 levels for ANOVA")
  }

  # Prepare data
  log_expr_mat <- .log_transform_expr_mat(expr_mat)
  data <- .prepare_multi_group_data(log_expr_mat, groups, is_logged = TRUE)

  # Generate raw results
  raw_main_test <- .generate_raw_main_results(data, stats::aov, ...)
  raw_post_hoc_test <- .generate_raw_posthoc_results(
    raw_main_test,
    data,
    stats::aov,
    p_adj_method
  )

  # Tibblify results
  main_test <- .tibblify_main_test_results(
    raw_main_test,
    stats::aov,
    p_adj_method,
    effect_size_method = "eta_squared",
    expr_mat = log_expr_mat,
    groups = groups
  )
  post_hoc_vec <- .format_posthoc_results(
    raw_post_hoc_test,
    stats::aov,
    main_test$variable
  )
  main_test$post_hoc <- post_hoc_vec
  post_hoc_test <- .tibblify_posthoc_results(raw_post_hoc_test, stats::aov)
  post_hoc_test <- .add_fold_change(post_hoc_test, expr_mat, groups)

  # Assemble tidy results
  tidy_result <- list(
    main_test = main_test,
    post_hoc_test = post_hoc_test
  )

  # Assemble raw results
  raw_result <- list(
    main_test = raw_main_test,
    post_hoc_test = raw_post_hoc_test
  )

  # Assemble final result
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result
    ),
    class = c("glystats_anova_res", "glystats_res")
  )
}

#' Analysis of Covariance (ANCOVA) for Differential Expression Analysis
#'
#' Perform ANCOVA for glycomics or glycoproteomics data with sample-level covariates.
#' The function supports parametric comparison of multiple groups while adjusting for covariates.
#' For significant results, emmeans post-hoc comparisons (Tukey adjustment) are automatically performed.
#' P-values are adjusted for multiple testing using the method specified by `p_adj_method`.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing an expression matrix and
#'   sample information.
#' @param group_col A character string specifying the column name of the grouping variable
#'  in the sample information. Default is `"group"`.
#' @param covariate_cols A character vector specifying column names in sample information
#'  to include as covariates. At least one covariate must be provided.
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'  See `p.adjust.methods` for available methods. Default is "BH".
#'  If NULL, no adjustment is performed.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#' @param ... Additional arguments passed to `stats::aov()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1e-6)) before
#' statistical testing. At least 2 groups and at least 1 covariate are required.
#'
#' For any variable failed to fit a `stats::aov()` model,
#' NAs will be assigned to the results in both main test and post-hoc test.
#'
#' **Post-hoc Test:**
#' emmeans pairwise comparisons with Tukey adjustment (`emmeans::emmeans()`) are performed
#' for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' A list containing two elements:
#'   - `tidy_result`: A list containing:
#'     - `main_test`: A tibble with ANCOVA results containing the following columns:
#'       - `variable`: Variable name
#'       - `term`: ANCOVA term (usually "group")
#'       - `df`: Degrees of freedom
#'       - `sumsq`: Sum of squares
#'       - `meansq`: Mean squares
#'       - `statistic`: F-statistic
#'       - `p_val`: Raw p-value from ANCOVA
#'       - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'       - `post_hoc`: Significant group pairs from post-hoc test, in the format of "ref_vs_test".
#'     - `post_hoc_test`: A tibble with pairwise comparison results containing the following columns:
#'       - `variable`: Variable name
#'       - `ref_group`: Reference group
#'       - `test_group`: Test/treatment/case group
#'       - `p_val`: Adjusted p-value from emmeans pairwise test
#'       - `p_adj`: Adjusted p-value from emmeans pairwise test
#'   - `raw_result`: A list containing:
#'     - `main_test`: A list of raw `aov` model objects.
#'     - `post_hoc_test`: A list of raw emmeans results.
#'
#' @section Required packages:
#' This function requires the `emmeans` package for post-hoc comparisons.
#'
#' @seealso [stats::aov()], [emmeans::emmeans()], [emmeans::contrast()]
#' @export
gly_ancova <- function(
  exp,
  group_col = "group",
  covariate_cols = NULL,
  p_adj_method = "BH",
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
  rlang::check_installed("emmeans")
  if (is.null(covariate_cols) || length(covariate_cols) == 0) {
    cli::cli_abort("covariate_cols must be provided for ANCOVA.")
  }
  checkmate::assert_character(covariate_cols, min.len = 1)

  # Extract data from experiment object
  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "ancova"
  )
  groups <- group_info$groups

  # Extract covariates
  covariates <- .extract_covariates_from_sample_info(
    sample_info,
    covariate_cols,
    group_col
  )
  if (is.null(covariates) || ncol(covariates) == 0) {
    cli::cli_abort("covariate_cols must be provided for ANCOVA.")
  }

  # Run the internal computation
  result <- .analyze_ancova(expr_mat, groups, covariates, p_adj_method, ...)
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  result
}

#' Run the internal computation for `gly_ancova()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_ancova <- function(
  expr_mat,
  groups,
  covariates,
  p_adj_method = "BH",
  ...
) {
  # Check package availability
  rlang::check_installed("emmeans")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )

  # Validate group count
  if (length(levels(groups)) < 2) {
    cli::cli_abort("groups must have at least 2 levels for ANCOVA")
  }

  # Normalize and validate covariates
  covariates <- .normalize_covariates(
    covariates,
    ncol(expr_mat),
    colnames(expr_mat)
  )
  if (is.null(covariates) || ncol(covariates) == 0) {
    cli::cli_abort("covariates must be provided for ANCOVA.")
  }
  if ("group" %in% colnames(covariates)) {
    cli::cli_abort(
      "covariates cannot include a column named {.field group}; rename it."
    )
  }

  # Prepare data
  log_expr_mat <- .log_transform_expr_mat(expr_mat)
  data <- .prepare_multi_group_data(
    log_expr_mat,
    groups,
    covariates,
    is_logged = TRUE
  )
  ancova_formula <- stats::reformulate(
    termlabels = c("group", colnames(covariates)),
    response = "log_value"
  )

  # Generate raw results
  raw_main_test <- .generate_raw_main_results(
    data,
    stats::aov,
    formula = ancova_formula,
    ...
  )
  raw_post_hoc_test <- .generate_raw_posthoc_results(
    raw_main_test,
    data,
    stats::aov,
    p_adj_method,
    formula = ancova_formula,
    posthoc_method = "emmeans"
  )

  # Tibblify results
  main_test <- .tibblify_main_test_results(
    raw_main_test,
    stats::aov,
    p_adj_method
  )
  post_hoc_vec <- .format_posthoc_results(
    raw_post_hoc_test,
    stats::aov,
    main_test$variable
  )
  main_test$post_hoc <- post_hoc_vec
  post_hoc_test <- .tibblify_posthoc_results(raw_post_hoc_test, stats::aov)
  post_hoc_test <- .add_fold_change(post_hoc_test, expr_mat, groups)

  # Assemble tidy results
  tidy_result <- list(
    main_test = main_test,
    post_hoc_test = post_hoc_test
  )

  # Assemble raw results
  raw_result <- list(
    main_test = raw_main_test,
    post_hoc_test = raw_post_hoc_test
  )

  # Assemble final result
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result
    ),
    class = c("glystats_ancova_res", "glystats_res")
  )
}

#' Kruskal-Wallis test for Differential Expression Analysis
#'
#' Perform Kruskal-Wallis test for glycomics or glycoproteomics data.
#' The function supports non-parametric comparison of multiple groups.
#' For significant results, Dunn's post-hoc test is automatically performed.
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
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#' @param ... Additional arguments passed to `stats::kruskal.test()`.
#'
#' @section Required packages:
#' This function requires the `FSA` package for Dunn's post-hoc test.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1e-6)) before
#' statistical testing. At least 2 groups are required in the grouping variable.
#'
#' For any variable failed to fit a `stats::kruskal.test()` model,
#' NAs will be assigned to the results in both main test and post-hoc test.
#'
#' **Post-hoc Test:**
#' Dunn's test with Holm correction for multiple comparisons (`FSA::dunnTest()`) is performed
#' for variables with significant main effects (p_adj < 0.05).
#'
#' @returns
#' A list containing three elements:
#'   - `tidy_result`: A list containing:
#'     - `main_test`: A tibble with Kruskal-Wallis test results containing the following columns:
#'       - `variable`: Variable name
#'       - `statistic`: Kruskal-Wallis test statistic
#'       - `p_val`: Raw p-value from Kruskal-Wallis test
#'       - `parameter`: Degrees of freedom
#'       - `method`: Statistical method used
#'       - `effect_size`: Epsilon-squared
#'       - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'       - `post_hoc`: Significant group pairs from post-hoc test, in the format of "ref_vs_test".
#'     - `post_hoc_test`: A tibble with pairwise comparison results containing the following columns:
#'       - `variable`: Variable name
#'       - `ref_group`: Reference group
#'       - `test_group`: Test/treatment/case group
#'       - `p_val`: Raw p-value from Dunn's test
#'       - `p_adj`: Adjusted p-value from Dunn's test
#'   - `raw_result`: A list containing:
#'     - `main_test`: A list of raw `kruskal.test` objects.
#'     - `post_hoc_test`: A list of raw `dunnTest` objects. Post-hoc
#'       comparison labels and direction-sensitive statistics follow the package
#'       direction, i.e. `test_group - ref_group`.
#'   - `meta_data`: A list containing metadata from the input experiment.
#'
#' @seealso [stats::kruskal.test()], [FSA::dunnTest()]
#' @export
gly_kruskal <- function(
  exp,
  group_col = "group",
  p_adj_method = "BH",
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

  # Check package availability
  rlang::check_installed("FSA")

  # Extract data from experiment object
  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "kruskal"
  )
  groups <- group_info$groups

  # Run the internal computation
  result <- .analyze_kruskal(expr_mat, groups, p_adj_method, ...)
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)
  result
}

#' Run the internal computation for `gly_kruskal()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_kruskal <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  ...
) {
  # Check package availability
  rlang::check_installed("FSA")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )

  # Validate group count
  if (length(levels(groups)) < 2) {
    cli::cli_abort("groups must have at least 2 levels for Kruskal-Wallis test")
  }

  # Prepare data
  log_expr_mat <- .log_transform_expr_mat(expr_mat)
  data <- .prepare_multi_group_data(log_expr_mat, groups, is_logged = TRUE)

  # Generate raw results
  raw_main_test <- .generate_raw_main_results(data, stats::kruskal.test, ...)
  raw_post_hoc_test <- .generate_raw_posthoc_results(
    raw_main_test,
    data,
    stats::kruskal.test,
    p_adj_method
  )

  # Tibblify results
  main_test <- .tibblify_main_test_results(
    raw_main_test,
    stats::kruskal.test,
    p_adj_method,
    effect_size_method = "epsilon_squared",
    expr_mat = log_expr_mat,
    groups = groups
  )
  post_hoc_vec <- .format_posthoc_results(
    raw_post_hoc_test,
    stats::kruskal.test,
    main_test$variable
  )
  main_test$post_hoc <- post_hoc_vec
  post_hoc_test <- .tibblify_posthoc_results(
    raw_post_hoc_test,
    stats::kruskal.test
  )
  post_hoc_test <- .add_fold_change(post_hoc_test, expr_mat, groups)

  # Assemble tidy results
  tidy_result <- list(
    main_test = main_test,
    post_hoc_test = post_hoc_test
  )

  # Assemble raw results
  raw_result <- list(
    main_test = raw_main_test,
    post_hoc_test = raw_post_hoc_test
  )

  # Assemble final result
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result
    ),
    class = c("glystats_kruskal_res", "glystats_res")
  )
}

# Internal helper functions for multi-group analysis (ANOVA and Kruskal-Wallis)

# Helper function to perform raw post-hoc tests (returns raw objects)
.perform_raw_posthoc_test <- function(
  data_nested,
  .f,
  formula = NULL,
  posthoc_method = "tukey"
) {
  if (identical(.f, stats::aov)) {
    # Post-hoc for ANOVA/ANCOVA
    if (is.null(formula)) {
      formula <- log_value ~ group
    }
    aov_model <- stats::aov(formula, data = data_nested)
    if (identical(posthoc_method, "emmeans")) {
      emm <- emmeans::emmeans(aov_model, specs = "group")
      contrast <- emmeans::contrast(
        emm,
        method = "pairwise",
        adjust = "tukey",
        reverse = FALSE
      )
      as.data.frame(summary(contrast))
    } else {
      suppressWarnings(stats::TukeyHSD(aov_model, which = "group"))
    }
  } else {
    # Dunn test for Kruskal-Wallis
    # Check number of groups - FSA::dunnTest fails with only 2 groups
    n_groups <- length(unique(data_nested$group))

    if (n_groups == 2) {
      # Use pairwise Wilcoxon test for 2 groups
      wilcox_result <- stats::pairwise.wilcox.test(
        data_nested$log_value,
        data_nested$group,
        p.adjust.method = "holm"
      )

      # Convert to dunnTest-like format for consistency
      # In pairwise.wilcox.test, column name is the first group (ref_group),
      # row name is the second group (test_group)
      p_matrix <- wilcox_result$p.value
      comparisons <- c()
      p_values <- c()

      for (i in seq_len(nrow(p_matrix))) {
        for (j in seq_len(ncol(p_matrix))) {
          if (!is.na(p_matrix[i, j])) {
            group1 <- colnames(p_matrix)[j] # ref_group (column name)
            group2 <- rownames(p_matrix)[i] # test_group (row name)
            comparisons <- c(comparisons, paste(group2, "-", group1))
            p_values <- c(p_values, p_matrix[i, j])
          }
        }
      }

      # Create dunnTest-like structure
      list(
        res = data.frame(
          Comparison = comparisons,
          P.unadj = p_values,
          P.adj = p_values,
          stringsAsFactors = FALSE
        ),
        method = wilcox_result$p.adjust.method
      )
    } else {
      # Use FSA::dunnTest for 3+ groups
      .standardize_dunn_posthoc_result(
        FSA::dunnTest(log_value ~ group, data = data_nested, method = "holm")
      )
    }
  }
}

# Helper function to generate raw main test results
.generate_raw_main_results <- function(data, .f, formula = NULL, ...) {
  dots <- rlang::list2(...)
  disallowed_args <- intersect(names(dots), c("formula", "data"))
  if (length(disallowed_args) > 0) {
    cli::cli_abort(
      "Arguments {cli::format_inline(disallowed_args)} should not be supplied through `...`; they are controlled internally."
    )
  }
  if (is.null(formula)) {
    formula <- log_value ~ group
  }
  safe_f <- purrr::possibly(function(...) rlang::exec(.f, ...), otherwise = NA)
  main_test_raw <- data %>%
    dplyr::nest_by(.data$variable) %>%
    dplyr::mutate(
      test_result = list(safe_f(formula, data = .data$data, !!!dots))
    )
  main_test_list <- main_test_raw$test_result
  n_na <- sum(is.na(main_test_list))
  if (n_na > 0) {
    cli::cli_warn("{.val {n_na}} variables failed to fit the model")
  }
  names(main_test_list) <- main_test_raw$variable
  main_test_list
}

.generate_raw_posthoc_results <- function(
  main_test_raw,
  data,
  .f,
  p_adj_method,
  formula = NULL,
  posthoc_method = "tukey"
) {
  p_fn <- ifelse(
    identical(.f, stats::aov),
    function(x) summary(x)[[1]][["Pr(>F)"]][1],
    function(x) x$p.value
  )
  valid_main_test_mods <- main_test_raw[!is.na(main_test_raw)]
  main_test_p_vals <- purrr::map_dbl(valid_main_test_mods, p_fn)
  if (!is.null(p_adj_method)) {
    main_test_p_vals <- stats::p.adjust(main_test_p_vals, method = p_adj_method)
  }
  sig_vars <- names(main_test_p_vals)[main_test_p_vals < 0.05]
  if (length(sig_vars) > 0) {
    posthoc_raw_results <- data %>%
      dplyr::filter(.data$variable %in% sig_vars) %>%
      dplyr::nest_by(.data$variable) %>%
      dplyr::mutate(
        posthoc_raw = list(.perform_raw_posthoc_test(
          .data$data,
          .f,
          formula = formula,
          posthoc_method = posthoc_method
        ))
      ) %>%
      dplyr::select(all_of(c("variable", "posthoc_raw")))
    raw_posthoc_list <- posthoc_raw_results$posthoc_raw
    names(raw_posthoc_list) <- posthoc_raw_results$variable
    return(raw_posthoc_list)
  } else {
    return(list())
  }
}

# Helper function to prepare data for analysis
.prepare_multi_group_data <- function(
  expr_mat,
  groups,
  covariates = NULL,
  is_logged = FALSE
) {
  data <- expr_mat %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(group = groups)

  if (!is.null(covariates)) {
    covariate_tbl <- tibble::as_tibble(covariates)
    data <- dplyr::bind_cols(data, covariate_tbl)
  }

  keep_cols <- c("sample", "group")
  if (!is.null(covariates)) {
    keep_cols <- c(keep_cols, colnames(covariates))
  }

  data <- data %>%
    tidyr::pivot_longer(
      cols = -all_of(keep_cols),
      names_to = "variable",
      values_to = "log_value"
    )

  if (!is_logged) {
    data <- data %>%
      dplyr::mutate(log_value = log2(.data$log_value + 1e-6))
  }

  data
}

# Helper function to convert main test raw results to tibble
.tibblify_main_test_results <- function(
  main_test_raw,
  .f,
  p_adj_method,
  effect_size_method = NULL,
  expr_mat = NULL,
  groups = NULL
) {
  var_names <- names(main_test_raw)
  result_tbl <- tibble::tibble(
    variable = var_names,
    test_result = main_test_raw
  ) %>%
    dplyr::mutate(
      params = purrr::map(
        .data$test_result,
        ~ {
          if (rlang::is_na(.x)) NULL else broom::tidy(.x)
        }
      )
    ) %>%
    dplyr::select(all_of(c("variable", "params"))) %>%
    tidyr::unnest(all_of("params"), keep_empty = TRUE) %>%
    dplyr::ungroup() %>%
    janitor::clean_names()

  # For ANOVA, filter to keep only group effects (not residuals)
  # For Kruskal-Wallis, keep all results (there's only one row per variable anyway)
  if (identical(.f, stats::aov)) {
    result_tbl <- result_tbl %>%
      dplyr::filter(.data$term == "group")
  }

  # Rename p_value to p_val for consistency
  if ("p_value" %in% colnames(result_tbl)) {
    result_tbl <- dplyr::rename(result_tbl, p_val = "p_value")
  }

  if (!is.null(p_adj_method)) {
    result_tbl <- dplyr::mutate(
      result_tbl,
      p_adj = stats::p.adjust(.data$p_val, method = p_adj_method)
    )
  }

  if (!is.null(effect_size_method)) {
    if (is.null(expr_mat) || is.null(groups)) {
      cli::cli_abort(
        "{.arg expr_mat} and {.arg groups} must be supplied when {.arg effect_size_method} is used."
      )
    }

    result_tbl <- .add_effect_size_to_main_test(
      result_tbl,
      effect_size_method,
      expr_mat,
      groups
    )
  }

  result_tbl
}

#' Add effect sizes to multi-group main test results
#'
#' @param result_tbl A tibble containing the tidy main test results.
#' @param effect_size_method A character string naming the effect size to compute.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor specifying group membership for each sample.
#'
#' @return A tibble with an `effect_size` column added.
#' @noRd
.add_effect_size_to_main_test <- function(
  result_tbl,
  effect_size_method,
  expr_mat,
  groups
) {
  checkmate::assert_choice(
    effect_size_method,
    c("eta_squared", "epsilon_squared")
  )

  effect_size_values <- purrr::map2_dbl(
    result_tbl$variable,
    result_tbl$statistic,
    function(var_name, statistic) {
      if (identical(effect_size_method, "eta_squared")) {
        .calculate_eta_squared(expr_mat, groups, var_name)
      } else {
        .calculate_epsilon_squared(expr_mat, groups, var_name, statistic)
      }
    }
  )

  result_tbl$effect_size <- effect_size_values
  result_tbl
}

#' Calculate eta-squared for one-way ANOVA
#'
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor specifying group membership for each sample.
#' @param var_name A character string specifying the variable to evaluate.
#'
#' @return A numeric scalar containing eta-squared.
#' @noRd
.calculate_eta_squared <- function(expr_mat, groups, var_name) {
  log_values <- expr_mat[var_name, ]
  valid_idx <- !is.na(log_values) & !is.na(groups)
  log_values <- log_values[valid_idx]
  groups <- droplevels(groups[valid_idx])

  if (length(log_values) == 0 || nlevels(groups) < 2) {
    return(NA_real_)
  }

  grand_mean <- mean(log_values)
  group_means <- purrr::map_dbl(levels(groups), function(group) {
    mean(log_values[groups == group])
  })
  group_sizes <- purrr::map_int(levels(groups), function(group) {
    sum(groups == group)
  })
  ss_between <- sum(group_sizes * (group_means - grand_mean)^2)
  ss_total <- sum((log_values - grand_mean)^2)

  if (!is.finite(ss_total) || ss_total == 0) {
    return(NA_real_)
  }

  ss_between / ss_total
}

#' Calculate epsilon-squared for the Kruskal-Wallis test
#'
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor specifying group membership for each sample.
#' @param var_name A character string specifying the variable to evaluate.
#'
#' @return A numeric scalar containing epsilon-squared.
#' @noRd
.calculate_epsilon_squared <- function(expr_mat, groups, var_name, h_stat) {
  log_values <- expr_mat[var_name, ]
  valid_idx <- !is.na(log_values) & !is.na(groups)
  log_values <- log_values[valid_idx]
  groups <- droplevels(groups[valid_idx])

  if (
    length(log_values) == 0 ||
      nlevels(groups) < 2 ||
      length(log_values) <= nlevels(groups)
  ) {
    return(NA_real_)
  }

  if (!is.finite(h_stat)) {
    return(NA_real_)
  }

  (h_stat - nlevels(groups) + 1) / (length(log_values) - nlevels(groups))
}

# Helper function to convert post-hoc raw results to tibble
.format_posthoc_results <- function(posthoc_raw, .f, variables) {
  if (length(posthoc_raw) == 0) {
    return(rep(NA_character_, length(variables)))
  }
  posthoc_map <- purrr::imap(posthoc_raw, function(raw_result, var_name) {
    if (identical(.f, stats::aov)) {
      if (is.data.frame(raw_result) && "contrast" %in% colnames(raw_result)) {
        sig_pairs <- raw_result$contrast[raw_result$p.value < 0.05]
        sig_pairs <- purrr::map_chr(sig_pairs, function(x) {
          parts <- stringr::str_split(x, "\\s*-\\s*", simplify = TRUE)
          stringr::str_c(parts[, 1], "_vs_", parts[, 2])
        })
      } else {
        tukey_df <- as.data.frame(raw_result$group)
        sig_pairs <- rownames(tukey_df)[tukey_df$"p adj" < 0.05]
        sig_pairs <- purrr::map(sig_pairs, function(x) {
          # switch two groups and change the "-" to "_vs_"
          parts <- stringr::str_split(x, "-", simplify = TRUE)
          stringr::str_c(parts[, 2], "_vs_", parts[, 1])
        })
      }
      sig_str <- if (length(sig_pairs) == 0) {
        NA_character_
      } else {
        paste(sig_pairs, collapse = ";")
      }
    } else {
      dunn_df <- raw_result$res
      sig_pairs <- dunn_df$Comparison[dunn_df$P.adj < 0.05]
      # Raw Dunn results are standardized to "test - ref". Convert them back to
      # the package-facing "ref_vs_test" label used elsewhere.
      sig_pairs <- purrr::map_chr(sig_pairs, function(x) {
        parts <- stringr::str_split(x, " - ", simplify = TRUE)
        stringr::str_c(parts[, 2], "_vs_", parts[, 1])
      })
      sig_str <- if (length(sig_pairs) == 0) {
        NA_character_
      } else {
        paste(sig_pairs, collapse = ";")
      }
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
      ref_group = character(0),
      test_group = character(0),
      p_val = numeric(0),
      p_adj = numeric(0)
    ))
  }

  # Convert each post-hoc result to tibble and combine
  result_list <- purrr::imap(posthoc_raw, function(raw_result, var_name) {
    if (identical(.f, stats::aov)) {
      if (is.data.frame(raw_result) && "contrast" %in% colnames(raw_result)) {
        comparison_parts <- stringr::str_split(
          raw_result$contrast,
          "\\s*-\\s*",
          simplify = TRUE
        )
        tibble::tibble(
          variable = var_name,
          ref_group = comparison_parts[, 1],
          test_group = comparison_parts[, 2],
          p_val = raw_result$p.value,
          p_adj = raw_result$p.value
        )
      } else {
        # For TukeyHSD results
        tukey_df <- as.data.frame(raw_result$group)
        tukey_df$comparison <- rownames(tukey_df)

        # Parse group comparisons (format: "group2-group1")
        comparison_parts <- stringr::str_split(
          tukey_df$comparison,
          "-",
          simplify = TRUE
        )

        tibble::tibble(
          variable = var_name,
          ref_group = comparison_parts[, 2], # Second part is ref_group
          test_group = comparison_parts[, 1], # First part is test_group
          p_val = tukey_df$`p adj`, # TukeyHSD already provides adjusted p-values
          p_adj = tukey_df$`p adj`
        )
      }
    } else {
      # For Dunn test results
      dunn_df <- raw_result$res

      # Parse group comparisons (format: "test_group - ref_group")
      comparison_parts <- stringr::str_split(
        dunn_df$Comparison,
        " - ",
        simplify = TRUE
      )

      tibble::tibble(
        variable = var_name,
        ref_group = comparison_parts[, 2],
        test_group = comparison_parts[, 1],
        p_val = dunn_df$P.unadj, # Unadjusted p-value
        p_adj = dunn_df$P.adj # Adjusted p-value
      )
    }
  })

  # Combine all results
  dplyr::bind_rows(result_list)
}

.add_fold_change <- function(post_hoc_test, expr_mat, groups) {
  if (length(levels(groups)) == 2) {
    fc_res <- .fc_2groups(expr_mat, groups)
    return(dplyr::left_join(post_hoc_test, fc_res, by = "variable"))
  } else {
    fc_res <- .fc_multi_groups(expr_mat, groups)
    return(dplyr::left_join(
      post_hoc_test,
      fc_res,
      by = c("ref_group", "test_group", "variable")
    ))
  }
}

.standardize_dunn_posthoc_result <- function(raw_result) {
  if (!("res" %in% names(raw_result)) || nrow(raw_result$res) == 0) {
    return(raw_result)
  }

  comparison_parts <- stringr::str_split(
    raw_result$res$Comparison,
    " - ",
    simplify = TRUE
  )
  raw_result$res$Comparison <- stringr::str_c(
    comparison_parts[, 2],
    " - ",
    comparison_parts[, 1]
  )

  if ("Z" %in% colnames(raw_result$res)) {
    raw_result$res$Z <- -raw_result$res$Z
  }

  raw_result
}
