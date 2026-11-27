#' One-way or Repeated-Measures ANOVA for Differential Expression Analysis
#'
#' Perform one-way or repeated-measures ANOVA for glycomics or
#' glycoproteomics data.
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
#' @param contrasts Optional character vector selecting post-hoc comparisons to
#'   report, e.g. `"H_vs_C"` compares C to reference H. The form `"H-C"` is
#'   also accepted when group names contain no hyphens. The omnibus test uses
#'   all groups, and post-hoc p-values retain the all-pairs adjustment. Raw
#'   post-hoc results contain all comparisons. By default all pairs are reported.
#'   A named list is passed to `stats::aov()` for factor contrast coding.
#' @param ... Additional arguments passed to `stats::aov()`.
#' @param subset Optional logical or integer index selecting samples in input
#'   order. Selection precedes subject matching, model fitting, post-hoc tests,
#'   and effect-size and fold-change calculations.
#' @param subject_col An optional character string naming the subject identifier
#'   column in sample information. When supplied, a subject-blocked
#'   repeated-measures ANOVA is performed.
#'
#' @details
#' For raw data with automatic transformation, the function performs log2 transformation on the expression data (log2(x + 1e-6)) before
#' statistical testing. At least 2 groups are required in the grouping variable.
#' In paired analyses, only subjects observed in every group are used, with at
#' most one sample per subject and group. The reported effect size is partial
#' eta-squared for the group term.
#'
#' For any variable failed to fit a `stats::aov()` model,
#' NAs will be assigned to the results in both main test and post-hoc test.
#'
#' **Post-hoc Test:**
#' Tukey's HSD test for pairwise comparisons (`stats::TukeyHSD()`) is performed
#' for variables with significant main effects (p_adj < 0.05). Paired analyses
#' instead use paired t-tests with Holm correction.
#'
#' @returns
#' A list containing four elements:
#'   - `tidy_result`: A list containing:
#'     - `main_test`: A tibble with ANOVA results containing the following columns:
#'       - `variable`: Variable name
#'       - `term`: ANOVA term (usually "groups")
#'       - `df`: Degrees of freedom
#'       - `sumsq`: Sum of squares
#'       - `meansq`: Mean squares
#'       - `statistic`: F-statistic
#'       - `p_val`: Raw p-value from ANOVA
#'       - `effect_size`: Eta-squared, or partial eta-squared for paired analyses
#'       - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'       - `post_hoc`: Significant group pairs from post-hoc test, in the format of "ref_vs_test".
#'     - `post_hoc_test`: A tibble with pairwise comparison results containing the following columns:
#'       - `variable`: Variable name
#'       - `ref_group`: Reference group
#'       - `test_group`: Test/treatment/case group
#'       - `p_val`: Raw post-hoc p-value
#'       - `p_adj`: Adjusted post-hoc p-value
#'   - `raw_result`: A list containing:
#'     - `main_test`: A list of raw `aov` model objects.
#'     - `post_hoc_test`: A list of raw `TukeyHSD` objects or paired t-test
#'       summaries.
#'   - `meta_data`: A list containing metadata from the input experiment.
#'   - `preprocessing`: Input and analysis scales and the applied transformation; see [gly_scale].
#'
#' @seealso [stats::aov()], [stats::TukeyHSD()]
#' @inheritParams gly_scale
#' @section SE identifiers:
#' For SE inputs, the result also contains an `identifiers` mapping of result
#' IDs to input positions and original names; see [gly_identity].
#' In SE tidy results, identity columns are `.variable`/`.sample` (dotted
#' endpoints for correlations); undotted metadata columns remain annotations.
#' @seealso [gly_scale] for input scales, transformations, and fold changes.
#'   [gly_identity] for SE identifiers, annotation alignment, and unnamed inputs.
#' @export
gly_anova <- function(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  add_info = TRUE,
  ...,
  subject_col = NULL,
  log_transform = NULL,
  contrasts = NULL,
  subset = NULL
) {
  # Validate inputs
  checkmate::assert_logical(
    log_transform,
    len = 1,
    any.missing = FALSE,
    null.ok = TRUE
  )
  .assert_data_container(exp)
  meta_data <- .get_meta_data(exp)
  transformation <- .resolve_transformation(
    exp,
    log_transform,
    default = TRUE,
    method = "log2",
    meta_data = meta_data
  )
  log_transform <- transformation$log_transform
  checkmate::assert_string(group_col)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_string(subject_col, null.ok = TRUE)

  # Extract data from experiment object
  input <- .multigroup_input(exp, subset, rlang::list2(...))
  expr_mat <- input$expr_mat
  sample_info <- input$sample_info

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "anova",
    show_info = !is.character(contrasts)
  )
  groups <- group_info$groups
  subjects <- .extract_paired_subjects(sample_info, subject_col, group_col)

  # Run the internal computation
  result <- .analyze_anova(
    expr_mat,
    groups,
    p_adj_method,
    !!!input$dots,
    subjects = subjects,
    contrasts = contrasts,
    log_transform = log_transform,
    input_scale = transformation$input_scale
  )
  if (!input$keep_qr) {
    result$raw_result$main_test <- .drop_aov_qr(result$raw_result$main_test)
  }
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  # Add meta_data from experiment
  result$meta_data <- meta_data
  result$preprocessing <- transformation
  .record_se_identity(result, exp)
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
  ...,
  subjects = NULL,
  log_transform = TRUE,
  input_scale = "raw",
  contrasts = NULL
) {
  checkmate::assert_flag(log_transform)
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

  model_args <- if (is.list(contrasts)) list(contrasts = contrasts)
  comparisons <- if (is.null(model_args)) {
    .parse_posthoc_comparisons(contrasts, groups)
  }

  design <- .match_paired_samples(
    expr_mat,
    groups,
    subjects,
    method = "ANOVA"
  )
  fit_args <- rlang::list2(..., !!!model_args)
  if (design$paired) {
    fit_args <- .index_sample_args(fit_args, design$indices, ncol(expr_mat))
  }
  expr_mat <- design$expr_mat
  groups <- design$groups

  # Prepare data
  log_expr_mat <- .log_transform_expr_mat(expr_mat, log_transform)
  data <- .prepare_multi_group_data(
    log_expr_mat,
    groups,
    subjects = design$subjects,
    is_logged = TRUE
  )
  anova_formula <- if (design$paired) {
    log_value ~ subject + group
  } else {
    NULL
  }

  # Generate raw results
  raw_main_test <- .generate_raw_main_results(
    data,
    stats::aov,
    formula = anova_formula,
    !!!fit_args,
    paired = design$paired
  )
  main_test <- .tibblify_main_test_results(
    raw_main_test,
    stats::aov,
    p_adj_method,
    effect_size_method = if (design$paired) {
      "partial_eta_squared"
    } else {
      "eta_squared"
    },
    expr_mat = log_expr_mat,
    groups = groups,
    effect_models = raw_main_test
  )
  raw_post_hoc_test <- .generate_raw_posthoc_results(
    raw_main_test,
    data,
    stats::aov,
    main_test,
    posthoc_method = if (design$paired) "paired_t" else "tukey",
    paired = design$paired
  )
  post_hoc_test <- .tibblify_posthoc_results(raw_post_hoc_test, stats::aov)
  post_hoc_test <- .add_fold_change(
    post_hoc_test,
    expr_mat,
    groups,
    paired_data = if (design$paired) data else NULL,
    input_scale = input_scale
  )

  tidy_result <- .report_posthoc(
    main_test,
    post_hoc_test,
    comparisons,
    levels(groups)
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

.multigroup_input <- function(exp, subset, dots, rank_test = FALSE) {
  backend_args <- if (!rank_test) {
    unique(c(names(formals(stats::aov)), names(formals(stats::lm)), "tol"))
  } else {
    c("subset", "na.action", "drop.unused.levels", "xlev")
  }
  dots <- .match_backend_args(
    dots,
    setdiff(backend_args, c("formula", "data", "..."))
  )
  if (!is.null(dots$subset)) {
    if (!is.null(subset)) {
      cli::cli_abort("{.arg subset} was supplied twice.")
    }
    subset <- dots$subset
    dots$subset <- NULL
  }
  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)
  if (!is.null(subset)) {
    if (
      !is.logical(subset) &&
        !checkmate::test_integerish(subset, any.missing = FALSE)
    ) {
      cli::cli_abort(
        "{.arg subset} must index samples with logical or integer values."
      )
    }
    indices <- seq_len(ncol(expr_mat))[subset]
    if (anyNA(indices)) {
      cli::cli_abort(
        "{.arg subset} contains missing or out-of-range sample indices."
      )
    }
    dots <- .index_sample_args(dots, indices, ncol(expr_mat))
    expr_mat <- expr_mat[, indices, drop = FALSE]
    sample_info <- sample_info[indices, , drop = FALSE]
  }
  keep_qr <- !identical(dots$qr, FALSE)
  if (!keep_qr) {
    dots$qr <- TRUE
  }
  list(
    expr_mat = expr_mat,
    sample_info = sample_info,
    dots = dots,
    keep_qr = keep_qr
  )
}

# Summary and Tukey calculations require QR even when the caller omits it
# from the returned fits.
.drop_aov_qr <- function(fits) {
  lapply(fits, function(fit) {
    if (inherits(fit, "aov")) {
      fit$qr <- NULL
      fit$call$qr <- FALSE
    }
    fit
  })
}

.index_sample_args <- function(dots, indices, n_samples) {
  for (name in intersect(names(dots), c("weights", "offset"))) {
    values <- dots[[name]]
    if (!is.null(values)) {
      if (length(values) != n_samples) {
        cli::cli_abort("{.arg {name}} must have one value per input sample.")
      }
      dots[[name]] <- values[indices]
    }
  }
  dots
}

.parse_posthoc_comparisons <- function(contrasts, groups) {
  if (is.null(contrasts)) {
    return(NULL)
  }
  comparisons <- .parse_custom_contrasts(contrasts, levels(groups))
  cli::cli_alert_info("Reporting post-hoc comparisons: {.val {contrasts}}")
  unique(unname(comparisons))
}

.report_posthoc <- function(
  main_test,
  post_hoc_test,
  comparisons,
  group_levels
) {
  if (!is.null(comparisons)) {
    post_hoc_test <- .select_posthoc_comparisons(
      post_hoc_test,
      comparisons,
      group_levels
    )
  }
  significant <- post_hoc_test |>
    dplyr::filter(.data$p_adj < 0.05) |>
    dplyr::group_by(.data$variable) |>
    dplyr::summarise(
      post_hoc = paste(
        paste0(.data$ref_group, "_vs_", .data$test_group),
        collapse = ";"
      )
    )
  main_test$post_hoc <- significant$post_hoc[match(
    main_test$variable,
    significant$variable
  )]

  list(
    main_test = main_test,
    post_hoc_test = post_hoc_test
  )
}

# The shared contrast parser returns each pair in test, reference order.
.select_posthoc_comparisons <- function(post_hoc, comparisons, group_levels) {
  purrr::map(comparisons, function(pair) {
    test <- pair[1]
    ref <- pair[2]
    reverse <- match(ref, group_levels) > match(test, group_levels)
    ordered <- if (reverse) pair else rev(pair)
    rows <- post_hoc$ref_group == ordered[1] & post_hoc$test_group == ordered[2]
    result <- post_hoc[rows, , drop = FALSE]
    if (reverse) {
      result <- dplyr::mutate(
        result,
        ref_group = ref,
        test_group = test,
        log2fc = -.data$log2fc
      )
    }
    result
  }) |>
    dplyr::bind_rows()
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
#' @inheritParams gly_anova
#' @param ... Additional arguments passed to `stats::aov()`.
#'
#' @details
#' For raw data with automatic transformation, the function performs log2 transformation on the expression data (log2(x + 1e-6)) before
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
#' A list containing four elements:
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
#'   - `meta_data`: Metadata copied from the input experiment.
#'   - `preprocessing`: Input and analysis scales and the applied transformation; see [gly_scale].
#'
#' @section Required packages:
#' This function requires the `emmeans` package for post-hoc comparisons.
#'
#' @seealso [stats::aov()], [emmeans::emmeans()], [emmeans::contrast()]
#' @inheritParams gly_scale
#' @section SE identifiers:
#' For SE inputs, the result also contains an `identifiers` mapping of result
#' IDs to input positions and original names; see [gly_identity].
#' In SE tidy results, identity columns are `.variable`/`.sample` (dotted
#' endpoints for correlations); undotted metadata columns remain annotations.
#' @seealso [gly_scale] for input scales, transformations, and fold changes.
#'   [gly_identity] for SE identifiers, annotation alignment, and unnamed inputs.
#' @export
gly_ancova <- function(
  exp,
  group_col = "group",
  covariate_cols = NULL,
  p_adj_method = "BH",
  add_info = TRUE,
  ...,
  log_transform = NULL,
  contrasts = NULL,
  subset = NULL
) {
  # Validate inputs
  checkmate::assert_logical(
    log_transform,
    len = 1,
    any.missing = FALSE,
    null.ok = TRUE
  )
  .assert_data_container(exp)
  meta_data <- .get_meta_data(exp)
  transformation <- .resolve_transformation(
    exp,
    log_transform,
    default = TRUE,
    method = "log2",
    meta_data = meta_data
  )
  log_transform <- transformation$log_transform
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
  input <- .multigroup_input(exp, subset, rlang::list2(...))
  expr_mat <- input$expr_mat
  sample_info <- input$sample_info

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "ancova",
    show_info = !is.character(contrasts)
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
  result <- .analyze_ancova(
    expr_mat,
    groups,
    covariates,
    p_adj_method,
    !!!input$dots,
    contrasts = contrasts,
    log_transform = log_transform,
    input_scale = transformation$input_scale
  )
  if (!input$keep_qr) {
    result$raw_result$main_test <- .drop_aov_qr(result$raw_result$main_test)
  }
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  result$meta_data <- meta_data
  result$preprocessing <- transformation
  .record_se_identity(result, exp)
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
  ...,
  log_transform = TRUE,
  input_scale = "raw",
  contrasts = NULL
) {
  checkmate::assert_flag(log_transform)
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

  model_args <- if (is.list(contrasts)) list(contrasts = contrasts)
  comparisons <- if (is.null(model_args)) {
    .parse_posthoc_comparisons(contrasts, groups)
  }

  # Normalize and validate covariates
  covariates <- .normalize_covariates(
    covariates,
    ncol(expr_mat),
    colnames(expr_mat),
    allow_subject = TRUE
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
  log_expr_mat <- .log_transform_expr_mat(expr_mat, log_transform)
  data <- .prepare_multi_group_data(
    log_expr_mat,
    groups,
    covariates,
    is_logged = TRUE
  )
  # aov uses sequential sums of squares. Enter the group after every covariate
  # to test its additional contribution to the fitted model.
  ancova_formula <- stats::reformulate(
    termlabels = c(colnames(covariates), "group"),
    response = "log_value"
  )

  # Generate raw results
  raw_main_test <- .generate_raw_main_results(
    data,
    stats::aov,
    formula = ancova_formula,
    ...,
    !!!model_args
  )
  main_test <- .tibblify_main_test_results(
    raw_main_test,
    stats::aov,
    p_adj_method
  )
  raw_post_hoc_test <- .generate_raw_posthoc_results(
    raw_main_test,
    data,
    stats::aov,
    main_test,
    posthoc_method = "emmeans"
  )
  post_hoc_test <- .tibblify_posthoc_results(raw_post_hoc_test, stats::aov)
  post_hoc_test <- .add_fold_change(
    post_hoc_test,
    expr_mat,
    groups,
    input_scale = input_scale
  )

  tidy_result <- .report_posthoc(
    main_test,
    post_hoc_test,
    comparisons,
    levels(groups)
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

#' Kruskal-Wallis or Friedman Test for Differential Expression Analysis
#'
#' Perform Kruskal-Wallis or Friedman tests for glycomics or glycoproteomics
#' data.
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
#' @inheritParams gly_anova
#' @param contrasts Optional character vector selecting post-hoc comparisons,
#'   e.g. `"H_vs_C"` with H as reference. As in [gly_anova()], the main test
#'   and the all-pairs adjustment use all groups. Raw results retain all pairs.
#' @param ... Additional arguments passed to `stats::kruskal.test()` or, for
#'   paired analyses, `stats::friedman.test()`.
#' @param subject_col An optional character string naming the subject identifier
#'   column in sample information. When supplied, a Friedman test is performed.
#'
#' @section Required packages:
#' Unpaired analyses require the `FSA` package for Dunn's post-hoc test.
#'
#' @details
#' For raw data with automatic transformation, the function performs log2 transformation on the expression data (log2(x + 1e-6)) before
#' statistical testing. At least 2 groups are required in the grouping variable.
#' In paired analyses, only subjects observed in every group are used, with at
#' most one sample per subject and group. The reported effect size is Kendall's
#' \eqn{W}.
#'
#' For any variable failed to fit a `stats::kruskal.test()` or
#' `stats::friedman.test()` model,
#' NAs will be assigned to the results in both main test and post-hoc test.
#'
#' **Post-hoc Test:**
#' Dunn's test with Holm correction for multiple comparisons (`FSA::dunnTest()`) is performed
#' for variables with significant main effects (p_adj < 0.05). Paired analyses
#' instead use paired Wilcoxon signed-rank tests with Holm correction.
#'
#' @returns
#' A list containing four elements:
#'   - `tidy_result`: A list containing:
#'     - `main_test`: A tibble with Kruskal-Wallis or Friedman test results
#'       containing the following columns:
#'       - `variable`: Variable name
#'       - `statistic`: Kruskal-Wallis test statistic
#'       - `p_val`: Raw p-value from Kruskal-Wallis test
#'       - `parameter`: Degrees of freedom
#'       - `method`: Statistical method used
#'       - `effect_size`: Epsilon-squared, or Kendall's W for paired analyses
#'       - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'       - `post_hoc`: Significant group pairs from post-hoc test, in the format of "ref_vs_test".
#'     - `post_hoc_test`: A tibble with pairwise comparison results containing the following columns:
#'       - `variable`: Variable name
#'       - `ref_group`: Reference group
#'       - `test_group`: Test/treatment/case group
#'       - `p_val`: Raw post-hoc p-value
#'       - `p_adj`: Holm-adjusted post-hoc p-value
#'   - `raw_result`: A list containing:
#'     - `main_test`: A list of raw `kruskal.test` or `friedman.test` objects.
#'     - `post_hoc_test`: A list of raw Dunn or paired Wilcoxon summaries.
#'   - `meta_data`: A list containing metadata from the input experiment.
#'   - `preprocessing`: Input and analysis scales and the applied transformation; see [gly_scale].
#'
#' @seealso [stats::kruskal.test()], [stats::friedman.test()],
#'   [FSA::dunnTest()]
#' @inheritParams gly_scale
#' @section SE identifiers:
#' For SE inputs, the result also contains an `identifiers` mapping of result
#' IDs to input positions and original names; see [gly_identity].
#' In SE tidy results, identity columns are `.variable`/`.sample` (dotted
#' endpoints for correlations); undotted metadata columns remain annotations.
#' @seealso [gly_scale] for input scales, transformations, and fold changes.
#'   [gly_identity] for SE identifiers, annotation alignment, and unnamed inputs.
#' @export
gly_kruskal <- function(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  add_info = TRUE,
  ...,
  subject_col = NULL,
  log_transform = NULL,
  contrasts = NULL,
  subset = NULL
) {
  # Validate inputs
  checkmate::assert_logical(
    log_transform,
    len = 1,
    any.missing = FALSE,
    null.ok = TRUE
  )
  .assert_data_container(exp)
  meta_data <- .get_meta_data(exp)
  transformation <- .resolve_transformation(
    exp,
    log_transform,
    default = TRUE,
    method = "log2",
    meta_data = meta_data
  )
  log_transform <- transformation$log_transform
  checkmate::assert_string(group_col)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_string(subject_col, null.ok = TRUE)

  # Extract data from experiment object
  input <- .multigroup_input(exp, subset, rlang::list2(...), rank_test = TRUE)
  expr_mat <- input$expr_mat
  sample_info <- input$sample_info

  # Extract and validate groups
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "kruskal",
    show_info = !is.character(contrasts)
  )
  groups <- group_info$groups
  subjects <- .extract_paired_subjects(sample_info, subject_col, group_col)

  # Run the internal computation
  result <- .analyze_kruskal(
    expr_mat,
    groups,
    p_adj_method,
    !!!input$dots,
    subjects = subjects,
    contrasts = contrasts,
    log_transform = log_transform,
    input_scale = transformation$input_scale
  )
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  # Add meta_data from experiment
  result$meta_data <- meta_data
  result$preprocessing <- transformation
  .record_se_identity(result, exp)
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
  ...,
  subjects = NULL,
  log_transform = TRUE,
  input_scale = "raw",
  contrasts = NULL
) {
  checkmate::assert_flag(log_transform)
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

  comparisons <- .parse_posthoc_comparisons(contrasts, groups)

  design <- .match_paired_samples(
    expr_mat,
    groups,
    subjects,
    method = "rank test"
  )
  if (!design$paired) {
    rlang::check_installed("FSA")
  }
  expr_mat <- design$expr_mat
  groups <- design$groups

  # Prepare data
  log_expr_mat <- .log_transform_expr_mat(expr_mat, log_transform)
  data <- .prepare_multi_group_data(
    log_expr_mat,
    groups,
    subjects = design$subjects,
    is_logged = TRUE
  )
  test_function <- if (design$paired) {
    stats::friedman.test
  } else {
    stats::kruskal.test
  }
  test_formula <- if (design$paired) {
    log_value ~ group | subject
  } else {
    NULL
  }

  # Generate raw results
  raw_main_test <- .generate_raw_main_results(
    data,
    test_function,
    formula = test_formula,
    ...,
    paired = design$paired
  )
  main_test <- .tibblify_main_test_results(
    raw_main_test,
    test_function,
    p_adj_method,
    effect_size_method = if (design$paired) {
      "kendalls_w"
    } else {
      "epsilon_squared"
    },
    expr_mat = log_expr_mat,
    groups = groups,
    effect_models = raw_main_test
  )
  raw_post_hoc_test <- .generate_raw_posthoc_results(
    raw_main_test,
    data,
    test_function,
    main_test,
    posthoc_method = if (design$paired) "paired_wilcox" else "dunn",
    paired = design$paired
  )
  post_hoc_test <- .tibblify_posthoc_results(
    raw_post_hoc_test,
    stats::kruskal.test
  )
  post_hoc_test <- .add_fold_change(
    post_hoc_test,
    expr_mat,
    groups,
    paired_data = if (design$paired) data else NULL,
    input_scale = input_scale
  )

  tidy_result <- .report_posthoc(
    main_test,
    post_hoc_test,
    comparisons,
    levels(groups)
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
  main_model,
  posthoc_method = "tukey",
  paired = FALSE
) {
  if (paired) {
    data_nested <- .complete_paired_variable_data(data_nested)
  }
  if (identical(posthoc_method, "paired_t")) {
    return(.paired_posthoc_test(data_nested, method = "t"))
  }
  if (identical(posthoc_method, "paired_wilcox")) {
    return(.paired_posthoc_test(data_nested, method = "wilcox"))
  }

  if (identical(.f, stats::aov)) {
    # Post-hoc for ANOVA/ANCOVA
    if (identical(posthoc_method, "emmeans")) {
      emm <- emmeans::emmeans(main_model, specs = "group")
      contrast <- emmeans::contrast(
        emm,
        method = "pairwise",
        adjust = "tukey",
        reverse = FALSE
      )
      result <- as.data.frame(summary(contrast))
      # emmeans wraps some group labels in parentheses when printing contrasts.
      # Use its marginal-mean grid to recover the original pairwise endpoints.
      pairs <- utils::combn(as.character(as.data.frame(emm)$group), 2L)
      result$ref_group <- pairs[1, ]
      result$test_group <- pairs[2, ]
      result
    } else {
      tukey <- suppressWarnings(stats::TukeyHSD(main_model, which = "group"))
      # Tukey orders rows by factor level. Keep the original names, which may
      # themselves contain the hyphen used in its printed comparison labels.
      attr(tukey, "group_pairs") <- utils::combn(main_model$xlevels$group, 2L)
      tukey
    }
  } else {
    # Dunn test for Kruskal-Wallis
    # Check number of groups - FSA::dunnTest fails with only 2 groups
    observed_groups <- droplevels(data_nested$group[
      !is.na(data_nested$log_value)
    ])
    group_levels <- levels(observed_groups)
    n_groups <- length(group_levels)

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
      observed <- which(!is.na(p_matrix), arr.ind = TRUE)
      ref <- colnames(p_matrix)[observed[, 2]]
      test <- rownames(p_matrix)[observed[, 1]]
      list(
        res = data.frame(
          Comparison = stringr::str_c(test, " - ", ref),
          ref_group = ref,
          test_group = test,
          P.unadj = p_matrix[observed],
          P.adj = p_matrix[observed]
        ),
        method = wilcox_result$p.adjust.method
      )
    } else {
      # Use FSA::dunnTest for 3+ groups
      .standardize_dunn_posthoc_result(
        FSA::dunnTest(log_value ~ group, data = data_nested, method = "holm"),
        group_levels
      )
    }
  }
}

# Helper function to generate raw main test results
.generate_raw_main_results <- function(
  data,
  .f,
  formula = NULL,
  ...,
  paired = FALSE
) {
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
  fit_one <- function(df) {
    args <- dots
    if (paired) {
      complete <- .complete_paired_variable_data(df)
      args <- .index_sample_args(
        args,
        match(complete$sample, df$sample),
        nrow(df)
      )
      df <- complete
    }
    rlang::exec(.f, formula, data = df, !!!args)
  }
  safe_fit <- purrr::possibly(fit_one, otherwise = NA)
  main_test_raw <- data |>
    dplyr::nest_by(.data$variable) |>
    dplyr::mutate(test_result = list(safe_fit(.data$data)))
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
  main_test,
  posthoc_method = "tukey",
  paired = FALSE
) {
  p <- main_test[["p_adj"]] %||% main_test$p_val
  sig_vars <- main_test$variable[!is.na(p) & p < 0.05]
  if (length(sig_vars) > 0) {
    posthoc_raw_results <- data %>%
      dplyr::filter(.data$variable %in% sig_vars) %>%
      dplyr::nest_by(.data$variable) %>%
      dplyr::mutate(
        posthoc_raw = list(.perform_raw_posthoc_test(
          .data$data,
          .f,
          main_model = main_test_raw[[.data$variable]],
          posthoc_method = posthoc_method,
          paired = paired
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
  subjects = NULL,
  is_logged = FALSE,
  log_transform = TRUE
) {
  checkmate::assert_flag(log_transform)
  # Downstream fits and effect sizes identify features by name. Retain the
  # former wide-table rejection instead of pooling distinct rows together.
  if (anyDuplicated(rownames(expr_mat)) > 0L) {
    cli::cli_abort("Assay row names must be unique for differential testing.")
  }
  if (!is_logged) {
    expr_mat <- .log_transform_expr_mat(expr_mat, log_transform)
  }

  # Keep feature identifiers as values, never as temporary column names.
  # Otherwise names such as "sample", "group", or a covariate can collide
  # with the analysis columns and be dropped, overwritten, or repaired.
  sample_indices <- rep(seq_len(ncol(expr_mat)), each = nrow(expr_mat))
  sample_names <- colnames(expr_mat)
  if (is.null(sample_names)) {
    sample_names <- as.character(seq_len(ncol(expr_mat)))
  }
  variable_names <- rownames(expr_mat)
  if (is.null(variable_names)) {
    variable_names <- paste0("V", seq_len(nrow(expr_mat)))
  }
  data <- tibble::tibble(
    sample = sample_names[sample_indices],
    group = groups[sample_indices],
    variable = rep(variable_names, times = ncol(expr_mat)),
    log_value = as.vector(expr_mat)
  )

  if (!is.null(subjects)) {
    data$subject <- subjects[sample_indices]
  }

  if (!is.null(covariates)) {
    covariate_tbl <- tibble::as_tibble(covariates)
    data <- dplyr::bind_cols(
      data,
      covariate_tbl[sample_indices, , drop = FALSE]
    )
  }

  data
}

# Keep complete subject blocks for one variable in a paired multi-group test
.complete_paired_variable_data <- function(data) {
  if (!"subject" %in% colnames(data)) {
    return(data)
  }

  group_levels <- levels(data$group)
  subject_levels <- unique(as.character(data$subject))
  complete <- purrr::map_lgl(subject_levels, function(subject) {
    subject_data <- data[as.character(data$subject) == subject, , drop = FALSE]
    identical(
      sort(unique(as.character(subject_data$group))),
      sort(group_levels)
    ) &&
      all(is.finite(subject_data$log_value))
  })
  complete_subjects <- subject_levels[complete]
  result <- data[
    as.character(data$subject) %in% complete_subjects,
    ,
    drop = FALSE
  ]
  result$subject <- droplevels(result$subject)
  result
}

# Run paired post-hoc tests with Holm correction
.paired_posthoc_test <- function(data, method) {
  checkmate::assert_choice(method, c("t", "wilcox"))
  group_levels <- levels(data$group)
  comparisons <- .make_comparisons(group_levels, reverse = TRUE)
  p_values <- purrr::map_dbl(comparisons, function(comparison) {
    test <- data$log_value[data$group == comparison[1]]
    ref <- data$log_value[data$group == comparison[2]]
    if (identical(method, "t")) {
      tryCatch(
        stats::t.test(test, ref, paired = TRUE)$p.value,
        error = function(cnd) NA_real_
      )
    } else {
      if (all(test == ref)) {
        return(NA_real_)
      }
      suppressWarnings(
        stats::wilcox.test(test, ref, paired = TRUE, exact = FALSE)$p.value
      )
    }
  })
  adjusted <- stats::p.adjust(p_values, method = "holm")
  ref_groups <- purrr::map_chr(comparisons, ~ .x[2])
  test_groups <- purrr::map_chr(comparisons, ~ .x[1])

  if (identical(method, "t")) {
    return(data.frame(
      contrast = paste(ref_groups, "-", test_groups),
      ref_group = ref_groups,
      test_group = test_groups,
      p.unadj = p_values,
      p.value = adjusted,
      stringsAsFactors = FALSE
    ))
  }

  list(
    res = data.frame(
      Comparison = paste(test_groups, "-", ref_groups),
      ref_group = ref_groups,
      test_group = test_groups,
      P.unadj = p_values,
      P.adj = adjusted,
      stringsAsFactors = FALSE
    ),
    method = "holm"
  )
}

# Helper function to convert main test raw results to tibble
.tibblify_main_test_results <- function(
  main_test_raw,
  .f,
  p_adj_method,
  effect_size_method = NULL,
  expr_mat = NULL,
  groups = NULL,
  effect_models = NULL
) {
  var_names <- names(main_test_raw)
  all_failed <- all(purrr::map_lgl(main_test_raw, rlang::is_na))

  if (all_failed && identical(.f, stats::aov)) {
    result_tbl <- tibble::tibble(
      variable = var_names,
      term = NA_character_,
      df = NA_real_,
      sumsq = NA_real_,
      meansq = NA_real_,
      statistic = NA_real_,
      p_val = NA_real_
    )
  } else if (
    all_failed &&
      (identical(.f, stats::kruskal.test) ||
        identical(.f, stats::friedman.test))
  ) {
    result_tbl <- tibble::tibble(
      variable = var_names,
      statistic = NA_real_,
      p_val = NA_real_,
      parameter = NA_integer_,
      method = NA_character_
    )
  } else {
    result_tbl <- tibble::tibble(
      variable = var_names,
      test_result = main_test_raw
    ) %>%
      dplyr::mutate(
        params = purrr::map(
          .data$test_result,
          function(fit) {
            if (rlang::is_na(fit)) {
              return(NULL)
            }
            terms <- broom::tidy(fit)
            if (identical(.f, stats::aov)) {
              # Keep one row even if the group term is aliased. Saturated
              # models have sums of squares but no F statistic or p-value.
              terms <- terms[match("group", terms$term), , drop = FALSE]
              if (!"statistic" %in% names(terms)) {
                terms$statistic <- NA_real_
                terms$p.value <- NA_real_
              }
            }
            terms
          }
        )
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
      groups,
      effect_models = effect_models
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
  groups,
  effect_models = NULL
) {
  checkmate::assert_choice(
    effect_size_method,
    c(
      "eta_squared",
      "epsilon_squared",
      "partial_eta_squared",
      "kendalls_w"
    )
  )

  effect_size_values <- purrr::map2_dbl(
    result_tbl$variable,
    result_tbl$statistic,
    function(var_name, statistic) {
      if (effect_size_method %in% c("eta_squared", "partial_eta_squared")) {
        .calculate_partial_eta_squared(effect_models[[var_name]])
      } else if (identical(effect_size_method, "epsilon_squared")) {
        .calculate_epsilon_squared(expr_mat, groups, var_name, statistic)
      } else {
        .calculate_kendalls_w(expr_mat, groups, var_name, statistic)
      }
    }
  )

  result_tbl$effect_size <- effect_size_values
  result_tbl
}

# Calculate partial eta-squared for the group term in a blocked ANOVA
.calculate_partial_eta_squared <- function(model) {
  if (rlang::is_na(model)) {
    return(NA_real_)
  }
  anova_table <- summary(model)[[1]]
  group_ss <- anova_table["group", "Sum Sq"]
  residual_ss <- anova_table["Residuals", "Sum Sq"]
  denominator <- group_ss + residual_ss
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  group_ss / denominator
}

# Calculate Kendall's W for a Friedman test
.calculate_kendalls_w <- function(expr_mat, groups, var_name, statistic) {
  if (!is.finite(statistic)) {
    return(NA_real_)
  }
  values <- purrr::map(
    levels(groups),
    ~ expr_mat[var_name, groups == .x]
  )
  complete <- Reduce(`&`, purrr::map(values, is.finite))
  n_subjects <- sum(complete)
  n_groups <- length(values)
  if (n_subjects < 2 || n_groups < 2) {
    return(NA_real_)
  }
  statistic / (n_subjects * (n_groups - 1))
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
      if (is.data.frame(raw_result)) {
        tibble::tibble(
          variable = var_name,
          ref_group = raw_result$ref_group,
          test_group = raw_result$test_group,
          p_val = raw_result$p.unadj %||% raw_result$p.value,
          p_adj = raw_result$p.value
        )
      } else {
        tukey_df <- as.data.frame(raw_result$group)
        group_pairs <- attr(raw_result, "group_pairs")
        tibble::tibble(
          variable = var_name,
          ref_group = group_pairs[1, ],
          test_group = group_pairs[2, ],
          p_val = tukey_df$`p adj`, # TukeyHSD already provides adjusted p-values
          p_adj = tukey_df$`p adj`
        )
      }
    } else {
      # For Dunn test results
      dunn_df <- raw_result$res

      tibble::tibble(
        variable = var_name,
        ref_group = dunn_df$ref_group,
        test_group = dunn_df$test_group,
        p_val = dunn_df$P.unadj, # Unadjusted p-value
        p_adj = dunn_df$P.adj # Adjusted p-value
      )
    }
  })

  # Combine all results
  dplyr::bind_rows(result_list)
}

.add_fold_change <- function(
  post_hoc_test,
  expr_mat,
  groups,
  paired_data = NULL,
  input_scale = "raw"
) {
  if (!is.null(paired_data)) {
    if (nrow(post_hoc_test) == 0) {
      post_hoc_test$log2fc <- numeric(0)
      return(post_hoc_test)
    }

    fc_res <- purrr::map_dfr(
      unique(post_hoc_test$variable),
      function(variable) {
        variable_data <- paired_data[
          paired_data$variable == variable,
          ,
          drop = FALSE
        ]
        complete_data <- .complete_paired_variable_data(variable_data)
        sample_idx <- match(complete_data$sample, colnames(expr_mat))
        variable_expr <- expr_mat[variable, sample_idx, drop = FALSE]
        variable_groups <- factor(
          complete_data$group,
          levels = levels(groups)
        )

        if (length(levels(groups)) == 2) {
          .fc_2groups(variable_expr, variable_groups, input_scale = input_scale)
        } else {
          .fc_multi_groups(
            variable_expr,
            variable_groups,
            input_scale = input_scale
          )
        }
      }
    )

    join_by <- if (length(levels(groups)) == 2) {
      "variable"
    } else {
      c("ref_group", "test_group", "variable")
    }
    return(dplyr::left_join(post_hoc_test, fc_res, by = join_by))
  }

  if (length(levels(groups)) == 2) {
    fc_res <- .fc_2groups(expr_mat, groups, input_scale = input_scale)
    return(dplyr::left_join(post_hoc_test, fc_res, by = "variable"))
  } else {
    fc_res <- .fc_multi_groups(expr_mat, groups, input_scale = input_scale)
    return(dplyr::left_join(
      post_hoc_test,
      fc_res,
      by = c("ref_group", "test_group", "variable")
    ))
  }
}

.standardize_dunn_posthoc_result <- function(raw_result, group_levels) {
  # Dunn lists all comparisons against each successive factor level. Match
  # whole labels because a level can itself contain the printed separator.
  indices <- seq_along(group_levels)
  pairs <- which(outer(indices, indices, `<`), arr.ind = TRUE)
  ref <- group_levels[pairs[, 1]]
  test <- group_levels[pairs[, 2]]
  rows <- match(
    make.unique(raw_result$res$Comparison),
    make.unique(paste(ref, "-", test))
  )
  if (anyNA(rows)) {
    cli::cli_abort(
      "Could not match Dunn comparisons to the observed group levels."
    )
  }
  raw_result$res$ref_group <- ref[rows]
  raw_result$res$test_group <- test[rows]
  raw_result$res$Comparison <- paste(test[rows], "-", ref[rows])
  raw_result$res$Z <- -raw_result$res$Z
  raw_result
}
