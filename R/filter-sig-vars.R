#' Filter Significant Variables
#'
#' Filtering the experiment to keep only significant variables is a common task.
#' This function provides a convenient way to do this.
#' It supports results from all glystats DEA functions including
#' [gly_anova()], [gly_ancova()], [gly_kruskal()], [gly_ttest()], [gly_wilcox()], and [gly_limma()].
#'
#' @param exp An [glyexp::experiment()]. Please use the same experiment used to generate the DEA result.
#' @param res A glystats result object from a glystats DEA function.
#' @param p_adj_cutoff The threshold for p-adjusted values. Default is 0.05.
#' @param p_val_cutoff The threshold for p-values. We don't recommend using this. Default is NULL.
#'   If you insist to use it, please set `p_adj_cutoff` to NULL.
#' @param fc_cutoff The threshold for fold changes. Only positive value is needed.
#'   For example, `2` means fold change > 2 or < 1/2.
#'   Default is `2` for glycoproteomics data and `NULL` for others.
#' @param ... Additional arguments passed to methods. See the method-specific documentation for details.
#'
#' @returns An new [glyexp::experiment()] object.
#'
#' @examples
#' library(glyexp)
#' library(glyclean)
#'
#' exp <- auto_clean(real_experiment) |>
#'   glyexp::slice_head_var(n = 10)
#' res <- gly_anova(exp)
#' sig_exp <- filter_sig_vars(exp, res)
#'
#' @export
filter_sig_vars <- function(exp, res, p_adj_cutoff = 0.05, p_val_cutoff = NULL, fc_cutoff = NULL, ...) {
  UseMethod("filter_sig_vars", res)
}

#' @rdname filter_sig_vars
#' @param on (For [gly_anova()] and [gly_kruskal()] results only)
#'   "main_test" or "post_hoc_test". Should the filter be applied on the main test results or the post-hoc test results?
#'   Default is "main_test". If "post_hoc_test", please set a `comparison` value.
#' @param comparison (For [gly_anova()], [gly_kruskal()], and [gly_limma()] results only)
#'   Specifies which comparison to filter on. A string with the format "group1_vs_group2".
#'   For [gly_anova()] and [gly_kruskal()] results, `comparison` is only used when `on` is "post_hoc_test".
#'   If not provided, filtering will be performed on the main test results for [gly_anova()] and [gly_kruskal()],
#'   and variables will be kept if they are significant in any comparison for [gly_limma()].
#' @export
filter_sig_vars.glystats_anova_res <- function(exp, res, p_adj_cutoff = 0.05, p_val_cutoff = NULL, fc_cutoff = NULL, on = "main_test", comparison = NULL, ...) {
  .check_filter_sig_vars_args_anova(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison)
  fc_cutoff <- .decide_fc(exp, fc_cutoff)
  .filter_sig_vars_anova(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison)
}

#' @rdname filter_sig_vars
#' @export
filter_sig_vars.glystats_ancova_res <- function(exp, res, p_adj_cutoff = 0.05, p_val_cutoff = NULL, fc_cutoff = NULL, on = "main_test", comparison = NULL, ...) {
  .check_filter_sig_vars_args_anova(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison)
  fc_cutoff <- .decide_fc(exp, fc_cutoff)
  .filter_sig_vars_anova(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison)
}

#' @rdname filter_sig_vars
#' @export
filter_sig_vars.glystats_kruskal_res <- function(exp, res, p_adj_cutoff = 0.05, p_val_cutoff = NULL, fc_cutoff = NULL, on = "main_test", comparison = NULL, ...) {
  .check_filter_sig_vars_args_anova(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison)
  fc_cutoff <- .decide_fc(exp, fc_cutoff)
  .filter_sig_vars_anova(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison)
}

#' @rdname filter_sig_vars
#' @export
filter_sig_vars.glystats_ttest_res <- function(exp, res, p_adj_cutoff = 0.05, p_val_cutoff = NULL, fc_cutoff = NULL, ...) {
  .check_filter_sig_vars_args_ttest(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff)
  fc_cutoff <- .decide_fc(exp, fc_cutoff)
  .filter_sig_vars_ttest(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff)
}

#' @rdname filter_sig_vars
#' @export
filter_sig_vars.glystats_wilcox_res <- function(exp, res, p_adj_cutoff = 0.05, p_val_cutoff = NULL, fc_cutoff = NULL, ...) {
  .check_filter_sig_vars_args_ttest(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff)
  fc_cutoff <- .decide_fc(exp, fc_cutoff)
  .filter_sig_vars_ttest(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff)
}

#' @rdname filter_sig_vars
#' @export
filter_sig_vars.glystats_limma_res <- function(exp, res, p_adj_cutoff = 0.05, p_val_cutoff = NULL, fc_cutoff = NULL, comparison = NULL, ...) {
  .check_filter_sig_vars_args_limma(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison)
  fc_cutoff <- .decide_fc(exp, fc_cutoff)
  .filter_sig_vars_limma(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison)
}

.filter_sig_vars_anova <- function(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison) {
  res_df <- res$tidy_result[[on]]

  if (on == "main_test") {
    sig_df <- .add_filters(res_df, p_adj_cutoff, p_val_cutoff, NULL)
  } else { # post_hoc_test
    sig_df <- .add_filters(res_df, p_adj_cutoff, NULL, fc_cutoff)
    if (!is.null(comparison)) {
      sig_df <- .filter_comparison(sig_df, comparison)
    }
  }

  sig_vars <- sig_df |> dplyr::pull("variable") |> unique()
  glyexp::filter_var(exp, .data$variable %in% sig_vars)
}

.filter_sig_vars_ttest <- function(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison) {
  sig_vars <- res$tidy_result |>
    .add_filters(p_adj_cutoff, p_val_cutoff, fc_cutoff) |>
    dplyr::pull("variable") |>
    unique()
  glyexp::filter_var(exp, .data$variable %in% sig_vars)
}

.filter_sig_vars_limma <- function(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison) {
  df <- res$tidy_result
  if (!is.null(comparison)) {
    df <- .filter_comparison(df, comparison)
  }
  df <- .add_filters(df, p_adj_cutoff, p_val_cutoff, fc_cutoff)
  sig_vars <- df |> dplyr::pull("variable") |> unique()
  glyexp::filter_var(exp, .data$variable %in% sig_vars)
}

.add_filters <- function(df, p_adj_cutoff, p_val_cutoff, fc_cutoff) {
  if (!is.null(p_val_cutoff)) {
    df <- dplyr::filter(df, .data$p_val < p_val_cutoff)
  }
  if (!is.null(fc_cutoff)) {
    df <- dplyr::filter(df, abs(.data$log2fc) > log2(fc_cutoff))
  }
  if (!is.null(p_adj_cutoff)) {
    df <- dplyr::filter(df, .data$p_adj < p_adj_cutoff)
  }
  df
}

.filter_comparison <- function(df, comparison) {
  df |>
    dplyr::mutate(comparison = paste0(.data$ref_group, "_vs_", .data$test_group)) |>
    dplyr::filter(.data$comparison == .env$comparison) |>
    dplyr::select(-dplyr::all_of("comparison"))
}

.check_filter_sig_vars_args <- function(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison) {
  # Type checking
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_number(p_adj_cutoff, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_number(p_val_cutoff, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_number(fc_cutoff, lower = 1, null.ok = TRUE)
  checkmate::assert_string(comparison, null.ok = TRUE)

  # Value checking
  if (is.null(p_adj_cutoff) && is.null(p_val_cutoff) && is.null(fc_cutoff)) {
    cli::cli_abort(c(
      "At least one of {.arg p_adj_cutoff}, {.arg p_val_cutoff}, or {.arg fc_cutoff} must be provided.",
      "x" = "All are NULL."
    ))
  }
  if (!is.null(p_val_cutoff) && !is.null(p_adj_cutoff)) {
    cli::cli_abort(c(
      "Only one of {.arg p_adj_cutoff} or {.arg p_val_cutoff} can be provided.",
      "x" = "Both are provided."
    ))
  }

  # Deal with comparison
  if (!is.null(comparison)) {
    groups <- stringr::str_split_1(comparison, "_vs_")
    if (length(groups) != 2) {
      cli::cli_abort(c(
        "{.arg comparison} must be in the format of {.val group1_vs_group2}.",
        "x" = "Invalid format: {.val {comparison}}."
      ))
    }
  }
}

.check_filter_sig_vars_args_anova <- function(exp, res, p_adj_cutoff, p_val_cutoff, fc_cutoff, on, comparison) {
  .check_filter_sig_vars_args(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison)
  checkmate::assert_choice(on, c("main_test", "post_hoc_test"))
  if (!is.null(p_val_cutoff) && on == "post_hoc_test") {
    cli::cli_abort(c(
      "{.arg p_val_cutoff} can't be used with {.fn gly_anova} or {.fn gly_kruskal} when {.arg on} is {.val post_hoc_test}.",
      "i" = "This is because post-hoc test results only have adjusted p-values.",
      "i" = "Please use {.arg p_adj_cutoff} instead."
    ))
  }
  if (!is.null(comparison) && on == "main_test") {
    cli::cli_abort(c(
      "{.arg comparison} can't be used with {.fn gly_anova} or {.fn gly_kruskal} when {.arg on} is {.val main_test}.",
      "i" = "Please set {.arg comparison} to NULL or {.arg on} to {.val post_hoc_test}."
    ))
  }
  if (!is.null(fc_cutoff) && on == "main_test") {
    cli::cli_abort(c(
      "{.arg fc_cutoff} can't be used with {.fn gly_anova} or {.fn gly_kruskal} when {.arg on} is {.val main_test}.",
      "i" = "Please set {.arg fc_cutoff} to NULL or {.arg on} to {.val post_hoc_test}."
    ))
  }
  if (is.null(comparison) && on == "post_hoc_test") {
    cli::cli_abort(c(
      "{.arg comparison} must be provided when {.arg on} is {.val post_hoc_test}.",
      "i" = "Please set {.arg comparison} to a string with the format of {.val group1_vs_group2}."
    ))
  }
  if (!is.null(comparison)) {
    unique_comparisons <- res$tidy_result$post_hoc_test |>
      dplyr::mutate(comparison = paste0(.data$ref_group, "_vs_", .data$test_group)) |>
      dplyr::pull("comparison") |>
      unique()
    if (!comparison %in% unique_comparisons) {
      maybe_correct <- local({
        groups <- stringr::str_split_1(comparison, "_vs_")
        paste0(groups[2], "_vs_", groups[1])
      })
      if (maybe_correct %in% unique_comparisons) {
        cli::cli_abort(c(
          "Can't find comparison: {.val {comparison}}.",
          "i" = "Available comparisons: {.val {unique_comparisons}}.",
          "i" = "Did you mean {.val {maybe_correct}}?"
        ))
      } else {
        cli::cli_abort(c(
          "Can't find comparison: {.val {comparison}}.",
          "i" = "Available comparisons: {.val {unique_comparisons}}."
        ))
      }
    }
  }
}

.check_filter_sig_vars_args_ttest <- function(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff) {
  .check_filter_sig_vars_args(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff, NULL)
}

.check_filter_sig_vars_args_limma <- function(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison) {
  .check_filter_sig_vars_args(exp, p_adj_cutoff, p_val_cutoff, fc_cutoff, comparison)
}

.decide_fc <- function(exp, user_provided) {
  if (!is.null(user_provided)) {
    return(user_provided)
  }
  if (glyexp::get_exp_type(exp) == "glycoproteomics") {
    return(2)
  }
  NULL
}
