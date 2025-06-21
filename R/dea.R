#' Differential Expression Analysis (DEA)
#'
#' Differential expression analysis for glycomics or glycoproteomics data.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param method A character string specifying the method to use for DEA.
#'  Either "t-test" or "wilcoxon".
#' @param group_col A character string specifying the column name of the group variable.
#'  Default is "group".
#' @param ... Additional arguments passed to the underlying function.
#'
#' @returns The DEA result as a tibble.
#'  If `method` is "t-test", a column `log2fc` is added to the result.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#' @export
gly_dea <- function(exp, method = NULL, group_col = "group", ...) {
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_choice(method, c("t-test", "wilcoxon"))
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
  if (length(levels(groups)) != 2) {
    cli::cli_abort(c(
      "{.field {group_col}} must be a factor with exactly 2 levels",
      "i" = "Current levels: {.val {levels(groups)}}"
    ))
  }
  cli::cli_alert_info("Group 1: {.val {levels(groups)[1]}}")
  cli::cli_alert_info("Group 2: {.val {levels(groups)[2]}}")

  # Perform DEA
  switch(method,
    "t-test" = .gly_dea_2groups(expr_mat, groups, t.test, ...),
    "wilcoxon" = .gly_dea_2groups(expr_mat, groups, wilcox.test, ...),
    cli::cli_abort("Invalid method: {.val {method}}")
  )
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
      ttest = list(.f(log_value ~ group, data = data)),
      params = list(parameters::model_parameters(ttest)),
      params = list(parameters::standardize_names(params)),
      ) %>%
    dplyr::select(all_of(c("variable", "params"))) %>%
    tidyr::unnest(all_of("params")) %>%
    dplyr::ungroup() %>%
    janitor::clean_names()

  if (identical(.f, t.test)) {
    ttest_res <- dplyr::mutate(ttest_res, log2fc = .data$mean_group1 - .data$mean_group2)
  }

  ttest_res
}
