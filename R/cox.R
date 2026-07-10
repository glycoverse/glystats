#' Cox Proportional Hazards Model for Survival Analysis
#'
#' Fit a Cox proportional hazards model for each variable in the expression data,
#' and extract p-values and hazard ratios from it.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param time_col A character string specifying the column name in sample information
#'   that contains survival time. Default is "time".
#' @param event_col A character string specifying the column name in sample information
#'   that contains event indicator (1 for event, 0 for censoring). Default is "event".
#' @param p_adj_method A character string specifying the method to adjust p-values.
#'   See `p.adjust.methods` for available methods. Default is "BH".
#'   If NULL, no adjustment is performed.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'   will be added to the result tibble. If FALSE, only the Cox model results are returned.
#' @param ... Additional arguments passed to `survival::coxph()`.
#'
#' @details
#' The function fits an Cox proportional hazards model for each variable by:
#'
#' ```r
#' coxph(Surv(time, event) ~ expression_value)
#' ````
#'
#' P-values are adjusted by Benjamini-Hochberg method by default.
#'
#' @returns A list with three elements:
#'  - `tidy_result`: A tibble with Cox model results containing the following columns:
#'    - `variable`: Variable name
#'    - `coefficient`: Regression coefficient (log hazard ratio)
#'    - `std.error`: Standard error of the coefficient
#'    - `statistic`: Wald test statistic
#'    - `p_val`: Raw p-value from Wald test
#'    - `hr`: Hazard ratio (exp(coefficient))
#'    - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'  - `raw_result`: A list of raw `coxph` model objects.
#'  - `meta_data`: A list containing metadata from the input experiment
#'
#' @seealso [survival::coxph()], [survival::Surv()]
#' @export
gly_cox <- function(
  exp,
  time_col = "time",
  event_col = "event",
  p_adj_method = "BH",
  add_info = TRUE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(time_col)
  checkmate::assert_string(event_col)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract time and event vectors
  time <- sample_info[[time_col]]
  event <- sample_info[[event_col]]

  # Run the internal computation
  result <- .analyze_cox(expr_mat, time, event, p_adj_method, ...)
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_cox()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_cox <- function(
  expr_mat,
  time,
  event,
  p_adj_method = "BH",
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, min.rows = 1, min.cols = 1)
  checkmate::assert_numeric(time, min.len = 1)
  checkmate::assert_numeric(event, min.len = 1)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )

  # Check dimensions match
  if (ncol(expr_mat) != length(time)) {
    cli::cli_abort(
      "Number of samples in expression matrix ({ncol(expr_mat)}) must match length of time vector ({length(time)})"
    )
  }
  if (ncol(expr_mat) != length(event)) {
    cli::cli_abort(
      "Number of samples in expression matrix ({ncol(expr_mat)}) must match length of event vector ({length(event)})"
    )
  }

  # Check if survival package is available
  rlang::check_installed("survival")

  # Fit Cox model for each variable
  safe_f <- purrr::possibly(survival::coxph, otherwise = NA)
  cox_models <- apply(expr_mat, 1, function(x) {
    safe_f(survival::Surv(time, event) ~ x, ...)
  })

  n_na <- sum(is.na(cox_models))
  if (n_na > 0) {
    cli::cli_warn("{.val {n_na}} variables failed to fit the model")
  }

  # Extract results and convert to tibble
  tidy_result <- .cox_tibblify(cox_models, p_adj_method)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = cox_models
    ),
    class = c("glystats_cox_res", "glystats_res")
  )
}

.cox_tibblify <- function(cox_models, p_adj_method) {
  # Extract results and convert to tibble
  result_df <- tibble::tibble(
    variable = names(cox_models),
    model = cox_models
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      params = list(
        if (rlang::is_na(.data$model)) NULL else broom::tidy(.data$model)
      )
    ) %>%
    dplyr::select(-all_of("model")) %>%
    tidyr::unnest(all_of("params")) %>%
    dplyr::ungroup() %>%
    # Rename columns to match expected names
    dplyr::rename(all_of(c(
      "coefficient" = "estimate",
      "p_val" = "p.value"
    ))) %>%
    # Calculate hazard ratio
    dplyr::mutate(hr = exp(.data$coefficient)) %>%
    # Remove the term column as it's redundant with variable
    dplyr::select(-all_of("term"))

  if (!is.null(p_adj_method)) {
    result_df <- dplyr::mutate(
      result_df,
      p_adj = stats::p.adjust(.data$p_val, method = p_adj_method)
    )
  }

  result_df
}
