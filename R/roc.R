#' ROC Analysis for Glycomics and Glycoproteomics Data
#'
#' Perform Receiver Operating Characteristic (ROC) analysis for binary classification
#' of glycomics or glycoproteomics data. The function calculates ROC curves and
#' Area Under the Curve (AUC) values for each variable to assess their discriminatory
#' power between two groups.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing an expression matrix and
#'   sample information.
#' @param group_col A character string specifying the column name of the grouping variable
#'   in the sample information. Default is `"group"`. The grouping variable must have
#'   exactly 2 levels for binary classification.
#' @param pos_class A character string specifying which group level should be treated as
#'   the positive class. If `NULL` (default), the second level (alphabetically) will be
#'   used as the positive class.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the ROC analysis results are returned.
#'
#' @details
#' For each variable, a ROC curve is computed using the expression values as predictor
#' and the binary group labels as response.
#'
#' The function requires exactly 2 groups in the specified grouping variable. If more than
#' 2 groups are present, an error will be thrown.
#'
#' **Underlying Function:**
#' - ROC analysis is performed using `pROC::roc()`
#' - Coordinates are extracted using `pROC::coords()`
#'
#' @section Required packages:
#' This function requires the `pROC` package to be installed for ROC curve computation.
#'
#' @returns
#' A list with three elements:
#' - `tidy_result`: A list containing two tibbles:
#'   - `auc`: A tibble containing AUC values for each variable with the following columns:
#'     - `variable`: Variable name
#'     - `auc`: Area Under the Curve value
#'     - `auc_ci_low`: Lower bound of 95% confidence interval for AUC
#'     - `auc_ci_high`: Upper bound of 95% confidence interval for AUC
#'   - `coords`: A tibble containing ROC curve coordinates with the following columns:
#'     - `variable`: Variable name
#'     - `threshold`: Threshold value for classification
#'     - `specificity`: Specificity (True Negative Rate)
#'     - `sensitivity`: Sensitivity (True Positive Rate)
#' - `raw_result`: A list of `pROC` objects
#' - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_roc_res` and `glystats_res`.
#' @seealso [pROC::roc()], [pROC::coords()]
#' @export
gly_roc <- function(
  exp,
  group_col = "group",
  pos_class = NULL,
  add_info = TRUE
) {
  rlang::check_installed("pROC")

  # Validate inputs
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_string(pos_class, null.ok = TRUE)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)

  # Extract and validate groups using helper function
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "ROC analysis",
    show_info = TRUE
  )
  groups <- group_info$groups

  # Run the internal computation
  result <- .analyze_roc(expr_mat, groups, pos_class)
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )

  # Add meta_data from experiment
  result$meta_data <- .get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_roc()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_roc <- function(expr_mat, groups, pos_class = NULL) {
  rlang::check_installed("pROC")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_string(pos_class, null.ok = TRUE)

  # Validate group count
  if (length(levels(groups)) != 2) {
    cli::cli_abort("groups must have exactly 2 levels for ROC analysis")
  }

  # Set positive class
  if (is.null(pos_class)) {
    pos_class <- levels(groups)[2] # Use second level alphabetically as default
    cli::cli_alert_info("Using {.val {pos_class}} as positive class")
  } else {
    if (!pos_class %in% levels(groups)) {
      cli::cli_abort(c(
        "Positive class {.val {pos_class}} not found in group levels",
        "i" = "Available levels: {.val {levels(groups)}}"
      ))
    }
  }

  # Prepare data for ROC analysis
  response <- as.numeric(groups == pos_class)
  safe_f <- purrr::possibly(pROC::roc, otherwise = NA)
  roc_objs <- purrr::map(
    rownames(expr_mat),
    ~ suppressMessages(safe_f(response, expr_mat[.x, ]))
  )
  names(roc_objs) <- rownames(expr_mat)
  n_na <- sum(is.na(roc_objs))
  if (n_na > 0) {
    cli::cli_warn("{.val {n_na}} variables failed to fit the model")
  }

  # Create tidy results
  safe_ci <- purrr::possibly(pROC::ci.auc, otherwise = NA)
  auc_ci <- purrr::map(roc_objs, function(roc_obj) {
    if (rlang::is_na(roc_obj)) {
      return(c(auc_ci_low = NA_real_, auc_ci_high = NA_real_))
    }
    ci <- suppressWarnings(suppressMessages(safe_ci(
      roc_obj,
      conf.level = 0.95
    )))
    if (rlang::is_na(ci)) {
      return(c(auc_ci_low = NA_real_, auc_ci_high = NA_real_))
    }
    ci <- as.numeric(ci)
    c(auc_ci_low = ci[1], auc_ci_high = ci[3])
  })
  roc_auc_tb <- tibble::tibble(
    variable = rownames(expr_mat),
    auc = purrr::map_dbl(roc_objs, ~ if (rlang::is_na(.x)) NA else .x$auc),
    auc_ci_low = purrr::map_dbl(auc_ci, 1),
    auc_ci_high = purrr::map_dbl(auc_ci, 2)
  )

  coords_tb <- roc_objs %>%
    purrr::map(
      ~ tibble::as_tibble(
        if (rlang::is_na(.x)) NULL else pROC::coords(.x, "all")
      )
    ) %>%
    rlang::set_names(rownames(expr_mat)) %>%
    dplyr::bind_rows(.id = "variable")

  # Create tidy result list
  tidy_result <- list(auc = roc_auc_tb, coords = coords_tb)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = roc_objs
    ),
    class = c("glystats_roc_res", "glystats_res")
  )
}
