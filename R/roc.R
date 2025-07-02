#' ROC Analysis for Glycomics and Glycoproteomics Data
#'
#' Perform Receiver Operating Characteristic (ROC) analysis for binary classification
#' of glycomics or glycoproteomics data. The function calculates ROC curves and 
#' Area Under the Curve (AUC) values for each variable to assess their discriminatory
#' power between two groups.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param group_col A character string specifying the column name of the grouping variable
#'   in the sample information. Default is `"group"`. The grouping variable must have
#'   exactly 2 levels for binary classification.
#' @param pos_class A character string specifying which group level should be treated as
#'   the positive class. If `NULL` (default), the second level (alphabetically) will be
#'   used as the positive class.
#' @param return_raw A logical value. If FALSE (default), returns processed results with
#'   AUC values and threshold coordinates. If TRUE, returns raw pROC objects as a list.
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
#' A list containing two elements:
#' - `auc`: A tibble containing AUC values for each variable, with columns:
#'   - `variable`: Variable name
#'   - `auc`: AUC value
#' - `coords`: A tibble containing ROC curve coordinates with columns:
#'   - `variable`: Variable name
#'   - `threshold`: Threshold value
#'   - `sensitivity`: Sensitivity (True Positive Rate)
#'   - `specificity`: Specificity (True Negative Rate)
#' If `return_raw` is TRUE, returns a list of `pROC` objects.
#' @seealso [pROC::roc()], [pROC::coords()]
#' @export
gly_roc <- function(exp, group_col = "group", pos_class = NULL, return_raw = FALSE) {
  .check_pkg_available("pROC")
  
  # Validate inputs
  checkmate::check_class(exp, "glyexp_experiment")
  checkmate::check_string(group_col)
  checkmate::check_string(pos_class, null.ok = TRUE)
  checkmate::check_logical(return_raw, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

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

  # Set positive class
  if (is.null(pos_class)) {
    pos_class <- levels(groups)[2]  # Use second level alphabetically as default
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
  # Convert to binary response (1 for positive class, 0 for negative class)
  response <- as.numeric(groups == pos_class)
  roc_objs <- purrr::map(
    rownames(expr_mat),
    ~ suppressMessages(pROC::roc(response, expr_mat[.x,]))
  )

  if (return_raw) {
    return(roc_objs)
  }

  roc_auc_tb <- tibble::tibble(
    variable = rownames(expr_mat),
    auc = purrr::map_dbl(roc_objs, ~ .x$auc)
  )

  coords_tb <- roc_objs %>%
    purrr::map(~ tibble::as_tibble(pROC::coords(.x, "all"))) %>%
    rlang::set_names(rownames(expr_mat)) %>%
    dplyr::bind_rows(.id = "variable")
  
  res <- list(auc = roc_auc_tb, coords = coords_tb)
  structure(res, class = c("glystats_roc_res", "glystats_res"))
}
