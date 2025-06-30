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
#' - `auc`: A named numeric vector of AUC values for each variable
#' - `thresholds`: A tibble containing ROC curve coordinates with columns:
#'   - `variable`: Variable name
#'   - `threshold`: Threshold value
#'   - `sensitivity`: Sensitivity (True Positive Rate)
#'   - `specificity`: Specificity (True Negative Rate)
#'
#' @examples
#' \dontrun{
#' # Create example experiment object with 2 groups
#' sample_info <- tibble::tibble(
#'   sample = paste0("S", 1:20),
#'   group = rep(c("Control", "Treatment"), each = 10)
#' )
#' 
#' var_info <- tibble::tibble(
#'   variable = paste0("Glycan_", 1:100)
#' )
#' 
#' expr_mat <- matrix(rnorm(2000), nrow = 100, ncol = 20)
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
#' # Perform ROC analysis
#' result <- gly_roc(exp, group_col = "group", pos_class = "Treatment")
#' 
#' # View AUC values
#' head(result$auc)
#' 
#' # View ROC curve coordinates for first variable
#' result$thresholds |> 
#'   dplyr::filter(variable == names(result$auc)[1])
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#'
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

  # Check if group column exists
  if (!group_col %in% colnames(sample_info)) {
    cli::cli_abort("Column {.field {group_col}} not found in sample information")
  }

  # Get group information
  groups <- sample_info[[group_col]]
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }

  # Check that we have exactly 2 groups
  n_groups <- length(levels(groups))
  if (n_groups != 2) {
    cli::cli_abort(c(
      "{.field {group_col}} must be a factor with exactly 2 levels for ROC analysis",
      "i" = "Current levels: {.val {levels(groups)}}"
    ))
  }

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

  cli::cli_alert_info("Group 1 (negative): {.val {levels(groups)[levels(groups) != pos_class]}}")
  cli::cli_alert_info("Group 2 (positive): {.val {pos_class}}")

  # Prepare data for ROC analysis
  # Convert to binary response (1 for positive class, 0 for negative class)
  response <- as.numeric(groups == pos_class)
  
  # Calculate ROC for each variable
  var_names <- rownames(expr_mat)
  n_vars <- length(var_names)
  
  cli::cli_alert_info("Performing ROC analysis for {.val {n_vars}} variables")
  
  # Helper function to compute ROC for a single variable
  .compute_roc_single <- function(var_name, predictor, return_raw = FALSE) {
    # Compute ROC curve
    roc_obj <- pROC::roc(response, predictor, quiet = TRUE)
    
    # Return raw ROC object if requested
    if (return_raw) {
      return(roc_obj)
    }
    
    # Extract coordinates
    coords <- pROC::coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
    
    # Return list with AUC and coordinates
    list(
      auc = as.numeric(roc_obj$auc),
      coords = tibble::tibble(
        variable = var_name,
        threshold = coords$threshold,
        sensitivity = coords$sensitivity,
        specificity = coords$specificity
      )
    )
  }
  
  # Use purrr to compute ROC for all variables
  roc_results <- purrr::map2(var_names, asplit(expr_mat, 1), 
                             ~ .compute_roc_single(.x, .y, return_raw))
  
  # Return raw results if requested
  if (return_raw) {
    names(roc_results) <- var_names
    return(roc_results)
  }
  
  # Extract AUC values using purrr
  auc_values <- purrr::map_dbl(roc_results, ~ .x$auc)
  names(auc_values) <- var_names
  
  # Extract and combine threshold data using purrr
  thresholds <- purrr::map_dfr(roc_results, ~ .x$coords)
  
  # Return results with S3 class
  result <- list(
    auc = auc_values,
    thresholds = thresholds
  )
  
  structure(result, class = c("glystats_roc_res", "glystats_res", class(result)))
}