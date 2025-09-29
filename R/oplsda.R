#' Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#'
#' Perform orthogonal partial least squares discriminant analysis on the expression data.
#' The function uses `ropls::opls()` to perform OPLS-DA and returns tidy results.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param groups A factor or character vector specifying group membership for each sample.
#'   Must have exactly 2 levels. Character vectors will be automatically converted to factors.
#' @param group_col A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param pred_i An integer indicating the number of predictive components to include. Default is 1.
#' @param ortho_i An integer indicating the number of orthogonal components to include. Default is NA (automatic).
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the OPLS-DA results are returned.
#'  Only applicable to `gly_oplsda()`.
#' @param ... Additional arguments passed to `ropls::opls()`.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `ropls` for OPLS-DA analysis
#'
#' @details
#' `gly_oplsda()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_oplsda_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' @return A list containing:
#'  - `tidy_result`: A list of tibbles with OPLS-DA results:
#'    - `samples`: OPLS-DA scores for each sample containing the following columns:
#'      - `sample`: Sample name
#'      - `group`: Group assignment
#'      - `p1`, `p2`, etc.: Predictive component scores
#'      - `o1`, `o2`, etc.: Orthogonal component scores
#'    - `variables`: OPLS-DA loadings for each variable containing the following columns:
#'      - `variable`: Variable name
#'      - `p1`, `p2`, etc.: Predictive component loadings
#'      - `o1`, `o2`, etc.: Orthogonal component loadings
#'      - `pcorr1`, `pcorr2`, etc.: Correlation between each variable and the corresponding predictive component
#'    - `variance`: OPLS-DA explained variance containing the following columns:
#'      - `component`: Component name (p1, o1, etc.)
#'      - `prop_var_explained`: Proportion of variance explained by each component
#'      - `cumulative_prop_var`: Cumulative proportion of variance explained
#'    - `vip`: Variable Importance in Projection scores containing the following columns:
#'      - `variable`: Variable name
#'      - `vip`: VIP score
#'    - `perm_test`: Permutation test results containing the following columns:
#'      - `model`: Model type ("Original" for the original model, "Permutation" for permuted models)
#'      - `perm_id`: Permutation ID (0 for original model, 1+ for permutations)
#'      - Additional columns from the permutation test matrix (e.g., R2X, R2Y, Q2, etc.)
#'  - `raw_result`: The raw ropls opls object from `ropls::opls()`
#' @seealso [ropls::opls()]
#' @export
gly_oplsda <- function(exp, group_col = "group", pred_i = 1, ortho_i = NA, scale = TRUE, add_info = TRUE, ...) {
  # Check package availability
  rlang::check_installed("ropls")

  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_int(pred_i, lower = 1)
  if (!is.na(ortho_i)) {
    checkmate::assert_int(ortho_i, lower = 0)
  }
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups (OPLS-DA only supports binary classification)
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "OPLS-DA"
  )
  groups <- group_info$groups

  # Call the underlying API
  result <- gly_oplsda_(expr_mat, groups, pred_i, ortho_i, scale, ...)
  result$tidy_result <- .process_results_add_info(result$tidy_result, exp, add_info)
  result
}

#' @rdname gly_oplsda
#' @export
gly_oplsda_ <- function(expr_mat, groups, pred_i = 1, ortho_i = NA, scale = TRUE, ...) {
  rlang::check_installed("ropls")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  groups <- .convert_groups_to_factor(groups)
  checkmate::assert_factor(groups, len = ncol(expr_mat))

  # Prepare data matrix (samples as rows, variables as columns)
  mat <- log(t(expr_mat) + 1)

  # Perform OPLS-DA
  # Set appropriate cross-validation folds based on sample size
  n_samples <- nrow(mat)
  crossval_i <- min(7, n_samples - 1)  # Default is 7, but must be less than sample size

  # Suppress plotting to prevent Rplots.pdf generation
  # Open a null device to capture any plotting output
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # For small datasets, reduce permutation tests to allow model building
  perm_i <- if (n_samples < 10) 0 else 20

  # Handle ortho_i parameter - if NA, let ropls decide automatically
  if (is.na(ortho_i)) {
    oplsda_res <- ropls::opls(x = mat, y = groups, predI = pred_i,
                              scaleC = if (scale) "standard" else "none",
                              crossvalI = crossval_i, permI = perm_i,
                              fig.pdfC = "none", info.txtC = "none", ...)
  } else {
    oplsda_res <- ropls::opls(x = mat, y = groups, predI = pred_i, orthoI = ortho_i,
                              scaleC = if (scale) "standard" else "none",
                              crossvalI = crossval_i, permI = perm_i,
                              fig.pdfC = "none", info.txtC = "none", ...)
  }

  # Extract and format results
  tidy_result <- .format_oplsda_results(oplsda_res, groups, NULL, FALSE)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = oplsda_res
    ),
    class = c("glystats_oplsda_res", "glystats_res")
  )
}

# Helper function to format OPLS-DA results
# Extract sample scores from OPLS-DA results
.format_oplsda_samples <- function(oplsda_res) {
  # Extract sample scores (predictive + orthogonal components)
  # ropls stores scores in @scoreMN for predictive and @orthoScoreMN for orthogonal
  pred_scores <- oplsda_res@scoreMN
  ortho_scores <- oplsda_res@orthoScoreMN

  # Combine predictive and orthogonal scores
  if (!is.null(ortho_scores) && ncol(ortho_scores) > 0 && nrow(ortho_scores) == nrow(pred_scores)) {
    all_scores <- cbind(pred_scores, ortho_scores)
    # Create column names: p1, p2, ... for predictive, o1, o2, ... for orthogonal
    pred_names <- paste0("p", seq_len(ncol(pred_scores)))
    ortho_names <- paste0("o", seq_len(ncol(ortho_scores)))
    colnames(all_scores) <- c(pred_names, ortho_names)
  } else {
    all_scores <- pred_scores
    colnames(all_scores) <- paste0("p", seq_len(ncol(pred_scores)))
  }

  samples_tbl <- tibble::as_tibble(all_scores, .name_repair = "minimal")
  samples_tbl$sample <- rownames(all_scores)

  # Note: We don't add group here because it will be handled by .process_results_add_info
  # This ensures consistent behavior across all functions

  samples_tbl
}

# Extract variable loadings and p(corr) from OPLS-DA results
.format_oplsda_variables <- function(oplsda_res) {
  # Extract variable loadings (predictive + orthogonal components)
  pred_loadings <- oplsda_res@loadingMN
  ortho_loadings <- oplsda_res@orthoLoadingMN

  # Combine predictive and orthogonal loadings
  if (!is.null(ortho_loadings) && ncol(ortho_loadings) > 0 && nrow(ortho_loadings) == nrow(pred_loadings)) {
    all_loadings <- cbind(pred_loadings, ortho_loadings)
    # Create column names: p1, p2, ... for predictive, o1, o2, ... for orthogonal
    pred_names <- paste0("p", seq_len(ncol(pred_loadings)))
    ortho_names <- paste0("o", seq_len(ncol(ortho_loadings)))
    colnames(all_loadings) <- c(pred_names, ortho_names)
  } else {
    all_loadings <- pred_loadings
    colnames(all_loadings) <- paste0("p", seq_len(ncol(pred_loadings)))
  }

  variables_tbl <- tibble::as_tibble(all_loadings, .name_repair = "minimal")
  variables_tbl$variable <- rownames(all_loadings)

  # Calculate p(corr) for each predictive component
  # Get the modeling matrix x_model (centered/scaled), rows are samples, columns are variables
  x_model <- oplsda_res@suppLs$xModelMN

  # Ensure sample alignment between x_model and scores
  if (!is.null(rownames(x_model)) && !is.null(rownames(oplsda_res@scoreMN))) {
    x_model <- x_model[rownames(oplsda_res@scoreMN), , drop = FALSE]
  }

  # Calculate p(corr) for each predictive component
  pred_scores <- oplsda_res@scoreMN
  n_pred_comp <- ncol(pred_scores)

  for (i in seq_len(n_pred_comp)) {
    # Get scores for component i
    t_i <- pred_scores[, i, drop = TRUE]

    # Calculate correlation between each variable and component i
    pcorr_i <- apply(x_model, 2, function(x) stats::cor(x, t_i, use = "pairwise.complete.obs"))

    # Add pcorr column to variables table
    pcorr_col_name <- paste0("pcorr", i)
    variables_tbl[[pcorr_col_name]] <- as.numeric(pcorr_i[variables_tbl$variable])
  }

  variables_tbl
}

# Extract variance explained information from OPLS-DA results
.format_oplsda_variance <- function(oplsda_res) {
  # Extract explained variance information from modelDF
  # modelDF contains individual R2X values for each component
  model_df <- oplsda_res@modelDF

  # Remove the 'sum' row if present
  if ("sum" %in% rownames(model_df)) {
    model_df <- model_df[rownames(model_df) != "sum", , drop = FALSE]
  }

  # Get component names and R2X values
  component_names <- rownames(model_df)
  r2x_values <- model_df[["R2X"]]

  # Remove NA values
  valid_idx <- !is.na(r2x_values)
  component_names <- component_names[valid_idx]
  r2x_values <- r2x_values[valid_idx]

  tibble::tibble(
    component = component_names,
    prop_var_explained = r2x_values,
    cumulative_prop_var = cumsum(r2x_values)
  )
}

# Extract VIP scores from OPLS-DA results
.format_oplsda_vip <- function(oplsda_res) {
  # Calculate VIP (Variable Importance in Projection) scores
  # ropls provides VIP scores in @vipVn
  vip_scores <- oplsda_res@vipVn

  tibble::tibble(
    variable = names(vip_scores),
    vip = as.numeric(vip_scores)
  )
}

# Extract permutation test results from OPLS-DA results
.format_oplsda_perm_test <- function(oplsda_res) {
  # Extract permutation test results
  # ropls stores permutation results in @suppLs$permMN
  perm_m <- oplsda_res@suppLs$permMN

  if (!is.null(perm_m) && nrow(perm_m) > 0) {
    tibble::as_tibble(perm_m, .name_repair = "minimal") |>
      dplyr::mutate(
        perm_id = dplyr::row_number() - 1,
        model   = dplyr::if_else(.data$perm_id == 0, "Original", "Permutation")
      ) |>
      dplyr::relocate(all_of(c("model", "perm_id")))
  } else {
    # If no permutation test was performed, return empty tibble
    tibble::tibble(
      model = character(0),
      perm_id = integer(0)
    )
  }
}

.format_oplsda_results <- function(oplsda_res, groups, sample_info, add_info = TRUE) {
  # Check if model was successfully built
  if (length(oplsda_res@scoreMN) == 0) {
    # Model building failed - provide informative error
    cli::cli_abort(c(
      "OPLS-DA model building failed.",
      "i" = "This usually happens when:",
      "*" = "The data doesn't have sufficient discriminatory power between groups",
      "*" = "The first predictive component is not statistically significant",
      "*" = "Sample size is too small relative to the number of variables",
      "!" = "Consider using PLS-DA instead, or check your data quality."
    ))
  }

  list(
    samples = .format_oplsda_samples(oplsda_res),
    variables = .format_oplsda_variables(oplsda_res),
    variance = .format_oplsda_variance(oplsda_res),
    vip = .format_oplsda_vip(oplsda_res),
    perm_test = .format_oplsda_perm_test(oplsda_res)
  )
}
