#' Partial Least Squares Discriminant Analysis (PLS-DA)
#'
#' Perform partial least squares discriminant analysis on the expression data.
#' The function uses `mixOmics::plsda()` to perform PLS-DA and extracts relevant results
#' for visualization and interpretation.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param group_col Character string specifying the column name in sample information 
#'   that contains group labels for discrimination.
#'   Default to "group".
#' @param ncomp Integer, number of components to extract. Default is 2.
#' @param center A logical indicating whether to center the data. Default is TRUE.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param ... Additional arguments passed to `mixOmics::plsda()`.
#'
#' @return A list containing three tibbles:
#'  - `samples`: PLS-DA scores for each sample with group information
#'  - `variables`: PLS-DA loadings for each variable
#'  - `components`: Component information including explained variance
#' 
#' @examples
#' \dontrun{
#' plsda_res <- gly_plsda(exp)
#' }
#'
#' @export
gly_plsda <- function(exp, group_col = "group", ncomp = 2, center = TRUE, scale = TRUE, ...) {
  
  # Check if mixOmics is available
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg mixOmics} is required for PLS-DA analysis.",
      "i" = "Install it with: {.code install.packages('mixOmics')}"
    ))
  }
  
  # Extract expression matrix and sample info
  mat <- t(exp$expr_mat)  # Samples as rows, variables as columns
  sample_info <- exp$sample_info
  
  # Determine group variable
  if (is.null(group_col)) {
    char_cols <- names(sample_info)[sapply(sample_info, function(x) is.character(x) || is.factor(x))]
    if (length(char_cols) == 0) {
      cli::cli_abort("No character or factor columns found in sample information for grouping.")
    }
    group_col <- char_cols[1]
    cli::cli_inform("Using {.field {group_col}} as grouping variable.")
  }
  
  # Check if group variable exists
  if (!group_col %in% names(sample_info)) {
    cli::cli_abort("Group variable {.field {group_col}} not found in sample information.")
  }
  
  # Extract group labels
  group_labels <- sample_info[[group_col]]
  if (is.character(group_labels)) {
    group_labels <- as.factor(group_labels)
  }
  
  # Apply centering and scaling if requested
  if (center) {
    mat <- scale(mat, center = TRUE, scale = FALSE)
  }
  if (scale) {
    mat <- scale(mat, center = FALSE, scale = TRUE)
  }
  
  # Perform PLS-DA
  plsda_res <- mixOmics::plsda(
    X = mat, 
    Y = group_labels, 
    ncomp = ncomp,
    ...
  )
  
  # Use broom-style helper functions to extract tidy results
  res <- list(
    "samples" = .tidy_plsda_samples(plsda_res, mat, group_labels, group_col, ncomp),
    "variables" = .tidy_plsda_variables(plsda_res, ncomp),
    "components" = .tidy_plsda_components(plsda_res, ncomp)
  )
  
  structure(res, class = "glystats_plsda_res")
}

# Helper functions to tidy PLS-DA results (broom-style) ----------------------

# Tidy sample scores from PLS-DA results
.tidy_plsda_samples <- function(plsda_res, mat, group_labels, group_col, ncomp) {
  samples_tbl <- tibble::as_tibble(plsda_res$variates$X) %>%
    dplyr::mutate(
      sample = rownames(mat),
      .before = 1
    )
  names(samples_tbl)[2:(1+ncomp)] <- paste0("comp", 1:ncomp)
  samples_tbl
}

# Tidy variable loadings from PLS-DA results
.tidy_plsda_variables <- function(plsda_res, ncomp) {
  variables_tbl <- tibble::as_tibble(plsda_res$loadings$X, rownames = "variable")
  names(variables_tbl)[2:(1+ncomp)] <- paste0("comp", 1:ncomp)
  variables_tbl
}

# Tidy component information from PLS-DA results
.tidy_plsda_components <- function(plsda_res, ncomp) {
  explained_var <- plsda_res$prop_expl_var$X
  tibble::tibble(
    component = paste0("comp", 1:ncomp),
    explained_variance = explained_var[1:ncomp],
    cumulative_variance = cumsum(explained_var[1:ncomp])
  )
}
