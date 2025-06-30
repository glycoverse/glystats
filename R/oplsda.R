#' Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#'
#' Perform orthogonal partial least squares discriminant analysis on the expression data.
#' The function uses `ropls::opls()` to perform OPLS-DA and extracts relevant results
#' for visualization and interpretation. OPLS-DA separates the predictive variation from
#' the orthogonal variation that is not correlated with the response.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param group_col Character string specifying the column name in sample information 
#'   that contains group labels for discrimination.
#'   Default to "group".
#' @param predI Integer, number of predictive components. Default is 1.
#' @param orthoI Integer, number of orthogonal components. If NULL (default), 
#'   the number is automatically determined.
#' @param center A logical indicating whether to center the data. Default is TRUE.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw ropls object.
#' @param ... Additional arguments passed to `ropls::opls()`.
#'
#' @section Required packages:
#' This function requires the `ropls` package to be installed for OPLS-DA analysis.
#'
#' @return A list containing three tibbles (when return_raw = FALSE):
#'  - `samples`: OPLS-DA scores for each sample with group information
#'  - `variables`: OPLS-DA loadings for each variable
#'  - `components`: Component information including explained variance
#' When return_raw = TRUE, returns the raw ropls object.
#' 
#' @examples
#' \dontrun{
#' oplsda_res <- gly_oplsda(exp)
#' }
#'
#' @export
gly_oplsda <- function(exp, group_col = "group", predI = 1, orthoI = NULL, 
                       center = TRUE, scale = TRUE, return_raw = FALSE, ...) {
  
  .check_pkg_available("ropls")
  
  checkmate::check_logical(return_raw, len = 1)
  
  # Extract expression matrix and sample info
  mat <- t(exp$expr_mat)  # Samples as rows, variables as columns
  sample_info <- exp$sample_info
  
  # Extract and validate groups using helper function
  .check_group_column_exists(sample_info, group_col)
  group_labels <- .extract_groups(sample_info, group_col)
  
  # Apply centering and scaling if requested
  if (center) {
    mat <- scale(mat, center = TRUE, scale = FALSE)
  }
  if (scale) {
    mat <- scale(mat, center = FALSE, scale = TRUE)
  }
  
  # Determine scaling option for ropls
  scaling_option <- "none"
  if (center && scale) {
    scaling_option <- "standard"
  } else if (center) {
    scaling_option <- "center"
  } else if (scale) {
    scaling_option <- "pareto"
  }
  
  # Determine orthoI if NULL
  if (is.null(orthoI)) {
    orthoI <- NA_integer_  # Let ropls determine automatically
  }
  
  # Perform OPLS-DA (suppress output and graphics to avoid R CMD check issues)
  oplsda_res <- suppressMessages(suppressWarnings({
    invisible(utils::capture.output({
      oplsda_result <- ropls::opls(
        x = mat, 
        y = group_labels, 
        predI = predI,
        orthoI = orthoI,
        scaleC = scaling_option,
        crossvalI = min(7, nrow(mat) - 1),
        fig.pdfC = "none",  # Explicitly disable PDF output
        info.txtC = "none", # Disable text info output
        ...
      )
    }, type = "output"))
    oplsda_result
  }))
  
  # Return raw results if requested
  if (return_raw) {
    return(oplsda_res)
  }

  # Use broom-style helper functions to extract tidy results
  res <- list(
    "samples" = .tidy_oplsda_samples(oplsda_res, mat, group_labels, group_col, predI, orthoI),
    "variables" = .tidy_oplsda_variables(oplsda_res, predI, orthoI),
    "components" = .tidy_oplsda_components(oplsda_res, predI, orthoI)
  )
  
  structure(res, class = c("glystats_oplsda_res", "glystats_res"))
}

# Helper functions to tidy OPLS-DA results (broom-style) ---------------------

# Utility function to safely extract data from ropls object
.safe_extract <- function(expr, default = NULL) {
  tryCatch(expr, error = function(e) default)
}

# Utility function to add components to tibble (DRY principle)
.add_components_to_tibble <- function(tbl, data_matrix, n_comp, prefix) {
  if (is.null(data_matrix) || n_comp <= 0 || ncol(data_matrix) < 1) {
    return(tbl)
  }
  
  actual_n_comp <- min(n_comp, ncol(data_matrix))
  
  if (actual_n_comp == 1) {
    tbl[[paste0(prefix, "1")]] <- data_matrix[, 1]
  } else {
    comp_df <- tibble::as_tibble(data_matrix[, 1:actual_n_comp, drop = FALSE])
    names(comp_df) <- paste0(prefix, 1:actual_n_comp)
    tbl <- dplyr::bind_cols(tbl, comp_df)
  }
  
  tbl
}

# Tidy sample scores from OPLS-DA results
.tidy_oplsda_samples <- function(oplsda_res, mat, group_labels, group_col, predI, orthoI) {
  # Initialize samples tibble
  samples_tbl <- tibble::tibble(sample = rownames(mat))
  
  # Extract predictive and orthogonal scores
  pred_scores <- .safe_extract(ropls::getScoreMN(oplsda_res))
  ortho_scores <- .safe_extract(ropls::getScoreMN(oplsda_res, orthoL = TRUE))
  
  # Add predictive components
  samples_tbl <- .add_components_to_tibble(samples_tbl, pred_scores, predI, "pred")
  
  # Add orthogonal components if available
  if (!is.null(orthoI) && !is.na(orthoI) && orthoI > 0) {
    samples_tbl <- .add_components_to_tibble(samples_tbl, ortho_scores, orthoI, "ortho")
  }
  
  samples_tbl
}

# Tidy variable loadings from OPLS-DA results
.tidy_oplsda_variables <- function(oplsda_res, predI, orthoI) {
  # Extract predictive and orthogonal loadings
  pred_loadings <- .safe_extract(ropls::getLoadingMN(oplsda_res))
  ortho_loadings <- .safe_extract(ropls::getLoadingMN(oplsda_res, orthoL = TRUE))
  
  # Initialize variables tibble
  if (!is.null(pred_loadings)) {
    variables_tbl <- tibble::tibble(variable = rownames(pred_loadings))
  } else {
    variables_tbl <- tibble::tibble(variable = character(0))
  }
  
  # Add predictive components
  variables_tbl <- .add_components_to_tibble(variables_tbl, pred_loadings, predI, "pred")
  
  # Add orthogonal components if available
  if (!is.null(orthoI) && !is.na(orthoI) && orthoI > 0) {
    variables_tbl <- .add_components_to_tibble(variables_tbl, ortho_loadings, orthoI, "ortho")
  }
  
  variables_tbl
}

# Utility function to create component info (DRY principle)
.create_component_info <- function(summary_df, n_comp, prefix, comp_type, col_name) {
  components_list <- list()
  
  if (!is.null(summary_df) && n_comp > 0) {
    for (i in 1:n_comp) {
      # For ropls getSummaryDF, the variance info is in the "Total" row
      variance_value <- .safe_extract({
        if ("Total" %in% rownames(summary_df) && col_name %in% colnames(summary_df)) {
          as.numeric(summary_df["Total", col_name])
        } else {
          NA_real_
        }
      }, default = NA_real_)
      
      components_list[[paste0(prefix, i)]] <- list(
        component = paste0(prefix, i),
        type = comp_type,
        explained_variance = variance_value,
        cumulative_variance = variance_value
      )
    }
  }
  
  components_list
}

# Tidy component information from OPLS-DA results
.tidy_oplsda_components <- function(oplsda_res, predI, orthoI) {
  # Extract summary information
  summary_df <- .safe_extract(ropls::getSummaryDF(oplsda_res))
  
  # Create component tibble
  components_list <- list()
  
  # Add predictive components
  pred_components <- .create_component_info(
    summary_df, predI, "pred", "predictive", "R2X(cum)"
  )
  components_list <- c(components_list, pred_components)
  
  # Add orthogonal components if available
  if (!is.null(orthoI) && !is.na(orthoI) && orthoI > 0) {
    ortho_components <- .create_component_info(
      summary_df, orthoI, "ortho", "orthogonal", "R2X(cum)"
    )
    components_list <- c(components_list, ortho_components)
  }
  
  # Convert to tibble
  if (length(components_list) > 0) {
    components_tbl <- dplyr::bind_rows(components_list)
  } else {
    components_tbl <- tibble::tibble(
      component = character(0),
      type = character(0),
      explained_variance = numeric(0),
      cumulative_variance = numeric(0)
    )
  }
  
  components_tbl
} 