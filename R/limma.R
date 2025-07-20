#' Linear Models for Differential Expression Analysis using limma
#'
#' Perform differential expression analysis using linear models with empirical Bayes
#' moderation from the limma package. Supports both two-group and multi-group comparisons.
#'
#' @param exp A `glyexp_experiment` object containing expression data and sample information.
#' @param group_col A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param p_adj_method A character string specifying the method for multiple testing correction.
#'   Must be one of the methods supported by `stats::p.adjust()`. Default is "BH" (Benjamini-Hochberg).
#'   Set to NULL to skip p-value adjustment.
#' @param ref_group A character string specifying the reference group.
#'  If NULL (default), the first level of the group factor is used as the reference.
#'  Only used for two-group comparisons.
#' @param contrasts A character vector specifying custom contrasts for multi-group comparisons.
#'   If NULL (default), all pairwise comparisons are automatically generated.
#'   Supports two formats: "group1-group2" or "group1_vs_group2".
#'   Use the second format if group names contain hyphens.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw limma fit objects as a list.
#' @param ... Additional arguments passed to `limma::lmFit()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. The analysis uses linear models with empirical Bayes moderation
#' to improve statistical power, especially for small sample sizes.
#'
#' For two-group comparisons, a simple contrast is performed between the groups.
#' For multi-group comparisons (3+ groups), all pairwise comparisons are automatically
#' generated and performed using contrast matrices unless custom contrasts are specified.
#'
#' When specifying custom contrasts, use either "group1-group2" or "group1_vs_group2" format.
#' If group names contain hyphens and you use the first format, an error will be raised
#' suggesting to use the second format.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `limma` for linear model fitting and empirical Bayes moderation
#'
#' @returns A tibble with limma results including log2 fold change (log2fc).
#'   For multi-group comparisons, includes a "contrast" column indicating the comparison.
#'   If `return_raw` is TRUE, returns raw limma fit objects as a list.
#' @seealso [limma::lmFit()], [limma::eBayes()], [limma::makeContrasts()]
#' @export
gly_limma <- function(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  ref_group = NULL,
  contrasts = NULL,
  add_info = TRUE,
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods, null.ok = TRUE)
  checkmate::assert_character(contrasts, null.ok = TRUE)
  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

  # Check package availability
  .check_pkg_available("limma")

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)
  sample_info <- glyexp::get_sample_info(exp)

  # Extract and validate groups (minimum 2 groups, no maximum limit)
  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = NULL,
    method = "limma"
  )
  groups <- group_info$groups
  n_groups <- length(levels(groups))

  # Validate ref_group parameter
  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Perform limma analysis based on number of groups
  if (n_groups == 2) {
    # Two-group comparison
    result <- .gly_limma_2groups(expr_mat, groups, p_adj_method, ref_group, return_raw, ...)
  } else {
    # Multi-group comparison
    result <- .gly_limma_multigroups(expr_mat, groups, p_adj_method, contrasts, return_raw, ...)
  }

  # Return raw results if requested
  if (return_raw) {
    return(result)
  }

  # Process results with add_info logic
  result <- .process_results_add_info(result, exp, add_info)

  # Add S3 class
  structure(result, class = c("glystats_dea_res_limma", "glystats_dea_res", "glystats_res", class(result)))
}

# Internal helper functions for limma analysis

# Limma-specific 2-group analysis function
.gly_limma_2groups <- function(expr_mat, groups, p_adj_method, ref_group, return_raw = FALSE, ...) {
  # Reorder groups if ref_group is specified
  if (!is.null(ref_group)) {
    groups <- .reorder_groups_for_ref(groups, ref_group)
  }
  cli::cli_alert_info("Reference group: {.val {levels(groups)[1]}}")

  # Prepare data for limma
  log_expr_mat <- log2(expr_mat + 1)

  # Create design matrix
  design <- stats::model.matrix(~ groups)
  colnames(design) <- c("Intercept", stringr::str_c(levels(groups)[2], "_vs_", levels(groups)[1]))

  # Fit linear model
  fit <- limma::lmFit(log_expr_mat, design, ...)

  # Apply empirical Bayes moderation
  fit <- limma::eBayes(fit)

  if (return_raw) {
    return(fit)
  }

  # Extract results and convert to tibble
  .gly_limma_tibblify(fit, p_adj_method, expr_mat, groups)
}

# Convert limma fit object to tibble
.gly_limma_tibblify <- function(fit, p_adj_method, expr_mat, groups) {
  # Use topTable to extract results - this is the standard limma approach
  # The coefficient of interest is the second column (group comparison)
  coef_idx <- 2

  # Extract all results using topTable
  top_table <- limma::topTable(
    fit,
    coef = coef_idx,
    number = Inf,  # Get all genes
    adjust.method = if (is.null(p_adj_method)) "none" else p_adj_method,
    sort.by = "none"  # Keep original order
  )

  # Convert to tibble and rename columns to match glystats conventions
  result_tbl <- tibble::as_tibble(top_table, rownames = "variable") %>%
    dplyr::rename(
      log2fc = "logFC",
      p = "P.Value",
      p_adj = "adj.P.Val",
      t = "t",
      b = "B"
    )

  # Remove p_adj column if p_adj_method was NULL
  if (is.null(p_adj_method)) {
    result_tbl <- dplyr::select(result_tbl, -"p_adj")
  }

  result_tbl
}

# Helper function to reorder groups so that ref_group becomes the first level
.reorder_groups_for_ref <- function(groups, ref_group) {
  current_levels <- levels(groups)

  # Move ref_group to the first position
  new_levels <- c(ref_group, setdiff(current_levels, ref_group))

  # Reorder the factor levels
  factor(groups, levels = new_levels)
}

# Limma-specific multi-group analysis function
.gly_limma_multigroups <- function(expr_mat, groups, p_adj_method, contrasts = NULL, return_raw = FALSE, ...) {
  group_levels <- levels(groups)
  n_groups <- length(group_levels)

  cli::cli_alert_info("Performing multi-group limma analysis with {n_groups} groups")
  cli::cli_alert_info("Groups: {.val {group_levels}}")

  # Prepare data for limma
  log_expr_mat <- log2(expr_mat + 1)

  # Create design matrix without intercept (means model)
  design <- stats::model.matrix(~ 0 + groups)
  # Use make.names to ensure valid R names for limma
  valid_names <- make.names(group_levels)
  colnames(design) <- valid_names

  # Create mapping from original names to valid names
  name_mapping <- rlang::set_names(valid_names, group_levels)

  # Generate contrasts based on user input or default pairwise
  if (is.null(contrasts)) {
    # Generate all pairwise contrasts
    contrast_pairs <- .generate_pairwise_contrasts(group_levels)
  } else {
    # Parse user-specified contrasts
    contrast_pairs <- .parse_custom_contrasts(contrasts, group_levels)
  }

  # Create contrast matrix using limma::makeContrasts
  # Use valid names for limma contrast matrix
  contrast_strings <- purrr::map_chr(contrast_pairs, ~ {
    valid_name1 <- name_mapping[.x[1]]
    valid_name2 <- name_mapping[.x[2]]
    stringr::str_c(valid_name1, "-", valid_name2)
  })

  # Build the makeContrasts call
  contrast_matrix <- do.call(limma::makeContrasts,
                             c(as.list(contrast_strings),
                               list(levels = colnames(design))))

  # Fit linear model
  fit <- limma::lmFit(log_expr_mat, design, ...)

  # Apply contrasts
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)

  # Apply empirical Bayes moderation
  fit2 <- limma::eBayes(fit2)

  if (return_raw) {
    return(fit2)
  }

  # Extract results and convert to tibble
  .gly_limma_multi_tibblify(fit2, p_adj_method, contrast_pairs)
}

# Helper function to generate all pairwise contrasts
.generate_pairwise_contrasts <- function(group_levels) {
  n_groups <- length(group_levels)
  contrast_pairs <- list()

  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      contrast_pairs <- append(contrast_pairs, list(c(group_levels[i], group_levels[j])))
    }
  }

  return(contrast_pairs)
}

# Convert multi-group limma fit object to tibble
.gly_limma_multi_tibblify <- function(fit, p_adj_method, contrast_pairs) {
  # Extract single contrast result as tibble
  .extract_contrast_result <- function(i) {
    # Extract results for this contrast using topTable
    top_table <- limma::topTable(
      fit,
      coef = i,
      number = Inf,  # Get all genes
      adjust.method = if (is.null(p_adj_method)) "none" else p_adj_method,
      sort.by = "none"  # Keep original order
    )

    # Convert to tibble and rename columns
    result_tbl <- tibble::as_tibble(top_table, rownames = "variable") %>%
      dplyr::rename(
        log2fc = "logFC",
        p = "P.Value",
        p_adj = "adj.P.Val",
        t = "t",
        b = "B"
      ) %>%
      dplyr::mutate(
        contrast = stringr::str_c(contrast_pairs[[i]][1], "_vs_", contrast_pairs[[i]][2])
      )

    # Remove p_adj column if p_adj_method was NULL
    if (is.null(p_adj_method)) {
      result_tbl <- dplyr::select(result_tbl, -"p_adj")
    }

    result_tbl
  }

  # Use purrr::map to extract all contrasts efficiently and combine
  seq_along(contrast_pairs) %>%
    purrr::map(.extract_contrast_result) %>%
    dplyr::bind_rows()
}

# Helper function to parse custom contrasts
.parse_custom_contrasts <- function(contrasts, group_levels) {
  # Helper function to parse a single contrast
  .parse_single_contrast <- function(contrast) {
    # Try to parse using both formats
    if (stringr::str_detect(contrast, "_vs_")) {
      # Format: "group1_vs_group2"
      parts <- stringr::str_split(contrast, "_vs_")[[1]]
      if (length(parts) != 2) {
        cli::cli_abort("Invalid contrast format: {.val {contrast}}. Expected format: {.val group1_vs_group2}")
      }
      group1 <- parts[1]
      group2 <- parts[2]
    } else if (stringr::str_detect(contrast, "-")) {
      # Format: "group1-group2"
      # Check if any group names contain hyphens
      if (any(stringr::str_detect(group_levels, "-"))) {
        cli::cli_abort(c(
          "Group names contain hyphens, but contrast uses hyphen format: {.val {contrast}}",
          "i" = "Use the format {.val group1_vs_group2} instead when group names contain hyphens"
        ))
      }
      parts <- stringr::str_split(contrast, "-")[[1]]
      if (length(parts) != 2) {
        cli::cli_abort("Invalid contrast format: {.val {contrast}}. Expected format: {.val group1-group2}")
      }
      group1 <- parts[1]
      group2 <- parts[2]
    } else {
      cli::cli_abort(c(
        "Invalid contrast format: {.val {contrast}}",
        "i" = "Use either {.val group1-group2} or {.val group1_vs_group2} format"
      ))
    }

    # Validate that both groups exist
    if (!group1 %in% group_levels) {
      cli::cli_abort("Group {.val {group1}} not found in data. Available groups: {.val {group_levels}}")
    }
    if (!group2 %in% group_levels) {
      cli::cli_abort("Group {.val {group2}} not found in data. Available groups: {.val {group_levels}}")
    }

    # Return as character vector
    c(group1, group2)
  }

  # Use purrr::map to parse all contrasts
  purrr::map(contrasts, .parse_single_contrast)
}
