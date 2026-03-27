#' Linear Models for Differential Expression Analysis using limma
#'
#' Perform differential expression analysis using linear models with empirical Bayes
#' moderation from the limma package. Supports both two-group and multi-group comparisons.
#'
#' @param exp A `glyexp_experiment` object containing expression data and sample information.
#' @param expr_mat (Only for [gly_limma_()]) A numeric matrix with variables as rows and samples as columns.
#' @param groups (Only for [gly_limma_()]) A factor or character vector specifying group membership for each sample.
#'   Must have at least 2 levels. Character vectors will be automatically converted to factors.
#'   If `contrasts` is not provided,
#'   the levels coming first in the factor will be used as the reference group.
#' @param group_col (Only for [gly_limma()]) A character string specifying the column name in sample information
#'   that contains group labels. Default is "group".
#' @param covariate_cols (Only for [gly_limma()]) A character vector specifying column names in sample information
#'   to include as covariates in the limma model. Default is NULL.
#' @param subject_col (Only for [gly_limma()]) A character string specifying the column name in sample information
#'   that contains subject identifiers for paired comparisons. Default is NULL.
#' @param covariates (Only for [gly_limma_()]) A data frame, matrix, or vector of sample-level covariates.
#'   Must have the same number of rows as `expr_mat` has columns. If row names are provided and
#'   match `colnames(expr_mat)`, they will be aligned automatically. Default is NULL.
#' @param subjects (Only for [gly_limma_()]) A vector or factor of subject identifiers for paired comparisons.
#'   Must have length equal to `ncol(expr_mat)`. If names or row names are provided and match
#'   `colnames(expr_mat)`, they will be aligned automatically. Default is NULL.
#' @param p_adj_method A character string specifying the method for multiple testing correction.
#'   Must be one of the methods supported by `stats::p.adjust()`. Default is "BH" (Benjamini-Hochberg).
#'   Set to NULL to skip p-value adjustment.
#' @param ref_group A character string specifying the reference group.
#'  If NULL (default), the first level of the group factor is used as the reference.
#'  Only used for two-group comparisons.
#' @param contrasts A character vector specifying custom contrasts.
#'   If NULL (default), all pairwise comparisons are automatically generated,
#'   and the levels coming first in the factor will be used as the reference group.
#'   Supports two formats: "group1-group2" or "group1_vs_group2".
#'   Use the second format if group names contain hyphens.
#'   "group1" will be used as the reference group.
#' @param add_info A logical value. If TRUE (default), variable information from the experiment
#'  will be added to the result tibble. If FALSE, only the statistical results are returned.
#'  Only applicable to `gly_limma()`.
#' @param ... Additional arguments passed to `limma::lmFit()`.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' statistical testing. The analysis uses linear models with empirical Bayes moderation
#' to improve statistical power, especially for small sample sizes.
#'
#' `gly_limma()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_limma_()` is the underlying API that works with matrices and factor vectors directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' For two or more groups, all pairwise comparisons are automatically generated
#' (one contrast for two groups) and performed using contrast matrices unless
#' custom contrasts are specified.
#'
#' If covariates are provided, they are added as additional terms in the design
#' matrix and are not part of the group contrasts.
#'
#' If subjects are provided, they are included as a blocking factor in the design
#' matrix using the formula `~ 0 + group + subject` (plus any covariates).
#'
#' When specifying custom contrasts, use either "group1-group2" or "group1_vs_group2" format.
#' If group names contain hyphens and you use the first format, an error will be raised
#' suggesting to use the second format.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `limma` for linear model fitting and empirical Bayes moderation
#'
#' @returns A list with two elements:
#'  - `tidy_result`: A tibble with limma results containing the following columns:
#'    - `variable`: Variable name
#'    - `log2fc`: Log2 fold change
#'    - `AveExpr`: Average expression level
#'    - `t`: t-statistic
#'    - `p_val`: Raw p-value
#'    - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)
#'    - `b`: B-statistic (log-odds of differential expression)
#'    - `ref_group`: Reference group
#'    - `test_group`: Test group
#'
#'  - `raw_result`: The raw limma fit object(s).
#' @seealso [limma::lmFit()], [limma::eBayes()], [limma::makeContrasts()]
#' @export
gly_limma <- function(
  exp,
  group_col = "group",
  covariate_cols = NULL,
  subject_col = NULL,
  p_adj_method = "BH",
  ref_group = NULL,
  contrasts = NULL,
  add_info = TRUE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(group_col)
  if (length(covariate_cols) == 0) {
    covariate_cols <- NULL
  }
  if (length(subject_col) == 0) {
    subject_col <- NULL
  }
  checkmate::assert_character(covariate_cols, null.ok = TRUE)
  checkmate::assert_character(subject_col, null.ok = TRUE)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_character(contrasts, null.ok = TRUE)
  checkmate::assert_logical(add_info, len = 1)

  # Check package availability
  rlang::check_installed("limma")

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

  covariates <- .extract_covariates_from_sample_info(
    sample_info,
    covariate_cols,
    group_col
  )
  subjects <- .extract_subjects_from_sample_info(
    sample_info,
    subject_col,
    group_col,
    covariate_cols
  )

  # Validate ref_group parameter
  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  # Call the underlying API
  result <- gly_limma_(
    expr_mat,
    groups,
    p_adj_method,
    ref_group,
    contrasts,
    covariates,
    subjects = subjects,
    ...
  )

  # Process results with add_info logic
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  result
}

#' @rdname gly_limma
#' @export
gly_limma_ <- function(
  expr_mat,
  groups,
  p_adj_method = "BH",
  ref_group = NULL,
  contrasts = NULL,
  covariates = NULL,
  subjects = NULL,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  groups <- .convert_groups_to_factor(groups)
  checkmate::assert_factor(groups, len = ncol(expr_mat))
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_character(contrasts, null.ok = TRUE)

  # Check package availability
  rlang::check_installed("limma")

  # Validate groups
  n_groups <- length(levels(groups))
  if (n_groups < 2) {
    cli::cli_abort("groups must have at least 2 levels for limma analysis")
  }

  # Validate ref_group parameter
  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)

  covariates <- .normalize_covariates(
    covariates,
    ncol(expr_mat),
    colnames(expr_mat)
  )
  subjects <- .normalize_subjects(subjects, ncol(expr_mat), colnames(expr_mat))

  # Align reference group for 2-group comparisons
  if (n_groups == 2 && !is.null(ref_group)) {
    groups <- .reorder_groups_for_ref(groups, ref_group)
  }

  # Use the multi-group workflow for all cases
  result <- .gly_limma_multigroups(
    expr_mat,
    groups,
    p_adj_method,
    contrasts,
    covariates,
    subjects,
    ...
  )

  result
}

# Internal helper functions for limma analysis

# Helper function to reorder groups so that ref_group becomes the first level
.reorder_groups_for_ref <- function(groups, ref_group) {
  current_levels <- levels(groups)

  # Move ref_group to the first position
  new_levels <- c(ref_group, setdiff(current_levels, ref_group))

  # Reorder the factor levels
  factor(groups, levels = new_levels)
}

# Extract covariates from sample_info for gly_limma
.extract_covariates_from_sample_info <- function(
  sample_info,
  covariate_cols,
  group_col
) {
  if (is.null(covariate_cols) || length(covariate_cols) == 0) {
    return(NULL)
  }

  if (anyDuplicated(covariate_cols) > 0) {
    cli::cli_abort("covariate_cols must be unique.")
  }
  if (group_col %in% covariate_cols) {
    cli::cli_abort(
      "covariate_cols cannot include the group column {.field {group_col}}."
    )
  }

  missing_cols <- setdiff(covariate_cols, colnames(sample_info))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "covariate_cols not found in sample information: {.field {missing_cols}}.",
      "i" = "Available columns: {.field {colnames(sample_info)}}"
    ))
  }

  sample_info[, covariate_cols, drop = FALSE]
}

# Extract subjects from sample_info for gly_limma
.extract_subjects_from_sample_info <- function(
  sample_info,
  subject_col,
  group_col,
  covariate_cols
) {
  if (is.null(subject_col) || length(subject_col) == 0) {
    return(NULL)
  }

  if (anyDuplicated(subject_col) > 0) {
    cli::cli_abort("subject_col must be unique.")
  }
  if (length(subject_col) != 1) {
    cli::cli_abort("subject_col must contain exactly 1 column name.")
  }
  if (group_col %in% subject_col) {
    cli::cli_abort(
      "subject_col cannot include the group column {.field {group_col}}."
    )
  }
  if (!is.null(covariate_cols) && any(subject_col %in% covariate_cols)) {
    cli::cli_abort("subject_col cannot overlap with covariate_cols.")
  }

  missing_cols <- setdiff(subject_col, colnames(sample_info))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "subject_col not found in sample information: {.field {missing_cols}}.",
      "i" = "Available columns: {.field {colnames(sample_info)}}"
    ))
  }

  subjects <- sample_info[[subject_col]]
  if (!is.factor(subjects)) {
    subjects <- factor(subjects)
  }

  subjects
}

# Normalize covariates for gly_limma_
.normalize_covariates <- function(covariates, n_samples, sample_names = NULL) {
  if (is.null(covariates)) {
    return(NULL)
  }

  if (is.vector(covariates) && !is.list(covariates)) {
    covariates <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
  } else if (is.matrix(covariates)) {
    covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
  } else if (is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
  } else {
    cli::cli_abort("covariates must be a data.frame, matrix, or vector.")
  }

  if (ncol(covariates) == 0) {
    return(NULL)
  }

  if (nrow(covariates) != n_samples) {
    cli::cli_abort(
      "covariates must have {.val {n_samples}} rows to match expr_mat columns."
    )
  }

  if (!is.null(sample_names)) {
    covariate_rows <- rownames(covariates)
    if (!is.null(covariate_rows)) {
      if (setequal(covariate_rows, sample_names)) {
        covariates <- covariates[sample_names, , drop = FALSE]
      } else if (!identical(covariate_rows, as.character(seq_len(n_samples)))) {
        cli::cli_abort("covariates row names must match expr_mat column names.")
      }
    }
  }

  if ("groups" %in% colnames(covariates)) {
    cli::cli_abort(
      "covariates cannot include a column named {.field groups}; rename it."
    )
  }
  if ("subject" %in% colnames(covariates)) {
    cli::cli_abort(
      "covariates cannot include a column named {.field subject}; rename it."
    )
  }

  covariates[] <- lapply(covariates, function(col) {
    if (is.character(col)) {
      factor(col)
    } else {
      col
    }
  })

  covariates
}

# Normalize subjects for gly_limma_
.normalize_subjects <- function(subjects, n_samples, sample_names = NULL) {
  if (is.null(subjects)) {
    return(NULL)
  }

  if (is.factor(subjects) || (is.vector(subjects) && !is.list(subjects))) {
    # keep as is
  } else if (is.matrix(subjects) || is.data.frame(subjects)) {
    subjects <- as.data.frame(subjects, stringsAsFactors = FALSE)
    if (ncol(subjects) == 0) {
      return(NULL)
    }
    if (ncol(subjects) != 1) {
      cli::cli_abort("subjects must have exactly 1 column.")
    }
    if (!is.null(sample_names)) {
      subject_rows <- rownames(subjects)
      if (!is.null(subject_rows)) {
        if (setequal(subject_rows, sample_names)) {
          subjects <- subjects[sample_names, , drop = FALSE]
        } else if (!identical(subject_rows, as.character(seq_len(n_samples)))) {
          cli::cli_abort("subjects row names must match expr_mat column names.")
        }
      }
    }
    subjects <- subjects[[1]]
  } else {
    cli::cli_abort("subjects must be a vector, matrix, or data.frame.")
  }

  if (length(subjects) != n_samples) {
    cli::cli_abort(
      "subjects must have {.val {n_samples}} values to match expr_mat columns."
    )
  }

  if (!is.null(sample_names)) {
    subject_names <- names(subjects)
    if (!is.null(subject_names)) {
      if (setequal(subject_names, sample_names)) {
        subjects <- subjects[sample_names]
      } else if (!identical(subject_names, as.character(seq_len(n_samples)))) {
        cli::cli_abort("subjects names must match expr_mat column names.")
      }
    }
  }

  if (!is.factor(subjects)) {
    subjects <- factor(subjects)
  }

  subjects
}

# Limma-specific multi-group analysis function
.gly_limma_multigroups <- function(
  expr_mat,
  groups,
  p_adj_method,
  contrasts = NULL,
  covariates = NULL,
  subjects = NULL,
  ...
) {
  group_levels <- levels(groups)
  n_groups <- length(group_levels)

  # Prepare data for limma
  log_expr_mat <- log2(expr_mat + 1e-6)
  row_has_data <- rowSums(is.finite(log_expr_mat)) > 0
  missing_vars <- character(0)
  if (!all(row_has_data)) {
    missing_vars <- rownames(log_expr_mat)[!row_has_data]
    log_expr_mat <- log_expr_mat[row_has_data, , drop = FALSE]
  }
  if (nrow(log_expr_mat) == 0) {
    cli::cli_abort(
      "No finite expression values remain after log2 transformation."
    )
  }

  # Create design matrix without intercept (means model)
  if (is.null(covariates) && is.null(subjects)) {
    design_data <- data.frame(groups = groups, check.names = FALSE)
    design_formula <- stats::as.formula("~ 0 + groups")
  } else {
    design_data <- data.frame(groups = groups, check.names = FALSE)
    if (!is.null(subjects)) {
      design_data$subject <- subjects
    }
    if (!is.null(covariates)) {
      design_data <- data.frame(design_data, covariates, check.names = FALSE)
    }
    design_formula <- stats::as.formula("~ 0 + groups + .")
  }
  design <- stats::model.matrix(design_formula, data = design_data)
  # Use make.names to ensure valid R names for limma group terms
  valid_names <- make.names(group_levels)
  design_terms <- stats::terms(design_formula, data = design_data)
  term_labels <- attr(design_terms, "term.labels")
  group_term_index <- which(term_labels == "groups")
  if (length(group_term_index) != 1) {
    cli::cli_abort("Failed to locate group terms in the limma design matrix.")
  }
  group_cols <- which(attr(design, "assign") == group_term_index)
  if (length(group_cols) != length(valid_names)) {
    cli::cli_abort(
      "Group columns in the design matrix do not match group levels."
    )
  }
  colnames(design)[group_cols] <- valid_names

  # Create mapping from original names to valid names
  name_mapping <- rlang::set_names(valid_names, group_levels)

  # Generate contrasts based on user input or default pairwise
  if (is.null(contrasts)) {
    # Generate all pairwise contrasts
    contrast_pairs <- .make_comparisons(group_levels, reverse = TRUE)
  } else {
    # Parse user-specified contrasts
    contrast_pairs <- .parse_custom_contrasts(contrasts, group_levels)
  }

  # Create contrast matrix using limma::makeContrasts
  # Use valid names for limma contrast matrix
  contrast_strings <- purrr::map_chr(
    contrast_pairs,
    ~ {
      valid_name1 <- name_mapping[.x[1]]
      valid_name2 <- name_mapping[.x[2]]
      stringr::str_c(valid_name1, "-", valid_name2)
    }
  )

  # Build the makeContrasts call
  contrast_matrix <- do.call(
    limma::makeContrasts,
    c(as.list(contrast_strings), list(levels = colnames(design)))
  )

  # Fit linear model
  fit <- limma::lmFit(log_expr_mat, design, ...)
  if (is.null(rownames(fit$coefficients)) && !is.null(rownames(log_expr_mat))) {
    rownames(fit$coefficients) <- rownames(log_expr_mat)
    rownames(fit$stdev.unscaled) <- rownames(log_expr_mat)
    names(fit$sigma) <- rownames(log_expr_mat)
    names(fit$Amean) <- rownames(log_expr_mat)
  }

  # Apply contrasts
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)

  # Apply empirical Bayes moderation
  fit2 <- tryCatch(
    limma::eBayes(fit2, trend = TRUE),
    error = function(err) {
      msg <- conditionMessage(err)
      if (grepl("covariate", msg, ignore.case = TRUE)) {
        cli::cli_warn(
          "Limma trend fitting failed due to covariate issues; re-running eBayes with trend = FALSE."
        )
        limma::eBayes(fit2, trend = FALSE)
      } else {
        stop(err)
      }
    }
  )

  # Extract results and convert to tibble
  tidy_result <- .gly_limma_multi_tibblify(fit2, p_adj_method, contrast_pairs)
  if (length(missing_vars) > 0) {
    missing_tbl <- purrr::map_dfr(contrast_pairs, function(pair) {
      base_tbl <- tibble::tibble(
        variable = missing_vars,
        log2fc = NA_real_,
        AveExpr = NA_real_,
        t = NA_real_,
        p_val = NA_real_,
        b = NA_real_,
        ref_group = pair[2],
        test_group = pair[1]
      )
      if (!is.null(p_adj_method)) {
        base_tbl$p_adj <- NA_real_
      }
      base_tbl
    })
    tidy_result <- dplyr::bind_rows(tidy_result, missing_tbl)
  }

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = fit2
    ),
    class = c("glystats_limma_res", "glystats_res")
  )
}

# Convert multi-group limma fit object to tibble
.gly_limma_multi_tibblify <- function(fit, p_adj_method, contrast_pairs) {
  # Extract single contrast result as tibble
  .extract_contrast_result <- function(i) {
    # Extract results for this contrast using topTable
    top_table <- limma::topTable(
      fit,
      coef = i,
      number = Inf, # Get all genes
      adjust.method = if (is.null(p_adj_method)) "none" else p_adj_method,
      sort.by = "none" # Keep original order
    )

    # Convert to tibble and rename columns
    result_tbl <- tibble::as_tibble(top_table, rownames = "variable") %>%
      dplyr::rename(
        log2fc = "logFC",
        p_val = "P.Value",
        p_adj = "adj.P.Val",
        t = "t",
        b = "B"
      ) %>%
      dplyr::mutate(
        ref_group = contrast_pairs[[i]][2],
        test_group = contrast_pairs[[i]][1],
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
        cli::cli_abort(
          "Invalid contrast format: {.val {contrast}}. Expected format: {.val group1_vs_group2}"
        )
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
        cli::cli_abort(
          "Invalid contrast format: {.val {contrast}}. Expected format: {.val group1-group2}"
        )
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
      cli::cli_abort(
        "Group {.val {group1}} not found in data. Available groups: {.val {group_levels}}"
      )
    }
    if (!group2 %in% group_levels) {
      cli::cli_abort(
        "Group {.val {group2}} not found in data. Available groups: {.val {group_levels}}"
      )
    }

    # Return as character vector
    c(group2, group1)
  }

  # Use purrr::map to parse all contrasts
  purrr::map(contrasts, .parse_single_contrast)
}
