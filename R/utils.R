# Data-container helpers -------------------------------------------------------

#' Check a glystats data-container input
#'
#' @param exp A `glyexp::GlycomicSE`, `glyexp::GlycoproteomicSE`, or another
#'   `SummarizedExperiment` object.
#'
#' @returns `NULL` invisibly, or an error if `exp` is unsupported.
#' @noRd
.assert_data_container <- function(exp) {
  if (
    inherits(exp, "glyexp_experiment") ||
      methods::is(exp, "SummarizedExperiment")
  ) {
    return(invisible(NULL))
  }

  cli::cli_abort(
    "{.arg exp} must be a {.cls glyexp_experiment} or {.cls SummarizedExperiment} object."
  )
}

#' Extract the expression matrix from a glystats data container
#'
#' @inheritParams .assert_data_container
#'
#' @returns A numeric matrix with variables in rows and samples in columns.
#' @noRd
.get_expr_mat <- function(exp) {
  if (methods::is(exp, "SummarizedExperiment")) {
    expr_mat <- as.matrix(SummarizedExperiment::assay(exp, 1))
    dimnames(expr_mat) <- list(
      .se_axis_ids(exp, "variable"),
      .se_axis_ids(exp, "sample")
    )
    return(expr_mat)
  }
  glyexp::get_expr_mat(exp)
}

#' Extract sample information from a glystats data container
#'
#' @inheritParams .assert_data_container
#'
#' @returns A tibble of sample metadata in input-column order.
#' @noRd
.get_sample_info <- function(exp) {
  if (methods::is(exp, "SummarizedExperiment")) {
    return(tibble::as_tibble(SummarizedExperiment::colData(exp)))
  }
  glyexp::get_sample_info(exp)
}

#' Extract variable information from a glystats data container
#'
#' @inheritParams .assert_data_container
#'
#' @returns A tibble of variable metadata in input-row order.
#' @noRd
.get_var_info <- function(exp) {
  if (methods::is(exp, "SummarizedExperiment")) {
    return(tibble::as_tibble(SummarizedExperiment::rowData(exp)))
  }
  glyexp::get_var_info(exp)
}

# Resolve IDs from SE axes, never from ordinary annotation columns.
# Temporary IDs are local to the supplied object's current row/column order.
.se_axis_map <- function(exp, axis) {
  original <- if (axis == "variable") rownames(exp) else colnames(exp)
  n <- if (axis == "variable") nrow(exp) else ncol(exp)
  generated <- is.null(original) ||
    anyNA(original) ||
    any(!nzchar(original)) ||
    anyDuplicated(original) > 0L
  ids <- if (n == 0L) {
    character()
  } else if (generated) {
    paste0(".glystats_", axis, "_", seq_len(n))
  } else {
    original
  }
  list(
    data = tibble::tibble(
      position = seq_len(n),
      id = ids,
      name = if (is.null(original)) rep(NA_character_, n) else original
    ),
    generated = generated,
    names_missing = is.null(original)
  )
}

.se_axis_ids <- function(exp, axis) {
  .se_axis_map(exp, axis)$data$id
}

# Retain the ID-to-position mapping independently of add_info. Generated IDs
# are local to this input, not persistent identifiers across arbitrary subsets.
.record_se_identity <- function(result, exp) {
  if (!methods::is(exp, "SummarizedExperiment")) {
    return(result)
  }
  variable <- .se_axis_map(exp, "variable")
  sample <- .se_axis_map(exp, "sample")
  result$identifiers <- list(
    variable = variable,
    sample = sample,
    columns = c(variable = ".variable", sample = ".sample")
  )
  # These analyses do not attach metadata through the common finalizer.
  if (inherits(result, c("glystats_cor_res", "glystats_set_test_res"))) {
    result$tidy_result <- .namespace_se_result_ids(result$tidy_result)
  }
  filterable <- inherits(
    result,
    c(
      "glystats_ttest_res",
      "glystats_wilcox_res",
      "glystats_anova_res",
      "glystats_ancova_res",
      "glystats_kruskal_res",
      "glystats_limma_res",
      "glystats_linear_model_res"
    )
  )
  if (variable$generated && filterable) {
    result$identifiers$input_hash <- rlang::hash(exp)
  }
  result
}

# Position-keyed results may only be filtered against their unmodified input.
# Names alone cannot detect a permutation of duplicate/unnamed rows.
.validate_filter_identity <- function(exp, res) {
  recorded <- res$identifiers
  if (!methods::is(exp, "SummarizedExperiment")) {
    if (isTRUE(recorded$variable$generated)) {
      cli::cli_abort(
        "Position-keyed results require the same unmodified experiment used for analysis."
      )
    }
    return(invisible(NULL))
  }
  variable <- .se_axis_map(exp, "variable")
  if (variable$generated || isTRUE(recorded$variable$generated)) {
    if (
      is.null(recorded$input_hash) ||
        !identical(recorded$input_hash, rlang::hash(exp))
    ) {
      cli::cli_abort(c(
        "Position-keyed results require the same unmodified experiment used for analysis.",
        "i" = "Re-run the analysis after changing or reordering the experiment."
      ))
    }
  }
  invisible(NULL)
}

# Attach annotations by the position of each result ID on the input SE axis.
# Matching against metadata values would reintroduce the legacy ID contract.
.add_se_info <- function(tbl, exp, axis) {
  key <- paste0(".", axis)
  if (!key %in% names(tbl)) {
    return(tbl)
  }
  info <- if (axis == "variable") .get_var_info(exp) else .get_sample_info(exp)
  if (ncol(info) == 0L) {
    return(tbl)
  }
  ids <- .se_axis_ids(exp, axis)
  indices <- match(tbl[[key]], ids)
  if (anyNA(indices)) {
    cli::cli_abort(
      "Result {.field {axis}} IDs do not match the input experiment."
    )
  }
  # Retain the historical add_info ordering: input-axis order, preserving
  # the relative order of repeated result rows (e.g. contrasts/components).
  ordering <- order(indices)
  tbl <- tbl[ordering, , drop = FALSE]
  indices <- indices[ordering]
  # Preserve all ordinary metadata columns without overwriting result fields.
  # Reserve existing names first so an original *_metadata column is kept.
  used <- union(names(tbl), names(info))
  info_names <- names(info)
  for (i in seq_along(info_names)) {
    if (info_names[i] %in% names(tbl)) {
      candidate <- paste0(info_names[i], "_metadata")
      while (candidate %in% used) {
        candidate <- paste0(candidate, "_metadata")
      }
      info_names[i] <- candidate
      used <- c(used, candidate)
    }
  }
  names(info) <- info_names
  dplyr::bind_cols(
    tbl[, key, drop = FALSE],
    info[indices, , drop = FALSE],
    tbl[, setdiff(names(tbl), key), drop = FALSE]
  )
}

#' Extract metadata from a glystats data container
#'
#' @inheritParams .assert_data_container
#'
#' @returns A metadata list.
#' @noRd
.get_meta_data <- function(exp) {
  if (methods::is(exp, "SummarizedExperiment")) {
    return(S4Vectors::metadata(exp))
  }
  glyexp::get_meta_data(exp)
}

#' Extract the experiment type from a glystats data container
#'
#' @inheritParams .assert_data_container
#'
#' @returns A character scalar or `NULL`.
#' @noRd
.get_exp_type <- function(exp) {
  .get_meta_data(exp)$exp_type
}

#' Filter variables in a glystats data container
#'
#' @inheritParams .assert_data_container
#' @param variables A character vector of variable identifiers to retain.
#'
#' @returns The filtered container, preserving its input class.
#' @noRd
.filter_variables <- function(exp, variables) {
  if (methods::is(exp, "SummarizedExperiment")) {
    return(exp[.se_axis_ids(exp, "variable") %in% variables, , drop = FALSE])
  }
  glyexp::filter_var(exp, .data$variable %in% .env$variables)
}

# Identify variables without any observed expression values
.all_na_variables <- function(expr_mat) {
  rownames(expr_mat)[rowSums(!is.na(expr_mat)) == 0]
}

# Remove variables without any observed expression values
.remove_all_na_variables <- function(expr_mat) {
  expr_mat[rowSums(!is.na(expr_mat)) > 0, , drop = FALSE]
}

# Restore one-row-per-variable results for variables omitted while fitting
.restore_all_na_variable_rows <- function(
  tbl,
  missing_variables,
  variable_order
) {
  if (length(missing_variables) == 0) {
    return(tbl)
  }

  missing_tbl <- tbl[
    rep(NA_integer_, length(missing_variables)),
    ,
    drop = FALSE
  ]
  missing_tbl$variable <- missing_variables

  dplyr::bind_rows(tbl, missing_tbl) |>
    dplyr::arrange(match(.data$variable, variable_order))
}

# Restore long-format component results for variables omitted while fitting
.restore_all_na_variable_components <- function(
  tbl,
  missing_variables,
  variable_order,
  component_col
) {
  if (length(missing_variables) == 0) {
    return(tbl)
  }

  components <- unique(tbl[[component_col]])
  missing_tbl <- tidyr::expand_grid(
    variable = missing_variables,
    component = components
  )
  colnames(missing_tbl)[2] <- component_col
  missing_tbl$value <- NA_real_

  dplyr::bind_rows(tbl, missing_tbl) |>
    dplyr::arrange(
      match(.data$variable, variable_order),
      match(.data[[component_col]], components)
    )
}

# Group extraction, validation, and conversion helper functions ---------------

# Check if group column exists in sample information
.check_group_column_exists <- function(sample_info, group_col) {
  if (!group_col %in% colnames(sample_info)) {
    cli::cli_abort(c(
      "Column {.field {group_col}} not found in sample information.",
      "i" = "Available columns: {.field {colnames(sample_info)}}",
      "i" = "Did you mistype the column name or use column other than {.field {group_col}}?"
    ))
  }
}

# Extract and convert groups to factor
.extract_groups <- function(sample_info, group_col) {
  groups <- sample_info[[group_col]]
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }
  droplevels(groups)
}

# Helper function to generate validation error message
.generate_validation_error <- function(
  group_col,
  method,
  groups,
  message_type,
  count
) {
  if (is.null(method)) {
    cli::cli_abort(c(
      "{.field {group_col}} must have {message_type} {.val {count}} levels",
      "i" = "Current levels: {.val {levels(groups)}}"
    ))
  } else {
    cli::cli_abort(c(
      "{.field {group_col}} must have {message_type} {.val {count}} levels for {.val {method}}",
      "i" = "Current levels: {.val {levels(groups)}}"
    ))
  }
}

# Validate group count for specific analysis methods
.validate_group_count <- function(
  groups,
  group_col,
  min_count = NULL,
  max_count = NULL,
  method = NULL
) {
  # Skip validation if both limits are NULL
  if (is.null(min_count) && is.null(max_count)) {
    return()
  }

  n_groups <- length(levels(groups))

  # Handle exact count case (when min_count == max_count)
  if (!is.null(min_count) && !is.null(max_count) && min_count == max_count) {
    if (n_groups != min_count) {
      .generate_validation_error(
        group_col,
        method,
        groups,
        "exactly",
        min_count
      )
    }
    return() # Early return, no need to check min/max separately
  }

  # Check minimum count
  if (!is.null(min_count) && n_groups < min_count) {
    .generate_validation_error(group_col, method, groups, "at least", min_count)
  }

  # Check maximum count
  if (!is.null(max_count) && n_groups > max_count) {
    .generate_validation_error(group_col, method, groups, "at most", max_count)
  }
}

# Display group information based on group count
.display_group_info <- function(groups) {
  n_groups <- length(levels(groups))
  if (n_groups == 2) {
    .display_two_group_info(groups)
  } else {
    .display_multi_group_info(groups)
  }
}

# Display group information for two-group analysis
.display_two_group_info <- function(groups) {
  cli::cli_alert_info("Ref Group: {.val {levels(groups)[1]}}")
  cli::cli_alert_info("Test Group: {.val {levels(groups)[2]}}")
}

# Display group information for multi-group analysis
.display_multi_group_info <- function(groups) {
  n_groups <- length(levels(groups))
  cli::cli_alert_info("Number of groups: {.val {n_groups}}")
  cli::cli_alert_info("Groups: {.val {levels(groups)}}")
  cli::cli_alert_info(
    "Pairwise comparisons will be performed, with levels coming first as reference groups."
  )
}

# Extract and validate groups (comprehensive function)
.extract_and_validate_groups <- function(
  sample_info,
  group_col,
  min_count = NULL,
  max_count = NULL,
  method = NULL,
  show_info = TRUE
) {
  .check_group_column_exists(sample_info, group_col)
  groups <- .extract_groups(sample_info, group_col)
  .validate_group_count(groups, group_col, min_count, max_count, method)
  if (show_info) {
    .display_group_info(groups)
  }
  list(groups = groups, group_col = group_col)
}

# Add info helper functions ---------------

# Check if a tibble contains variable or sample columns
.has_variable_column <- function(tbl) {
  "variable" %in% colnames(tbl) || "column" %in% colnames(tbl)
}

.has_sample_column <- function(tbl) {
  "sample" %in% colnames(tbl) || "row" %in% colnames(tbl)
}

# Add variable information to a tibble
.add_variable_info <- function(tbl, exp) {
  if (methods::is(exp, "SummarizedExperiment")) {
    return(.add_se_info(tbl, exp, "variable"))
  }
  if (!.has_variable_column(tbl)) {
    return(tbl)
  }

  var_info <- .get_var_info(exp)
  # Remove the variable column from var_info to avoid duplication
  var_info_subset <- var_info[,
    !colnames(var_info) %in% "variable",
    drop = FALSE
  ]

  # Only join if there are columns other than variable
  if (ncol(var_info_subset) > 0) {
    # Join var_info to tbl, keeping all rows from tbl
    tbl <- dplyr::right_join(var_info, tbl, by = "variable")
  }

  tbl
}

# Add sample information to a tibble
.add_sample_info <- function(tbl, exp) {
  if (methods::is(exp, "SummarizedExperiment")) {
    return(.add_se_info(tbl, exp, "sample"))
  }
  if (!.has_sample_column(tbl)) {
    return(tbl)
  }

  sample_info <- .get_sample_info(exp)
  # Remove the sample column from sample_info to avoid duplication
  sample_info_subset <- sample_info[,
    !colnames(sample_info) %in% "sample",
    drop = FALSE
  ]

  # Only join if there are columns other than sample
  if (ncol(sample_info_subset) > 0) {
    # Join sample_info to tbl, keeping all rows from tbl
    tbl <- dplyr::right_join(sample_info, tbl, by = "sample")
  }

  tbl
}

# Process a single tibble with add_info logic
# Run once, before adding annotations, so ordinary metadata is never renamed.
.namespace_se_result_ids <- function(x) {
  if (tibble::is_tibble(x)) {
    identity <- c(
      "variable",
      "sample",
      "variable1",
      "variable2",
      "sample1",
      "sample2"
    )
    matched <- names(x) %in% identity
    names(x)[matched] <- paste0(".", names(x)[matched])
    return(x)
  }
  if (is.list(x)) {
    return(lapply(x, .namespace_se_result_ids))
  }
  x
}

.process_tibble_add_info <- function(tbl, exp, add_info) {
  if (methods::is(exp, "SummarizedExperiment")) {
    tbl <- .namespace_se_result_ids(tbl)
    if (add_info) {
      has_variable <- ".variable" %in% names(tbl)
      has_sample <- ".sample" %in% names(tbl)
      if (has_variable) {
        tbl <- .add_se_info(tbl, exp, "variable")
      }
      if (has_sample) tbl <- .add_se_info(tbl, exp, "sample")
    }
    return(tbl)
  }
  if (!add_info) {
    return(tbl)
  }
  has_variable <- .has_variable_column(tbl)
  has_sample <- .has_sample_column(tbl)
  if (has_variable) {
    tbl <- .add_variable_info(tbl, exp)
  }
  if (has_sample) {
    tbl <- .add_sample_info(tbl, exp)
  }
  tbl
}

# Process results with add_info logic (handles both single tibbles and lists)
.process_results_add_info <- function(results, exp, add_info) {
  if (!add_info && !methods::is(exp, "SummarizedExperiment")) {
    return(results)
  }

  old_class <- class(results)

  # If results is a tibble, process it directly
  if (tibble::is_tibble(results)) {
    new_results <- .process_tibble_add_info(results, exp, add_info)
  } else if (is.list(results)) {
    new_results <- purrr::map(results, function(item) {
      if (tibble::is_tibble(item)) {
        .process_tibble_add_info(item, exp, add_info)
      } else {
        item
      }
    })
  } else {
    cli::cli_abort("Results must be a tibble or a list of tibbles.")
  }

  class(new_results) <- old_class
  return(new_results)
}

#' Make comparisons pairs for multi-group analysis
#'
#' @param levels A character vector of group levels.
#' @param reverse A logical value. Default is FALSE. If TRUE, reference group will be the second level.
#' @return A list of character vectors, each containing two group levels.
#' @noRd
.make_comparisons <- function(levels, reverse = FALSE) {
  n_pairs <- length(levels) * (length(levels) - 1) / 2
  pairs <- vector("list", n_pairs)
  count <- 1
  for (i in 1:(length(levels) - 1)) {
    for (j in (i + 1):length(levels)) {
      if (reverse) {
        pairs[[count]] <- c(levels[j], levels[i])
      } else {
        pairs[[count]] <- c(levels[i], levels[j])
      }
      count <- count + 1
    }
  }
  pairs
}

# Extract subject identifiers for paired analyses
.extract_paired_subjects <- function(sample_info, subject_col, group_col) {
  if (is.null(subject_col)) {
    return(NULL)
  }
  if (identical(subject_col, group_col)) {
    cli::cli_abort("{.arg subject_col} cannot equal {.arg group_col}.")
  }
  if (!subject_col %in% colnames(sample_info)) {
    cli::cli_abort(c(
      "Column {.field {subject_col}} not found in sample information.",
      "i" = "Available columns: {.field {colnames(sample_info)}}."
    ))
  }

  subjects <- sample_info[[subject_col]]
  if (any(is.na(subjects))) {
    cli::cli_abort(
      "Column {.field {subject_col}} cannot contain missing values."
    )
  }
  as.character(subjects)
}

# Match samples into complete subject blocks for paired analyses
.match_paired_samples <- function(expr_mat, groups, subjects, method) {
  if (is.null(subjects)) {
    return(list(
      expr_mat = expr_mat,
      groups = groups,
      subjects = NULL,
      paired = FALSE
    ))
  }
  if (length(subjects) != ncol(expr_mat)) {
    cli::cli_abort(
      "subjects must have {.val {ncol(expr_mat)}} values to match expr_mat columns."
    )
  }

  groups <- droplevels(groups)

  subject_group <- paste(subjects, groups, sep = "\r")
  if (anyDuplicated(subject_group) > 0) {
    cli::cli_abort(
      "Each subject must have at most one sample in each group for paired {method}."
    )
  }

  group_levels <- levels(groups)
  subjects_by_group <- purrr::map(
    group_levels,
    ~ subjects[groups == .x]
  )
  complete_subjects <- Reduce(intersect, subjects_by_group)
  if (length(complete_subjects) < 2) {
    cli::cli_abort(
      "At least 2 subjects must have samples in every group for paired {method}."
    )
  }

  indices <- purrr::map(group_levels, function(group) {
    group_indices <- which(groups == group)
    group_indices[match(complete_subjects, subjects[group_indices])]
  })
  sample_indices <- unlist(indices, use.names = FALSE)
  n_subjects <- length(complete_subjects)

  list(
    expr_mat = expr_mat[, sample_indices, drop = FALSE],
    groups = factor(
      rep(group_levels, each = n_subjects),
      levels = group_levels
    ),
    subjects = factor(
      rep(complete_subjects, times = length(group_levels)),
      levels = complete_subjects
    ),
    paired = TRUE
  )
}

# Return aligned complete values for one variable in a two-group paired design
.paired_variable_values <- function(expr_mat, groups, var_name) {
  group_levels <- levels(groups)
  ref <- expr_mat[var_name, groups == group_levels[1]]
  test <- expr_mat[var_name, groups == group_levels[2]]
  complete <- is.finite(ref) & is.finite(test)
  list(ref = ref[complete], test = test[complete])
}
