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
  groups
}

# Convert groups to factor (for gly_xxx_() functions)
.convert_groups_to_factor <- function(groups) {
  if (is.character(groups)) {
    groups <- factor(groups)
  } else if (!is.factor(groups)) {
    cli::cli_abort("groups must be a factor or character vector")
  }
  groups
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
  if (!.has_variable_column(tbl)) {
    return(tbl)
  }

  var_info <- glyexp::get_var_info(exp)
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
  if (!.has_sample_column(tbl)) {
    return(tbl)
  }

  sample_info <- glyexp::get_sample_info(exp)
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
.process_tibble_add_info <- function(tbl, exp, add_info) {
  if (!add_info) {
    return(tbl)
  }

  # Add variable info if tibble has variable column
  if (.has_variable_column(tbl)) {
    tbl <- .add_variable_info(tbl, exp)
  }

  # Add sample info if tibble has sample column
  if (.has_sample_column(tbl)) {
    tbl <- .add_sample_info(tbl, exp)
  }

  return(tbl)
}

# Process results with add_info logic (handles both single tibbles and lists)
.process_results_add_info <- function(results, exp, add_info) {
  if (!add_info) {
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
