# Check if required packages are available
.check_pkg_available <- function(pkg, install_hint = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is.null(install_hint)) {
      install_hint <- stringr::str_glue("install.packages('{pkg}')")
    }
    cli::cli_abort(c(
      "Package {.pkg {pkg}} is required for this analysis.",
      "i" = "Install it with: {.code {install_hint}}"
    ))
  }
}

# Group extraction, validation, and conversion helper functions ---------------

# Check if group column exists in sample information
.check_group_column_exists <- function(sample_info, group_col) {
  if (!group_col %in% colnames(sample_info)) {
    cli::cli_abort("Column {.field {group_col}} not found in sample information")
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

# Helper function to generate validation error message
.generate_validation_error <- function(group_col, method, groups, message_type, count) {
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
.validate_group_count <- function(groups, group_col, min_count = NULL, max_count = NULL, method = NULL) {
  # Skip validation if both limits are NULL
  if (is.null(min_count) && is.null(max_count)) {
    return()
  }
  
  n_groups <- length(levels(groups))
  
  # Handle exact count case (when min_count == max_count)
  if (!is.null(min_count) && !is.null(max_count) && min_count == max_count) {
    if (n_groups != min_count) {
      .generate_validation_error(group_col, method, groups, "exactly", min_count)
    }
    return()  # Early return, no need to check min/max separately
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
  cli::cli_alert_info("Group 1: {.val {levels(groups)[1]}}")
  cli::cli_alert_info("Group 2: {.val {levels(groups)[2]}}")
}

# Display group information for multi-group analysis
.display_multi_group_info <- function(groups) {
  n_groups <- length(levels(groups))
  cli::cli_alert_info("Number of groups: {.val {n_groups}}")
  cli::cli_alert_info("Groups: {.val {levels(groups)}}")
}

# Extract and validate groups (comprehensive function)
.extract_and_validate_groups <- function(sample_info, group_col, min_count = NULL, max_count = NULL,
                                        method = NULL, show_info = TRUE) {
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
  var_info_subset <- var_info[, !colnames(var_info) %in% "variable", drop = FALSE]

  # Only join if there are columns other than variable
  if (ncol(var_info_subset) > 0) {
    # Put var_info columns first, then the original tibble columns
    result <- dplyr::left_join(var_info, tbl, by = "variable")
    return(result)
  }

  return(tbl)
}

# Add sample information to a tibble
.add_sample_info <- function(tbl, exp) {
  if (!.has_sample_column(tbl)) {
    return(tbl)
  }

  sample_info <- glyexp::get_sample_info(exp)
  # Remove the sample column from sample_info to avoid duplication
  sample_info_subset <- sample_info[, !colnames(sample_info) %in% "sample", drop = FALSE]

  # Only join if there are columns other than sample
  if (ncol(sample_info_subset) > 0) {
    # Put sample_info columns first, then the original tibble columns
    result <- dplyr::left_join(sample_info, tbl, by = "sample")
    return(result)
  }

  return(tbl)
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

  # If results is a tibble, process it directly
  if (tibble::is_tibble(results)) {
    return(.process_tibble_add_info(results, exp, add_info))
  }

  # If results is a list, process each tibble in the list
  if (is.list(results)) {
    processed_results <- purrr::map(results, function(item) {
      if (tibble::is_tibble(item)) {
        .process_tibble_add_info(item, exp, add_info)
      } else {
        item
      }
    })
    return(processed_results)
  }

  # If results is neither tibble nor list, return as is
  return(results)
}
