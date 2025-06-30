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
