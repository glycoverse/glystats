#' Construct and Test Correlated Variable Sets
#'
#' Construct and jointly test correlated variable sets with Hotelling's
#' \eqn{T^2} statistic.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing expression data and sample
#'   information.
#' @param sets A uniquely named list of character vectors defining custom sets.
#'   When `NULL` (the default), correlated sets are constructed automatically
#'   from `exp` using `threshold`, `correlation`, `clustering`, and `within`.
#' @param group_col A character string naming the two-level grouping column in
#'   sample information. Default is `"group"`.
#' @param subject_col An optional character string naming subject identifiers.
#'   When supplied, one-sample Hotelling tests are applied to paired
#'   test-minus-reference difference vectors. Otherwise independent two-sample
#'   Hotelling tests are used.
#' @param method Statistical method. Currently only `"hotelling"` is
#'   supported.
#' @param p_adj_method A character string specifying the method used to adjust
#'   set-level p-values. See [stats::p.adjust()]. Default is `"BH"`. If `NULL`,
#'   no adjustment is performed.
#' @param threshold A numeric correlation threshold between 0 and 1 used when
#'   `sets` is `NULL`. Variables are linked only when their correlation is
#'   strictly greater than this value. Default is `0.9`.
#' @param correlation Correlation coefficient used when `sets` is `NULL`.
#'   Either `"pearson"` (default) or `"spearman"`.
#' @param clustering How correlated variables are grouped when `sets` is
#'   `NULL`. `"connected"` (default) allows transitive chaining; `"complete"`
#'   requires every pair in a set to exceed `threshold`.
#' @param within Optional variable-information columns defining strata within
#'   which automatic sets are constructed, such as `c("protein",
#'   "protein_site")`.
#'
#' @details
#' Expression values are transformed using `log2(x + 1e-6)`. The first factor
#' level of `group_col` is the reference group and the second is the test group.
#' The estimate is therefore `test - reference`.
#'
#' Every usable variable is assigned to a set, including a one-variable set when
#' it is not correlated above `threshold` with another variable. Structurally
#' rank-deficient covariance matrices are tested in their nonredundant subspace
#' using a Moore-Penrose inverse when the original set dimension satisfies the
#' classical sample-size requirement and the observed mean contrast lies in the
#' estimable subspace. Other ineligible sets are retained with
#' `status = "failed"`, `NA` statistics, and an explicit `failure_reason`.
#'
#' @returns A list with classes `glystats_set_test_res` and `glystats_res`
#' containing:
#' - `tidy_result$sets`: One row per set with estimates for every member, the
#'   effective covariance rank in `test_dimension`, Hotelling statistic,
#'   degrees of freedom, Mahalanobis effect size, p-value, adjusted p-value, and
#'   fit status.
#' - `tidy_result$members`: One row per set member with its marginal estimate
#'   and mean within-set correlation.
#' - `raw_result`: Set definitions, their correlation matrix, individual
#'   Hotelling test objects, and automatic set-construction details (when
#'   `sets` is `NULL`).
#' - `meta_data`: Metadata copied from `exp`.
#'
#' @examples
#' set.seed(1)
#' expression <- matrix(
#'   stats::rexp(36),
#'   nrow = 3,
#'   dimnames = list(c("A", "B", "C"), paste0("S", 1:12))
#' )
#' expression["B", ] <- expression["A", ] * exp(stats::rnorm(12, sd = 0.01))
#' exp <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(expression = expression),
#'   colData = S4Vectors::DataFrame(
#'     group = factor(
#'       rep(c("control", "case"), each = 6),
#'       levels = c("control", "case")
#'     )
#'   )
#' )
#' result <- gly_set_test(exp, threshold = 0.8)
#' result$tidy_result$sets
#'
#' @seealso [stats::cor()], [stats::hclust()], [stats::p.adjust()]
#' @export
gly_set_test <- function(
  exp,
  sets = NULL,
  group_col = "group",
  subject_col = NULL,
  method = "hotelling",
  p_adj_method = "BH",
  threshold = 0.9,
  correlation = "pearson",
  clustering = "connected",
  within = NULL
) {
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_string(subject_col, null.ok = TRUE)
  checkmate::assert_choice(method, "hotelling")
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_number(threshold, lower = 0, upper = 1)
  checkmate::assert_choice(correlation, c("pearson", "spearman"))
  checkmate::assert_choice(clustering, c("complete", "connected"))
  checkmate::assert_character(within, unique = TRUE, null.ok = TRUE)

  expr_mat <- .get_expr_mat(exp)
  sample_info <- .get_sample_info(exp)
  var_info <- .get_var_info(exp)
  .validate_set_variable_names(expr_mat)
  if (is.null(sets)) {
    set_construction <- .construct_correlated_sets(
      expr_mat = expr_mat,
      strata = .set_strata(var_info, within),
      threshold = threshold,
      correlation = correlation,
      clustering = clustering
    )
    set_construction["within"] <- list(within)
    set_info <- list(
      sets = set_construction$sets,
      correlation_matrix = set_construction$correlation_matrix
    )
  } else {
    set_construction <- NULL
    set_info <- .normalize_set_definitions(sets, rownames(expr_mat))
  }

  group_info <- .extract_and_validate_groups(
    sample_info = sample_info,
    group_col = group_col,
    min_count = 2,
    max_count = 2,
    method = "Hotelling set test"
  )
  groups <- group_info$groups
  if (any(is.na(groups))) {
    cli::cli_abort("Column {.field {group_col}} cannot contain missing values.")
  }

  subjects <- .set_test_subjects(sample_info, subject_col, group_col)
  result <- .analyze_set_test(
    expr_mat = expr_mat,
    sets = set_info$sets,
    groups = groups,
    subjects = subjects,
    p_adj_method = p_adj_method,
    correlation_matrix = set_info$correlation_matrix
  )
  result$raw_result$definitions <- set_info$sets
  result$raw_result$set_construction <- set_construction
  result$meta_data <- .get_meta_data(exp)
  result
}

.validate_set_variable_names <- function(expr_mat) {
  variable_names <- rownames(expr_mat)
  if (
    is.null(variable_names) ||
      any(is.na(variable_names)) ||
      any(variable_names == "") ||
      anyDuplicated(variable_names) > 0
  ) {
    cli::cli_abort(
      "Expression variables must have unique, non-missing row names."
    )
  }
  invisible(NULL)
}

.set_strata <- function(var_info, within) {
  if (is.null(within)) {
    return(rep("all", nrow(var_info)))
  }

  missing_columns <- setdiff(within, colnames(var_info))
  if (length(missing_columns) > 0) {
    cli::cli_abort(c(
      "Columns in {.arg within} were not found in variable information.",
      "x" = "Missing columns: {.field {missing_columns}}.",
      "i" = "Available columns: {.field {colnames(var_info)}}."
    ))
  }

  values <- var_info[, within, drop = FALSE]
  values[] <- lapply(values, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  })
  do.call(
    interaction,
    c(as.list(values), list(drop = TRUE, lex.order = TRUE, sep = "\r"))
  ) |>
    as.character()
}

.construct_correlated_sets <- function(
  expr_mat,
  strata,
  threshold,
  correlation,
  clustering
) {
  log_expr_mat <- suppressWarnings(.log_transform_expr_mat(expr_mat))
  usable <- .classify_set_profiles(log_expr_mat)
  usable_names <- usable$variable[is.na(usable$reason)]
  exclusions <- dplyr::filter(usable, !is.na(.data$reason))

  correlation_matrix <- matrix(
    NA_real_,
    nrow = length(usable_names),
    ncol = length(usable_names),
    dimnames = list(usable_names, usable_names)
  )
  set_groups <- list()

  for (stratum in unique(strata)) {
    stratum_variables <- intersect(
      rownames(expr_mat)[strata == stratum],
      usable_names
    )
    if (length(stratum_variables) == 0) {
      next
    }

    stratum_matrix <- log_expr_mat[stratum_variables, , drop = FALSE]
    stratum_correlation <- .profile_correlation(stratum_matrix, correlation)
    correlation_matrix[
      stratum_variables,
      stratum_variables
    ] <- stratum_correlation

    clusters <- .cluster_correlated_profiles(
      stratum_correlation,
      threshold,
      clustering
    )
    set_groups <- c(set_groups, clusters)
  }

  set_groups <- set_groups[
    order(vapply(
      set_groups,
      function(x) min(match(x, rownames(expr_mat))),
      numeric(1)
    ))
  ]
  if (length(set_groups) > 0) {
    names(set_groups) <- paste0(
      "set_",
      seq_along(set_groups)
    )
  }

  sets <- lapply(set_groups, function(members) {
    rownames(expr_mat)[rownames(expr_mat) %in% members]
  })

  if (length(sets) == 0) {
    membership <- tibble::tibble(
      set_id = character(),
      variable = character()
    )
  } else {
    membership <- purrr::imap_dfr(sets, function(members, set_id) {
      tibble::tibble(
        set_id = set_id,
        variable = members
      )
    })
  }

  list(
    sets = sets,
    membership = membership,
    correlation_matrix = correlation_matrix,
    excluded_variables = exclusions,
    threshold = threshold,
    correlation = correlation,
    clustering = clustering
  )
}

.classify_set_profiles <- function(log_expr_mat) {
  tibble::tibble(
    variable = rownames(log_expr_mat),
    reason = vapply(
      seq_len(nrow(log_expr_mat)),
      function(i) {
        values <- log_expr_mat[i, ]
        values <- values[is.finite(values)]
        if (length(values) == 0) {
          "all_missing"
        } else if (length(values) < 2) {
          "insufficient_observations"
        } else if (stats::var(values) <= sqrt(.Machine$double.eps)) {
          "zero_variance"
        } else {
          NA_character_
        }
      },
      character(1)
    )
  )
}

.profile_correlation <- function(profile_matrix, method) {
  if (nrow(profile_matrix) == 1) {
    return(matrix(
      1,
      nrow = 1,
      ncol = 1,
      dimnames = list(rownames(profile_matrix), rownames(profile_matrix))
    ))
  }
  stats::cor(
    t(profile_matrix),
    use = "pairwise.complete.obs",
    method = method
  )
}

.cluster_correlated_profiles <- function(
  correlation_matrix,
  threshold,
  method
) {
  variables <- rownames(correlation_matrix)
  if (length(variables) == 1) {
    return(stats::setNames(list(variables), variables))
  }

  linked <- !is.na(correlation_matrix) & correlation_matrix > threshold
  diag(linked) <- TRUE
  if (method == "connected") {
    return(.connected_components(linked))
  }

  distance <- 1 - correlation_matrix
  distance[!linked] <- 2
  diag(distance) <- 0
  tree <- stats::hclust(stats::as.dist(distance), method = "complete")
  groups <- stats::cutree(tree, h = 1 - threshold)
  split(variables, groups)
}

.connected_components <- function(adjacency) {
  variables <- rownames(adjacency)
  visited <- stats::setNames(rep(FALSE, length(variables)), variables)
  components <- list()

  for (variable in variables) {
    if (visited[[variable]]) {
      next
    }
    queue <- variable
    visited[[variable]] <- TRUE
    component <- character()
    while (length(queue) > 0) {
      current <- queue[[1]]
      queue <- queue[-1]
      component <- c(component, current)
      neighbors <- variables[adjacency[current, ] & !visited]
      if (length(neighbors) > 0) {
        visited[neighbors] <- TRUE
        queue <- c(queue, neighbors)
      }
    }
    components[[length(components) + 1L]] <- component
  }
  components
}

.normalize_set_definitions <- function(sets, variables) {
  if (inherits(sets, "glystats_correlated_sets")) {
    definitions <- sets$sets
    correlation_matrix <- sets$correlation_matrix
  } else if (is.list(sets)) {
    definitions <- sets
    correlation_matrix <- NULL
  } else {
    cli::cli_abort(
      "{.arg sets} must be a named list of character vectors."
    )
  }

  if (length(definitions) == 0 && inherits(sets, "glystats_correlated_sets")) {
    return(list(sets = definitions, correlation_matrix = correlation_matrix))
  }

  if (
    is.null(names(definitions)) ||
      any(is.na(names(definitions))) ||
      any(names(definitions) == "") ||
      anyDuplicated(names(definitions)) > 0
  ) {
    cli::cli_abort("{.arg sets} must have unique, non-missing names.")
  }
  if (!all(vapply(definitions, is.character, logical(1)))) {
    cli::cli_abort("Every element of {.arg sets} must be a character vector.")
  }
  if (any(lengths(definitions) < 1)) {
    cli::cli_abort(
      "Every element of {.arg sets} must contain at least 1 variable."
    )
  }
  duplicated_members <- vapply(
    definitions,
    function(x) anyDuplicated(x) > 0,
    logical(1)
  )
  if (any(duplicated_members)) {
    cli::cli_abort(
      "Variables within each element of {.arg sets} must be unique."
    )
  }

  missing_variables <- setdiff(unique(unlist(definitions)), variables)
  if (length(missing_variables) > 0) {
    cli::cli_abort(c(
      "Some variables in {.arg sets} were not found in {.arg exp}.",
      "x" = "Missing variables: {.val {missing_variables}}."
    ))
  }

  list(sets = definitions, correlation_matrix = correlation_matrix)
}

.set_test_subjects <- function(sample_info, subject_col, group_col) {
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

.analyze_set_test <- function(
  expr_mat,
  sets,
  groups,
  subjects = NULL,
  p_adj_method = "BH",
  correlation_matrix = NULL
) {
  log_expr_mat <- suppressWarnings(.log_transform_expr_mat(expr_mat))
  design <- .set_test_design(groups, subjects)

  if (length(sets) == 0) {
    set_results <- .empty_set_results(p_adj_method)
    member_results <- tibble::tibble(
      set_id = character(),
      variable = character(),
      marginal_estimate = numeric(),
      correlation_summary = numeric()
    )
    return(structure(
      list(
        tidy_result = list(sets = set_results, members = member_results),
        raw_result = list(
          tests = list(),
          correlation_matrix = correlation_matrix
        )
      ),
      class = c("glystats_set_test_res", "glystats_res")
    ))
  }

  tests <- purrr::imap(sets, function(members, set_id) {
    .test_one_set(
      log_expr_mat = log_expr_mat,
      members = members,
      design = design,
      set_id = set_id
    )
  })

  set_results <- purrr::map_dfr(tests, "tidy")
  if (!is.null(p_adj_method)) {
    set_results$p_adj <- stats::p.adjust(
      set_results$p_val,
      method = p_adj_method
    )
  }

  member_results <- purrr::imap_dfr(sets, function(members, set_id) {
    test <- tests[[set_id]]
    correlation <- .set_member_correlation(
      log_expr_mat,
      members,
      correlation_matrix
    )
    tibble::tibble(
      set_id = set_id,
      variable = members,
      marginal_estimate = unname(test$member_estimates[members]),
      correlation_summary = unname(correlation[members])
    )
  })

  structure(
    list(
      tidy_result = list(
        sets = set_results,
        members = member_results
      ),
      raw_result = list(
        tests = purrr::map(tests, "raw"),
        correlation_matrix = correlation_matrix
      )
    ),
    class = c("glystats_set_test_res", "glystats_res")
  )
}

.empty_set_results <- function(p_adj_method) {
  result <- tibble::tibble(
    set_id = character(),
    n_variables = integer(),
    test_dimension = integer(),
    variables = list(),
    estimate = list(),
    statistic = numeric(),
    df1 = numeric(),
    df2 = numeric(),
    effect_size = numeric(),
    p_val = numeric(),
    n_ref = integer(),
    n_test = integer(),
    ref_group = character(),
    test_group = character(),
    paired = logical(),
    status = character(),
    failure_reason = character()
  )
  if (!is.null(p_adj_method)) {
    result$p_adj <- numeric()
  }
  result
}

.set_test_design <- function(groups, subjects) {
  group_levels <- levels(groups)
  if (is.null(subjects)) {
    return(list(
      paired = FALSE,
      ref_group = group_levels[[1]],
      test_group = group_levels[[2]],
      ref_indices = which(groups == group_levels[[1]]),
      test_indices = which(groups == group_levels[[2]])
    ))
  }

  subject_group <- paste(subjects, groups, sep = "\r")
  if (anyDuplicated(subject_group) > 0) {
    cli::cli_abort(
      "Each subject must have at most one sample in each group for paired set testing."
    )
  }
  ref_indices <- which(groups == group_levels[[1]])
  test_indices <- which(groups == group_levels[[2]])
  complete_subjects <- intersect(subjects[ref_indices], subjects[test_indices])
  if (length(complete_subjects) == 0) {
    cli::cli_abort("No subjects have samples in both groups.")
  }

  list(
    paired = TRUE,
    ref_group = group_levels[[1]],
    test_group = group_levels[[2]],
    ref_indices = ref_indices[match(complete_subjects, subjects[ref_indices])],
    test_indices = test_indices[match(
      complete_subjects,
      subjects[test_indices]
    )],
    subjects = complete_subjects
  )
}

.test_one_set <- function(log_expr_mat, members, design, set_id) {
  member_matrix <- log_expr_mat[members, , drop = FALSE]
  test_matrix <- t(member_matrix)

  if (design$paired) {
    ref <- test_matrix[design$ref_indices, , drop = FALSE]
    test <- test_matrix[design$test_indices, , drop = FALSE]
    complete <- stats::complete.cases(ref) & stats::complete.cases(test)
    ref <- ref[complete, , drop = FALSE]
    test <- test[complete, , drop = FALSE]
    differences <- test - ref
    hotelling <- .hotelling_one_sample(differences)
    included_subjects <- design$subjects[complete]
  } else {
    ref <- test_matrix[design$ref_indices, , drop = FALSE]
    test <- test_matrix[design$test_indices, , drop = FALSE]
    ref <- ref[stats::complete.cases(ref), , drop = FALSE]
    test <- test[stats::complete.cases(test), , drop = FALSE]
    hotelling <- .hotelling_two_sample(ref, test)
    included_subjects <- NULL
  }

  member_estimates <- .set_member_estimates(member_matrix, design)
  if (is.null(hotelling$estimate)) {
    set_estimates <- member_estimates[members]
  } else {
    set_estimates <- stats::setNames(
      as.numeric(hotelling$estimate),
      members
    )
  }
  tidy <- tibble::tibble(
    set_id = set_id,
    n_variables = length(members),
    test_dimension = as.integer(hotelling$df1),
    variables = list(members),
    estimate = list(set_estimates),
    statistic = hotelling$statistic,
    df1 = hotelling$df1,
    df2 = hotelling$df2,
    effect_size = hotelling$effect_size,
    p_val = hotelling$p_val,
    n_ref = nrow(ref),
    n_test = nrow(test),
    ref_group = design$ref_group,
    test_group = design$test_group,
    paired = design$paired,
    status = hotelling$status,
    failure_reason = hotelling$failure_reason
  )

  list(
    tidy = tidy,
    member_estimates = member_estimates,
    raw = c(
      hotelling,
      list(
        included_subjects = included_subjects
      )
    )
  )
}

.set_member_estimates <- function(member_matrix, design) {
  if (design$paired) {
    differences <- member_matrix[, design$test_indices, drop = FALSE] -
      member_matrix[, design$ref_indices, drop = FALSE]
    result <- rowMeans(differences, na.rm = TRUE)
    result[is.nan(result)] <- NA_real_
    return(result)
  }

  result <- rowMeans(
    member_matrix[, design$test_indices, drop = FALSE],
    na.rm = TRUE
  ) -
    rowMeans(member_matrix[, design$ref_indices, drop = FALSE], na.rm = TRUE)
  result[is.nan(result)] <- NA_real_
  result
}

.hotelling_one_sample <- function(differences) {
  n <- nrow(differences)
  p <- ncol(differences)
  if (n < 2) {
    return(.failed_hotelling(
      NA_real_,
      NA_real_,
      "At least 2 complete pairs are required."
    ))
  }

  covariance <- stats::cov(differences)
  subspace <- .set_covariance_subspace(covariance, p)
  if (!is.null(subspace$failure_reason)) {
    return(.failed_hotelling(
      subspace$rank,
      if (is.na(subspace$rank)) NA_real_ else n - subspace$rank,
      paste("The paired covariance matrix", subspace$failure_reason)
    ))
  }
  rank <- subspace$rank
  if (n <= p) {
    return(.failed_hotelling(
      rank,
      n - p,
      "The number of variables must be smaller than the number of complete pairs."
    ))
  }

  estimate <- colMeans(differences)
  if (!.estimate_in_covariance_subspace(estimate, subspace$basis)) {
    return(.failed_hotelling(
      rank,
      n - rank,
      "The mean difference lies outside the estimable covariance subspace."
    ))
  }
  distance_squared <- as.numeric(
    t(estimate) %*% subspace$inverse %*% estimate
  )
  statistic <- n * distance_squared
  f_statistic <- (n - rank) * statistic / (rank * (n - 1))

  list(
    statistic = statistic,
    f_statistic = f_statistic,
    df1 = rank,
    df2 = n - rank,
    effect_size = sqrt(max(distance_squared, 0)),
    p_val = stats::pf(f_statistic, rank, n - rank, lower.tail = FALSE),
    status = "ok",
    failure_reason = NA_character_,
    covariance = subspace$covariance,
    estimate = estimate
  )
}

.hotelling_two_sample <- function(ref, test) {
  n_ref <- nrow(ref)
  n_test <- nrow(test)
  p <- ncol(ref)
  if (n_ref < 2 || n_test < 2) {
    return(.failed_hotelling(
      NA_real_,
      NA_real_,
      "At least 2 complete samples are required in each group."
    ))
  }

  covariance <- ((n_ref - 1) *
    stats::cov(ref) +
    (n_test - 1) * stats::cov(test)) /
    (n_ref + n_test - 2)
  subspace <- .set_covariance_subspace(covariance, p)
  if (!is.null(subspace$failure_reason)) {
    return(.failed_hotelling(
      subspace$rank,
      if (is.na(subspace$rank)) {
        NA_real_
      } else {
        n_ref + n_test - subspace$rank - 1
      },
      paste("The pooled covariance matrix", subspace$failure_reason)
    ))
  }
  rank <- subspace$rank
  sample_size_df2 <- n_ref + n_test - p - 1
  if (sample_size_df2 <= 0) {
    return(.failed_hotelling(
      rank,
      sample_size_df2,
      "The number of variables is too large for the complete group sample sizes."
    ))
  }

  df2 <- n_ref + n_test - rank - 1
  estimate <- colMeans(test) - colMeans(ref)
  if (!.estimate_in_covariance_subspace(estimate, subspace$basis)) {
    return(.failed_hotelling(
      rank,
      df2,
      "The mean difference lies outside the estimable covariance subspace."
    ))
  }

  distance_squared <- as.numeric(
    t(estimate) %*% subspace$inverse %*% estimate
  )
  statistic <- n_ref * n_test / (n_ref + n_test) * distance_squared
  f_statistic <- df2 * statistic / ((n_ref + n_test - 2) * rank)

  list(
    statistic = statistic,
    f_statistic = f_statistic,
    df1 = rank,
    df2 = df2,
    effect_size = sqrt(max(distance_squared, 0)),
    p_val = stats::pf(f_statistic, rank, df2, lower.tail = FALSE),
    status = "ok",
    failure_reason = NA_character_,
    covariance = subspace$covariance,
    estimate = estimate
  )
}

.failed_hotelling <- function(df1, df2, reason) {
  list(
    statistic = NA_real_,
    f_statistic = NA_real_,
    df1 = as.numeric(df1),
    df2 = as.numeric(df2),
    effect_size = NA_real_,
    p_val = NA_real_,
    status = "failed",
    failure_reason = reason,
    covariance = NULL,
    estimate = NULL
  )
}

.set_covariance_subspace <- function(covariance, dimension) {
  covariance <- matrix(
    covariance,
    nrow = dimension,
    ncol = dimension,
    dimnames = dimnames(covariance)
  )
  if (any(!is.finite(covariance))) {
    return(list(
      covariance = covariance,
      inverse = NULL,
      basis = NULL,
      rank = NA_integer_,
      failure_reason = "contains non-finite values."
    ))
  }

  covariance <- (covariance + t(covariance)) / 2
  decomposition <- eigen(covariance, symmetric = TRUE)
  eigenvalues <- decomposition$values
  max_eigenvalue <- max(abs(eigenvalues))
  tolerance <- sqrt(.Machine$double.eps) * max_eigenvalue
  if (any(eigenvalues < -tolerance)) {
    return(list(
      covariance = covariance,
      inverse = NULL,
      basis = NULL,
      rank = NA_integer_,
      failure_reason = "is not positive semidefinite."
    ))
  }

  retained <- eigenvalues > tolerance
  rank <- sum(retained)
  if (rank == 0) {
    return(list(
      covariance = covariance,
      inverse = NULL,
      basis = NULL,
      rank = 0L,
      failure_reason = "has no estimable dimensions."
    ))
  }

  vectors <- decomposition$vectors[, retained, drop = FALSE]
  inverse <- sweep(
    vectors,
    MARGIN = 2,
    STATS = eigenvalues[retained],
    FUN = "/"
  ) %*%
    t(vectors)
  dimnames(inverse) <- dimnames(covariance)
  list(
    covariance = covariance,
    inverse = inverse,
    basis = vectors,
    rank = as.integer(rank),
    failure_reason = NULL
  )
}

.estimate_in_covariance_subspace <- function(estimate, basis) {
  projected <- as.vector(basis %*% crossprod(basis, estimate))
  residual_norm <- sqrt(sum((estimate - projected)^2))
  estimate_norm <- sqrt(sum(estimate^2))
  if (estimate_norm == 0) {
    return(TRUE)
  }
  residual_norm <= sqrt(.Machine$double.eps) * estimate_norm
}

.set_member_correlation <- function(
  log_expr_mat,
  members,
  correlation_matrix
) {
  if (
    is.null(correlation_matrix) ||
      !all(members %in% rownames(correlation_matrix))
  ) {
    correlation_matrix <- .profile_correlation(
      log_expr_mat[members, , drop = FALSE],
      "pearson"
    )
  } else {
    correlation_matrix <- correlation_matrix[
      members,
      members,
      drop = FALSE
    ]
  }

  vapply(
    seq_along(members),
    function(i) {
      values <- correlation_matrix[i, -i]
      if (length(values) == 0 || all(is.na(values))) {
        return(NA_real_)
      }
      mean(values, na.rm = TRUE)
    },
    numeric(1)
  ) |>
    stats::setNames(members)
}
