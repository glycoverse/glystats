#' Match variables to a group template
#'
#' Assign each sample the template score for its experimental group and
#' correlate each variable's abundance with those scores. Positive and
#' negative coefficients identify matches in opposite directions.
#'
#' @inheritParams gly_cor
#' @inheritParams gly_scale
#' @param template A named numeric vector with one finite score for each
#'   observed group in `group_col`. Names must match the group labels exactly.
#'   Groups may share a score. At least two distinct scores are required.
#' @param group_col Name of the sample metadata column containing the group
#'   labels. Each sample must have a non-missing, non-empty label.
#' @param method Correlation coefficient, one of `"pearson"` (default),
#'   `"spearman"`, or `"kendall"`.
#' @param add_info Whether to include variable annotations in the tidy result.
#'
#' @details
#' Each finite sample observation contributes once. Groups with more observed
#' samples have more weight. For each variable, the test uses the groups and
#' scores represented among its finite abundances on the analysis scale. If a
#' group is entirely missing, the coefficient describes the remaining groups.
#' A test requires at least three such samples, two distinct template scores,
#' and variation in abundance. Other rows retain their counts and receive
#' `NA` for the coefficient and p-values.
#'
#' Raw abundances use `log2(x + 1e-6)` by default. See [gly_scale] for scale
#' overrides. With relative abundances, a match describes a change in relative
#' abundance. Even a strong positive correlation with an increasing template
#' can include a decrease between adjacent groups.
#'
#' The two-sided tests assume independent samples and a prespecified template.
#' For fixed scores, the Pearson test is the slope test in a linear regression
#' with independent normal errors of constant variance. The rank tests assume
#' exchangeability under independence. Spearman uses average ranks for ties
#' and Kendall uses tau-b. With ties, [stats::cor.test()] uses approximate
#' rank-test p-values, which can be inaccurate with few samples. Without ties,
#' it uses its usual sample-size-dependent calculation.
#'
#' Multiplicity adjustment covers the testable variables for this template.
#' If several templates are screened, account for the full set of tests.
#' A template chosen after inspecting the data needs independent validation.
#'
#' @returns A `glystats_pattern_hunter_res` object inheriting from
#'   `glystats_res`. [get_tidy_result()] returns one row per input variable
#'   in assay order.
#'
#'   * `.variable` gives the SE variable ID (`variable` for legacy inputs).
#'   * `cor` is the signed correlation coefficient.
#'   * `n_samples`, `n_groups`, and `n_scores` count the finite samples,
#'     represented groups, and distinct template scores for the variable.
#'   * `p_val` is the two-sided p-value and `p_adj` is its adjusted value.
#'     Setting `p_adj_method = NULL` omits `p_adj`.
#'
#'   [get_raw_result()] returns a list of [stats::cor.test()] results named by
#'   variable ID, with `NULL` for untestable variables. The result records
#'   `template`, `group_col`, `method`, `p_adj_method`, the input `meta_data`,
#'   and the `preprocessing` used. SE inputs also retain the `identifiers`
#'   mapping described in [gly_identity].
#'
#' @references
#' Pavlidis P, Noble WS (2001). Analysis of strain and regional variation in
#' gene expression in mouse brain. *Genome Biology*, 2, research0042.
#' \doi{10.1186/gb-2001-2-10-research0042}.
#'
#' Pavlidis P (2003). Using ANOVA for gene selection from microarray studies
#' of the nervous system. *Methods*, 31, 282-289.
#' \doi{10.1016/S1046-2023(03)00157-9}.
#'
#' For a related pattern search in metabolomics, see Xia J, Mandal R,
#' Sinelnikov IV, Broadhurst D, Wishart DS (2012). MetaboAnalyst 2.0---a
#' comprehensive server for metabolomic data analysis. *Nucleic Acids
#' Research*, 40, W127-W133. \doi{10.1093/nar/gks374}.
#'
#' @examples
#' exp <- glyexp::real_experiment2
#' groups <- SummarizedExperiment::colData(exp)$group
#' exp <- exp[1:5, groups %in% c("H", "C")]
#' res <- gly_pattern_hunter(exp, template = c(H = 0, C = 1), add_info = FALSE)
#' get_tidy_result(res)
#' @export
gly_pattern_hunter <- function(
  exp,
  template,
  group_col = "group",
  method = "pearson",
  p_adj_method = "BH",
  add_info = TRUE,
  log_transform = NULL
) {
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_choice(method, c("pearson", "spearman", "kendall"))
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_flag(add_info)
  meta_data <- .get_meta_data(exp)
  transformation <- .resolve_transformation(
    exp,
    log_transform,
    meta_data = meta_data
  )

  group_info <- .extract_and_validate_groups(
    .get_sample_info(exp),
    group_col,
    min_count = 2L,
    show_info = FALSE
  )
  groups <- as.character(group_info$groups)
  target <- .pattern_target(template, groups, group_col)
  expr_mat <- .get_expr_mat(exp)
  checkmate::assert_matrix(expr_mat, mode = "numeric", ncols = length(groups))
  expr_mat <- .log_transform_expr_mat(expr_mat, transformation$log_transform)

  tests <- .test_pattern(expr_mat, groups, target, method)
  tidy <- tests$tidy_result
  if (!is.null(p_adj_method)) {
    tidy$p_adj <- stats::p.adjust(tidy$p_val, method = p_adj_method)
  }
  result <- list(
    tidy_result = .process_results_add_info(tidy, exp, add_info),
    raw_result = tests$raw_result,
    template = template,
    group_col = group_col,
    method = method,
    p_adj_method = p_adj_method,
    meta_data = meta_data,
    preprocessing = transformation
  )
  class(result) <- c("glystats_pattern_hunter_res", "glystats_res")
  .record_se_identity(result, exp)
}

.pattern_target <- function(template, groups, group_col) {
  checkmate::assert_vector(template, strict = TRUE)
  checkmate::assert_numeric(
    template,
    finite = TRUE,
    any.missing = FALSE,
    names = "unique"
  )
  if (anyNA(groups) || any(!nzchar(groups))) {
    cli::cli_abort(
      "Column {.field {group_col}} contains missing or empty group labels."
    )
  }
  if (!setequal(names(template), groups)) {
    cli::cli_abort(c(
      "{.arg template} names must match the observed groups in {.field {group_col}}.",
      "i" = "Observed groups are {.val {unique(groups)}}.",
      "i" = "Template names are {.val {names(template)}}."
    ))
  }
  if (length(unique(template)) < 2L) {
    cli::cli_abort("{.arg template} must have at least two distinct scores.")
  }
  unname(template[groups])
}

.test_pattern <- function(expr_mat, groups, target, method) {
  ids <- as.character(rownames(expr_mat))
  fits <- stats::setNames(vector("list", length(ids)), ids)
  cor <- p_val <- rep(NA_real_, length(ids))
  n_samples <- n_groups <- n_scores <- integer(length(ids))

  for (i in seq_along(ids)) {
    keep <- is.finite(expr_mat[i, ])
    values <- expr_mat[i, keep]
    scores <- target[keep]
    n_samples[i] <- length(values)
    n_groups[i] <- length(unique(groups[keep]))
    n_scores[i] <- length(unique(scores))
    n_values <- length(unique(values))
    if (n_samples[i] < 3L || n_scores[i] < 2L || n_values < 2L) {
      next
    }
    ties <- n_scores[i] < n_samples[i] || n_values < n_samples[i]
    fit <- stats::cor.test(
      values,
      scores,
      method = method,
      exact = if (ties) FALSE else NULL
    )
    cor[i] <- fit$estimate
    p_val[i] <- fit$p.value
    fits[[i]] <- fit
  }

  tidy <- tibble::tibble(
    variable = ids,
    cor = cor,
    n_samples = n_samples,
    n_groups = n_groups,
    n_scores = n_scores,
    p_val = p_val
  )
  list(tidy_result = tidy, raw_result = fits)
}
