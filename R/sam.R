#' Test for differences between two independent groups using SAM
#'
#' Apply Significance Analysis of Microarrays (SAM) to glycan or glycopeptide
#' abundances. Positive scores indicate higher mean abundance in the test
#' group on the analysis scale.
#'
#' @inheritParams gly_limma
#' @inheritParams gly_scale
#' @param delta A positive threshold or a vector of distinct positive
#'   thresholds for the departure of ordered scores from their permutation
#'   expectations. The default, `NULL`, returns the ten grid entries chosen by
#'   [siggenes::samControl()]. Supply one value to call individual features.
#'   With a grid, `called` is `NA` for every feature.
#' @param n_perm Requested number of label permutations, at least two. Uses
#'   all allocations if fewer exist. For more than 30,000 possible allocations,
#'   siggenes can draw the same allocation more than once.
#' @param s0 Nonnegative regularization constant in analysis units. By
#'   default siggenes estimates it across features. This requires at least
#'   25 distinct testable standard errors and variation in the scores used
#'   to estimate it. A prespecified value also permits analysis of smaller
#'   panels with at least two testable features.
#' @param seed Nonnegative integer seed for the permutations. The caller's
#'   random-number state is restored after fitting.
#'
#' @details
#' Supply two groups with at least two independent samples each. The null
#' distribution assumes sample labels are exchangeable and preserves group
#' sizes. Paired, blocked, covariate-adjusted and multigroup designs are
#' unsupported. For blocking or covariate adjustment, see [gly_limma()].
#'
#' Install the optional Bioconductor backend with
#' `BiocManager::install("siggenes")`. Its dependency `multtest` supplies
#' the observed score components used to check which features are testable.
#'
#' [siggenes::sam()] uses the test-minus-reference mean difference divided
#' by the unpooled standard error plus `s0` (its two-class Welch statistic).
#' Raw assays use `log2(x + 1e-6)` by default. For relative abundances, the
#' sign describes a change relative to the normalization reference. With
#' `logx`, it describes a difference in the supplied coordinates.
#'
#' Only features with finite values in every sample and estimable
#' within-group variation enter the common permutation distribution. The
#' backend treats numerically zero standard errors as untestable, including
#' features constant within both groups. It also reserves the finite value
#' `-93074815` as a missing-data code on the analysis scale. A feature
#' containing that value is excluded before permutation. Filtering is reported, and each
#' excluded feature retains its ID and reason. The p-values, q-values and
#' threshold FDR estimates apply to the retained family. No imputation or
#' fold-change filter is applied.
#'
#' @returns A `glystats_sam_res` object inheriting from `glystats_res`.
#'   [get_tidy_result()] returns one row per input feature in assay order.
#'
#'   * `.variable` is the SE feature ID (`variable` for legacy experiments).
#'     `ref_group` and `test_group` identify the comparison.
#'   * `d` is the signed SAM score. `p_val` is the backend's two-sided
#'     p-value from the pooled permuted scores. `q_value` is its feature
#'     q-value computed from those p-values and the estimated null proportion
#'     (`q.version = 1`).
#'   * `called` uses the inclusive lower and upper score cutoffs for a
#'     single supplied delta. `n_finite` counts finite observations on the
#'     analysis scale. Excluded rows have missing statistics and calls, and
#'     an `exclusion` of `"non_finite"`, `"backend_missing_code"`, or
#'     `"unestimable_variance"`.
#'
#'   `thresholds` has one row per delta, with the lower and upper score
#'   cutoffs, number called (`n_called`), mean number of false calls under
#'   permutation (`expected_false`), estimated null proportion (`p0`), and
#'   the backend's permutation FDR (`fdr`, a fraction). siggenes reports an
#'   FDR of zero when a threshold calls zero features.
#'
#'   `diagnostics` records the input and tested feature counts, fitted
#'   `s0`, and the numbers of actual and distinct permutations. `settings`
#'   records the requested comparison and tuning parameters. [get_raw_result()]
#'   returns the fitted `SAM` object. Its feature labels are the original
#'   one-based row positions so even long SE identifiers map unambiguously.
#'   The result retains the input `meta_data`, the applied `preprocessing`
#'   and the SE `identifiers` mapping described in [gly_identity].
#'
#' @references
#' Tusher VG, Tibshirani R, Chu G (2001). Significance analysis of
#' microarrays applied to the ionizing radiation response. *Proceedings of
#' the National Academy of Sciences*, 98, 5116-5121.
#' \doi{10.1073/pnas.091062498}.
#'
#' [siggenes on Bioconductor](https://bioconductor.org/packages/siggenes).
#'
#' @examplesIf requireNamespace("siggenes", quietly = TRUE)
#' exp <- glyexp::real_experiment2
#' group <- SummarizedExperiment::colData(exp)$group
#' two_groups <- exp[, group %in% c("H", "C")]
#' res <- gly_sam(two_groups, ref_group = "H", delta = 1, seed = 42)
#' res$thresholds
#' head(subset(get_tidy_result(res), called, select = c(.variable, d, q_value)))
#' @export
gly_sam <- function(
  exp,
  group_col = "group",
  ref_group = NULL,
  delta = NULL,
  n_perm = 100L,
  s0 = NULL,
  seed = 1L,
  add_info = TRUE,
  log_transform = NULL
) {
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_flag(add_info)
  checkmate::assert_int(n_perm, lower = 2L)
  checkmate::assert_int(seed, lower = 0L)
  checkmate::assert_number(s0, lower = 0, finite = TRUE, null.ok = TRUE)
  if (!is.null(delta)) {
    checkmate::assert_vector(delta, strict = TRUE)
    checkmate::assert_numeric(
      delta,
      min.len = 1L,
      any.missing = FALSE,
      finite = TRUE,
      unique = TRUE
    )
    if (any(delta <= 0)) {
      cli::cli_abort("{.arg delta} must contain positive thresholds.")
    }
  }
  rlang::check_installed("siggenes", reason = "for two-group SAM")

  meta_data <- .get_meta_data(exp)
  preprocessing <- .resolve_transformation(
    exp,
    log_transform,
    meta_data = meta_data
  )
  mat <- .get_expr_mat(exp)
  groups <- .sam_groups(exp, mat, group_col, ref_group)
  cl <- as.integer(groups) - 1L
  mat <- .log_transform_expr_mat(mat, preprocessing$log_transform)
  features <- .sam_testable_features(mat, cl, s0)
  rows <- features$rows

  draws <- min(n_perm, round(choose(length(cl), sum(cl == 0L))))
  # d.stat() drops a dimension when accumulating false calls for a singleton
  # tail across permutation batches. One batch keeps the tail count intact.
  batch_size <- if (features$singleton_tail) draws else 10L
  withr::local_seed(seed)
  fit <- tryCatch(
    siggenes::sam(
      mat[rows, , drop = FALSE],
      cl,
      gene.names = as.character(rows),
      control = siggenes::samControl(delta = delta),
      B = draws,
      B.more = 0,
      n.subset = batch_size,
      s0 = s0 %||% NA_real_,
      R.fold = 0,
      na.replace = FALSE
    ),
    error = function(err) {
      cli::cli_abort(
        c(
          "siggenes could not fit SAM to the {length(rows)} eligible features.",
          if (is.null(s0)) {
            c(
              i = "If automatic regularization is undefined for this panel, supply a prespecified {.arg s0}."
            )
          }
        ),
        parent = err
      )
    }
  )
  distinct_permutations <- nrow(unique(fit@mat.samp))
  if (distinct_permutations < 2L) {
    cli::cli_abort(
      "SAM drew fewer than two distinct label permutations. Increase {.arg n_perm}."
    )
  }

  grid <- fit@mat.fdr
  thresholds <- tibble::tibble(
    delta = grid[, "Delta"],
    cutoff_lower = grid[, "cutlow"],
    cutoff_upper = grid[, "cutup"],
    n_called = as.integer(grid[, "Called"]),
    expected_false = grid[, "False"],
    p0 = grid[, "p0"],
    fdr = grid[, "FDR"]
  )
  # siggenes truncates feature names at 50 characters. All fit slots are in
  # retained row order, including p- and q-values.
  index <- match(seq_len(nrow(mat)), rows)
  d <- unname(fit@d[index])
  called <- rep(NA, nrow(mat))
  if (length(delta) == 1L) {
    called <- d <= thresholds$cutoff_lower | d >= thresholds$cutoff_upper
  }
  tidy <- tibble::tibble(
    variable = rownames(mat),
    ref_group = levels(groups)[1L],
    test_group = levels(groups)[2L],
    d = d,
    p_val = unname(fit@p.value[index]),
    q_value = unname(fit@q.value[index]),
    called = called,
    n_finite = features$n_finite,
    exclusion = features$exclusion
  )
  result <- list(
    tidy_result = .process_results_add_info(tidy, exp, add_info),
    raw_result = fit,
    thresholds = thresholds,
    diagnostics = list(
      n_input = nrow(mat),
      n_tested = length(rows),
      s0 = unname(fit@s0),
      n_permutations = nrow(fit@mat.samp),
      n_distinct_permutations = distinct_permutations
    ),
    settings = list(
      group_col = group_col,
      ref_group = levels(groups)[1L],
      test_group = levels(groups)[2L],
      delta = delta,
      n_perm = n_perm,
      s0 = s0,
      seed = seed
    ),
    meta_data = meta_data,
    preprocessing = preprocessing
  )
  class(result) <- c("glystats_sam_res", "glystats_res")
  .record_se_identity(result, exp)
}

.sam_groups <- function(exp, mat, group_col, ref_group) {
  info <- .get_sample_info(exp)
  checkmate::assert_matrix(mat, mode = "numeric", ncols = nrow(info))
  if (inherits(exp, "glyexp_experiment")) {
    ids <- as.character(info[["sample"]])
    if (
      !checkmate::test_character(
        ids,
        len = ncol(mat),
        unique = TRUE,
        min.chars = 1L,
        any.missing = FALSE
      ) ||
        !setequal(ids, colnames(mat))
    ) {
      cli::cli_abort(
        "Sample metadata IDs must match the assay columns uniquely."
      )
    }
    info <- info[match(colnames(mat), ids), , drop = FALSE]
  }
  .check_group_column_exists(info, group_col)
  groups <- .extract_groups(info, group_col)
  checkmate::assert_factor(groups, len = ncol(mat))
  if (
    anyNA(groups) ||
      anyNA(levels(groups)) ||
      any(!nzchar(trimws(levels(groups))))
  ) {
    cli::cli_abort("Every sample needs a non-missing, non-empty group label.")
  }
  .validate_group_count(
    groups,
    group_col,
    min_count = 2L,
    max_count = 2L,
    method = "unpaired SAM"
  )
  if (any(table(groups) < 2L)) {
    cli::cli_abort("Unpaired SAM requires at least two samples in each group.")
  }
  checkmate::assert_choice(ref_group, levels(groups), null.ok = TRUE)
  .reorder_groups_for_ref(groups, ref_group %||% levels(groups)[1L])
}

.sam_testable_features <- function(mat, cl, s0) {
  n_finite <- as.integer(rowSums(is.finite(mat)))
  # multtest and siggenes interpret this literal as missing even in a
  # complete numeric matrix, which would otherwise drop observations silently.
  reserved <- rowSums(mat == -93074815, na.rm = TRUE) > 0
  complete <- which(n_finite == ncol(mat) & !reserved)
  if (length(complete) < 2L) {
    cli::cli_abort(
      "SAM requires at least two complete features. Found {length(complete)} of {nrow(mat)}."
    )
  }
  # Use the same observed score components and zero check as siggenes::d.stat().
  observed <- multtest::mt.teststat.num.denum(
    mat[complete, , drop = FALSE],
    cl,
    test = "t"
  )
  se <- observed$teststat.denum
  usable <- is.finite(observed$teststat.num) & is.finite(se) & round(se, 10) > 0
  rows <- complete[usable]
  if (length(rows) < 2L) {
    cli::cli_abort(
      "SAM requires at least two features with estimable within-group variation. Found {length(rows)} of {nrow(mat)}."
    )
  }
  if (is.null(s0) && length(unique(se[usable])) < 25L) {
    cli::cli_abort(c(
      "Automatic {.arg s0} estimation requires at least 25 distinct testable standard errors.",
      "i" = "Supply a prespecified nonnegative {.arg s0} to analyze a smaller panel."
    ))
  }
  exclusion <- rep("non_finite", nrow(mat))
  exclusion[reserved] <- "backend_missing_code"
  exclusion[complete] <- "unestimable_variance"
  exclusion[rows] <- NA_character_
  if (length(rows) < nrow(mat)) {
    cli::cli_warn(
      c(
        "SAM inference uses {length(rows)} of {nrow(mat)} features.",
        "i" = "All input IDs are retained. See {.field exclusion} for the reason each feature was omitted."
      ),
      class = "glystats_sam_exclusions"
    )
  }
  negative <- observed$teststat.num[usable] < 0
  list(
    rows = rows,
    n_finite = n_finite,
    exclusion = exclusion,
    singleton_tail = sum(negative) == 1L || sum(!negative) == 1L
  )
}
