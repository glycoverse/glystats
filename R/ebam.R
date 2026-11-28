#' Compare two independent groups by empirical Bayes analysis of microarrays
#'
#' Estimate each feature's posterior probability of differential abundance
#' using EBAM. Positive scores indicate a higher mean in the test group on
#' the analysis scale.
#'
#' @inheritParams gly_limma
#' @inheritParams gly_scale
#' @inheritParams gly_sam
#' @param posterior_threshold One posterior probability in `(0, 1]` (default
#'   0.9), or a vector of distinct probabilities for comparing thresholds.
#'   For one threshold, `called` identifies features in the selected score
#'   tails. For a vector, `called` is `NA` and `thresholds` reports each choice.
#' @param a0 Positive regularization constant in analysis units. The default
#'   uses the median Welch standard error of the tested features
#'   (`quan.a0 = 0.5` in siggenes). A positive constant keeps scores defined
#'   when a permutation has zero within-group variance.
#'
#' @details
#' Supply two independent groups with at least two samples each. Label
#' permutations assume exchangeability and preserve the group sizes. Use
#' [gly_limma()] for paired, blocked, or covariate-adjusted comparisons.
#' Install the optional backend with `BiocManager::install("siggenes")`.
#'
#' [siggenes::ebam()] divides the test-minus-reference mean difference by
#' the Welch standard error plus `a0`. It compares observed and permuted
#' score densities to estimate posterior probabilities across the eligible
#' features. Raw abundances use `log2(x + 1e-6)` by default. On relative
#' abundance and `logx` scales the direction refers to the supplied
#' normalization reference or transformed coordinate, respectively.
#'
#' Features enter the fit when every analysis value is finite and the
#' within-group standard error is estimable and positive after rounding to
#' ten decimal places. The backend also reserves the finite value
#' `-93074815` to indicate missing data. Each exclusion is reported and
#' keeps its original ID. The posterior estimates and FDR apply to the
#' retained family. Fitting the score density requires variation across
#' features.
#'
#' The fit uses the full tail-count calculation (`fast = FALSE`). The
#' permutation budget passed to siggenes is capped at the number of possible
#' group allocations. For example, three samples per group give 20 allocations
#' even if `n_perm = 100`. This makes the denominator used for the null tail
#' counts equal the number of allocations actually evaluated.
#'
#' @returns A `glystats_ebam_res` object inheriting from `glystats_res`.
#'   [get_tidy_result()] returns one row per input feature in assay order.
#'
#'   * `.variable` gives the SE feature ID (`variable` for legacy experiments).
#'     `ref_group` and `test_group` identify the comparison.
#'   * `z` is the signed EBAM score, `posterior` is the estimated probability
#'     of differential abundance, and `local_fdr = 1 - posterior` is the
#'     estimated posterior probability of the null for that feature.
#'   * For a single threshold, `called` uses the inclusive lower and upper
#'     score cutoffs from siggenes. The backend follows each tail inward from
#'     its extreme. An isolated interior score with a high posterior can
#'     fall outside both tails.
#'   * `n_finite` counts finite analysis values. Excluded features have
#'     missing statistics and calls, with reason `"non_finite"`,
#'     `"backend_missing_code"`, or `"unestimable_variance"` in `exclusion`.
#'
#'   `thresholds` gives the posterior threshold, score cutoffs, number called
#'   (`n_called`), estimated null proportion (`p0`), and permutation FDR
#'   (`fdr`) for those calls. FDR is `p0` times the mean number of permuted
#'   scores in the selected tails divided by `n_called`, capped at one.
#'   An empty tail has a cutoff of infinite magnitude. When neither tail
#'   calls a feature, siggenes reports FDR zero.
#'   If siggenes leaves an empty upper tail undefined or splits tied scores
#'   in its reported count, the table uses the inclusive count and the
#'   fitted mean false calls at those same cutoffs to calculate FDR.
#'
#'   `settings$n_perm` records the request. `diagnostics` gives input and
#'   tested feature counts, fitted `a0`, actual permutation rows
#'   (`n_permutations`), and `n_distinct_permutations`. [get_raw_result()]
#'   returns the fitted `EBAM` object with original one-based assay row
#'   positions as its feature keys. The result also stores the input
#'   `meta_data`, `preprocessing` settings, and SE `identifiers` mapping.
#'
#' @references
#' Efron B, Tibshirani R, Storey JD, Tusher V (2001). Empirical Bayes
#' Analysis of a Microarray Experiment. *Journal of the American Statistical
#' Association*, 96, 1151-1160. \doi{10.1198/016214501753382129}.
#'
#' [siggenes on Bioconductor](https://bioconductor.org/packages/siggenes).
#'
#' @seealso [gly_sam()], [gly_identity], [gly_scale]
#' @examplesIf requireNamespace("siggenes", quietly = TRUE)
#' exp <- glyexp::real_experiment2
#' two_groups <- exp[, exp$group %in% c("H", "C")]
#' fit <- gly_ebam(two_groups, ref_group = "H", posterior_threshold = 0.9, seed = 2)
#' fit$thresholds
#' subset(get_tidy_result(fit), called, select = c(.variable, z, posterior))
#' @export
gly_ebam <- function(
  exp,
  group_col = "group",
  ref_group = NULL,
  posterior_threshold = 0.9,
  n_perm = 100L,
  a0 = NULL,
  seed = 1L,
  add_info = TRUE,
  log_transform = NULL
) {
  .assert_data_container(exp)
  checkmate::assert_string(group_col)
  checkmate::assert_flag(add_info)
  checkmate::assert_int(seed, lower = 0L)
  checkmate::assert_int(n_perm, lower = 2L)
  checkmate::assert_vector(posterior_threshold, strict = TRUE)
  checkmate::assert_numeric(
    posterior_threshold,
    min.len = 1L,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE,
    unique = TRUE
  )
  if (any(posterior_threshold <= 0)) {
    cli::cli_abort("{.arg posterior_threshold} must be in (0, 1].")
  }
  checkmate::assert_number(a0, lower = 0, finite = TRUE, null.ok = TRUE)
  if (isTRUE(a0 == 0)) {
    cli::cli_abort("{.arg a0} must be positive to regularize all permutations.")
  }
  withr::local_seed(seed)
  rlang::check_installed("siggenes", reason = "for two-group EBAM")

  meta_data <- .get_meta_data(exp)
  preprocessing <- .resolve_transformation(
    exp,
    log_transform,
    meta_data = meta_data
  )
  mat <- .get_expr_mat(exp)
  groups <- .siggenes_groups(exp, mat, group_col, ref_group, "EBAM")
  cl <- as.integer(groups) - 1L
  mat <- .log_transform_expr_mat(mat, preprocessing$log_transform)
  features <- .siggenes_testable_features(mat, cl, "EBAM")
  rows <- features$rows
  # siggenes truncates feature names to 50 characters in its fitted slots.
  keys <- as.character(seq_len(nrow(mat)))
  # z.ebam divides the tail counts by B even when it enumerates fewer labels.
  draws <- min(n_perm, round(choose(length(cl), sum(cl == 0L))))
  fit <- tryCatch(
    siggenes::ebam(
      mat[rows, , drop = FALSE],
      cl,
      gene.names = keys[rows],
      delta = posterior_threshold,
      a0 = a0,
      quan.a0 = if (is.null(a0)) 0.5,
      B = draws,
      B.more = 0,
      fast = FALSE
    ),
    error = function(err) {
      cli::cli_abort(
        "siggenes could not fit the EBAM score density to the {length(rows)} eligible features.",
        parent = err
      )
    }
  )
  actual <- nrow(fit@mat.samp)
  if (actual != draws) {
    cli::cli_abort(
      "EBAM evaluated {actual} permutations but used a denominator of {draws} for its null tail counts."
    )
  }
  distinct <- nrow(unique(fit@mat.samp))
  if (distinct < 2L) {
    cli::cli_abort(
      "EBAM drew fewer than two distinct label permutations. Increase {.arg n_perm}."
    )
  }

  z <- unname(fit@z[keys])
  posterior <- unname(fit@posterior[keys])
  if (anyNA(z[rows]) || anyNA(posterior[rows])) {
    cli::cli_abort(
      "EBAM did not return a score and posterior for every eligible feature key."
    )
  }
  thresholds <- .ebam_thresholds(fit)
  called <- rep(NA, nrow(mat))
  if (length(posterior_threshold) == 1L) {
    called <- z <= thresholds$cutoff_lower | z >= thresholds$cutoff_upper
  }
  tidy <- tibble::tibble(
    variable = rownames(mat),
    ref_group = levels(groups)[1L],
    test_group = levels(groups)[2L],
    z = z,
    posterior = posterior,
    local_fdr = 1 - posterior,
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
      a0 = unname(fit@a0),
      n_permutations = actual,
      n_distinct_permutations = distinct
    ),
    settings = list(
      group_col = group_col,
      ref_group = levels(groups)[1L],
      test_group = levels(groups)[2L],
      posterior_threshold = posterior_threshold,
      n_perm = n_perm,
      a0 = a0,
      seed = seed
    ),
    meta_data = meta_data,
    preprocessing = preprocessing
  )
  class(result) <- c("glystats_ebam_res", "glystats_res")
  .record_se_identity(result, exp)
}

.ebam_thresholds <- function(fit) {
  grid <- as.data.frame(fit@mat.fdr)
  # compNumber() can leave the upper cutoff undefined for all-negative scores.
  if (all(fit@z < 0)) {
    grid$CU[is.na(grid$CU)] <- Inf
  }
  # summary(EBAM) selects inclusive tails, including all ties at a cutoff.
  counts <- vapply(
    seq_len(nrow(grid)),
    function(i) {
      sum(fit@z <= grid$CL[i] | fit@z >= grid$CU[i])
    },
    integer(1)
  )
  if (anyNA(counts)) {
    cli::cli_abort("EBAM returned an undefined score cutoff.")
  }
  recount <- is.na(grid$FDR) | is.na(grid$Number) | counts != grid$Number
  if (any(recount)) {
    # The fitted tail vectors already contain the mean counts per permutation.
    lower <- match(grid$CL[recount], c(-Inf, fit@z))
    upper <- match(grid$CU[recount], c(Inf, fit@z))
    expected_false <- c(0, fit@vec.neg)[lower] + c(0, fit@vec.pos)[upper]
    grid$FDR[recount] <- pmin(
      1,
      fit@p0 * expected_false / pmax(1, counts[recount])
    )
  }
  if (anyNA(grid$FDR)) {
    cli::cli_abort(
      "The EBAM score cutoffs and permutation FDR do not define a consistent set of calls."
    )
  }
  tibble::tibble(
    posterior_threshold = grid$Delta,
    cutoff_lower = grid$CL,
    cutoff_upper = grid$CU,
    n_called = counts,
    p0 = fit@p0,
    fdr = grid$FDR
  )
}
