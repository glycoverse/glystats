ebam_test_se <- function(n = 400L, n_ref = 3L, n_test = 3L, seed = 613L) {
  withr::local_seed(seed)
  group <- factor(
    rep(c("H", "C"), c(n_ref, n_test)),
    levels = c("H", "C", "unused")
  )
  mat <- matrix(rnorm(n * length(group)), nrow = n)
  up <- seq_len(n %/% 10L)
  down <- length(up) + up
  mat[up, group == "C"] <- mat[up, group == "C"] + 4
  mat[down, group == "C"] <- mat[down, group == "C"] - 4
  dimnames(mat) <- list(paste0("g", seq_len(n)), paste0("s", seq_along(group)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = mat),
    colData = S4Vectors::DataFrame(group = group, sample = rev(colnames(mat))),
    rowData = S4Vectors::DataFrame(
      variable = paste0("annotation ", seq_len(n))
    ),
    metadata = list(scale = "log2")
  )
}

test_that("EBAM normalizes the null counts by all 20 allocations in a 3+3 design", {
  skip_if_not_installed("siggenes")
  se <- ebam_test_se()
  mat <- SummarizedExperiment::assay(se)
  cl <- as.integer(se$group == "C")
  thresholds <- c(0.9, 0.5, 1)
  direct <- siggenes::ebam(
    mat,
    cl,
    delta = thresholds,
    B = 20L,
    B.more = 0,
    fast = FALSE,
    quan.a0 = 0.5,
    gene.names = as.character(seq_len(nrow(mat)))
  )
  exact <- gly_ebam(se, posterior_threshold = thresholds, n_perm = 20L)
  capped <- gly_ebam(se, posterior_threshold = thresholds, n_perm = 100L)
  expect_s3_class(capped, c("glystats_ebam_res", "glystats_res"))
  expect_equal(get_raw_result(exact), direct)
  expect_equal(get_raw_result(capped), direct)
  expect_identical(capped$settings$n_perm, 100L)
  expect_identical(capped$diagnostics$n_permutations, 20L)
  expect_identical(capped$diagnostics$n_distinct_permutations, 20L)
  expect_identical(capped$thresholds, exact$thresholds)
  expect_equal(capped$diagnostics$a0, unname(direct@a0))
  tidy <- get_tidy_result(capped)
  expect_equal(tidy$z, unname(direct@z))
  expect_equal(tidy$posterior, unname(direct@posterior))
  expect_equal(tidy$local_fdr, unname(direct@local))
  expect_false(any(c("p_val", "p_adj", "q_value", "fdr") %in% names(tidy)))
  expect_true(all(is.na(tidy$called)))
  expect_identical(capped$thresholds$posterior_threshold, thresholds)
  expect_gt(capped$thresholds$n_called[1], 0L)
  expect_gt(capped$thresholds$fdr[1], 0)
  expect_equal(capped$thresholds$n_called[3], 0L)
  expect_equal(capped$thresholds$fdr[3], 0)

  # Enumerate labels independently and count every permuted score in each tail.
  allocations <- utils::combn(ncol(mat), sum(cl), FUN = function(test) {
    as.integer(seq_len(ncol(mat)) %in% test)
  })
  expect_setequal(
    apply(allocations, 2, paste, collapse = ""),
    apply(direct@mat.samp, 1, paste, collapse = "")
  )
  permuted <- as.vector(apply(allocations, 2L, function(labels) {
    parts <- multtest::mt.teststat.num.denum(mat, labels, test = "t")
    parts$teststat.num / (parts$teststat.denum + direct@a0)
  }))
  pos <- vapply(abs(tidy$z), function(z) sum(permuted >= z), integer(1))
  neg <- vapply(abs(tidy$z), function(z) sum(permuted <= -z), integer(1))
  expect_equal(unname(direct@vec.pos) * ncol(allocations), pos)
  expect_equal(unname(direct@vec.neg) * ncol(allocations), neg)
  grid <- capped$thresholds
  for (i in seq_len(nrow(grid))) {
    in_tail <- function(z) z <= grid$cutoff_lower[i] | z >= grid$cutoff_upper[i]
    n_called <- sum(in_tail(tidy$z))
    if (grid$posterior_threshold[i] == 0.5) {
      expect_true(any(tidy$z[1:40] >= grid$cutoff_upper[i]))
      expect_true(any(tidy$z[41:80] <= grid$cutoff_lower[i]))
    }
    false_calls <- sum(in_tail(permuted)) / ncol(allocations)
    expect_equal(grid$n_called[i], n_called)
    expect_equal(
      grid$fdr[i],
      min(1, direct@p0 * false_calls / max(1, n_called))
    )
  }
})

test_that("EBAM calls and FDR use the backend's inclusive score tails", {
  skip_if_not_installed("siggenes")
  se <- ebam_test_se(seed = 2L)
  res <- gly_ebam(se)
  tidy <- get_tidy_result(res)
  direct <- siggenes::summary(get_raw_result(res), delta = 0.9, what = "stats")
  expect_identical(which(tidy$called), as.integer(direct@row.sig.genes))
  expect_equal(sum(tidy$called), unname(direct@mat.fdr[, "Number"]))
  expect_equal(res$thresholds$fdr, unname(direct@mat.fdr[, "FDR"]))
  expect_true(any(tidy$posterior >= 0.9 & !tidy$called))
  expect_true(any(tidy$called[1:40] & tidy$z[1:40] > 0))
  expect_true(all(tidy$ref_group == "H" & tidy$test_group == "C"))
  reversed <- get_tidy_result(gly_ebam(se, ref_group = "C"))
  expect_equal(reversed$z, -tidy$z)
  expect_true(all(reversed$ref_group == "C" & reversed$test_group == "H"))
})

test_that("an empty upper tail and tied scores retain their inclusive calls", {
  skip_if_not_installed("siggenes")
  se <- ebam_test_se(n = 100L)
  mat <- SummarizedExperiment::assay(se)
  mat[, 4:6] <- mat[, 4:6] - 15
  mat[1, ] <- mat[2, ] <- c(-1, 0, 1, -61, -60, -59)
  SummarizedExperiment::assay(se) <- mat
  # The separated data trigger the backend's fallback estimate of p0.
  expect_warning(fit <- gly_ebam(se, posterior_threshold = 1), "pi0")
  raw <- get_raw_result(fit)
  called <- which(get_tidy_result(fit)$called)
  expect_identical(called, 1:2)
  expect_identical(
    called,
    as.integer(siggenes::summary(raw, delta = 1, what = "stats")@row.sig.genes)
  )
  expect_equal(fit$thresholds$cutoff_upper, Inf)
  expect_equal(fit$thresholds$n_called, length(called))
  expect_equal(
    fit$thresholds$fdr,
    unname(raw@p0 * raw@vec.neg[1] / length(called))
  )
})

test_that("EBAM retains exclusions and long, colliding feature prefixes", {
  skip_if_not_installed("siggenes")
  se <- ebam_test_se()
  rownames(se)[c(1, 41)] <- paste0(strrep("glycan", 10), c("-up", "-down"))
  SummarizedExperiment::rowData(se)$posterior <- seq_len(nrow(se))
  mat <- SummarizedExperiment::assay(se)
  mat[2, ] <- NA_real_
  mat[4, 1:2] <- c(Inf, NaN)
  mat[6, ] <- 1
  mat[8, ] <- as.integer(se$group == "C")
  mat[10, 1] <- -93074815
  SummarizedExperiment::assay(se) <- mat
  bad <- c(2L, 4L, 6L, 8L, 10L)
  expect_warning(
    res <- gly_ebam(se, a0 = 0.7),
    class = "glystats_ebam_exclusions"
  )
  tidy <- get_tidy_result(res)
  expect_identical(tidy$.variable, rownames(se))
  expect_identical(tidy$.variable, res$identifiers$variable$data$id)
  expect_identical(tidy$variable, paste0("annotation ", seq_len(nrow(se))))
  expect_identical(tidy$posterior_metadata, seq_len(nrow(se)))
  expect_identical(tidy$n_finite, as.integer(rowSums(is.finite(mat))))
  expect_identical(
    tidy$exclusion[bad],
    c(
      "non_finite",
      "non_finite",
      "unestimable_variance",
      "unestimable_variance",
      "backend_missing_code"
    )
  )
  expect_true(all(is.na(tidy[bad, c("z", "posterior", "local_fdr", "called")])))
  expect_true(all(is.na(tidy$exclusion[-bad])))
  expect_gt(tidy$z[1], 0)
  expect_lt(tidy$z[41], 0)
  expect_equal(res$diagnostics$n_input, nrow(se))
  expect_equal(res$diagnostics$n_tested, nrow(se) - length(bad))
  direct <- siggenes::ebam(
    mat[-bad, , drop = FALSE],
    as.integer(se$group == "C"),
    gene.names = as.character(seq_len(nrow(se))[-bad]),
    a0 = 0.7,
    B = 20,
    B.more = 0,
    fast = FALSE
  )
  expect_equal(get_raw_result(res), direct)
  expect_equal(tidy$z[-bad], unname(direct@z))
  expect_equal(tidy$posterior[-bad], unname(direct@posterior))
  subset <- se[-bad, rev(seq_len(ncol(se)))]
  dimnames(subset) <- list(NULL, NULL)
  no_info <- gly_ebam(subset, a0 = 0.7, add_info = FALSE)
  expect_true(no_info$identifiers$variable$generated)
  expect_true(no_info$identifiers$sample$generated)
  expect_identical(
    get_tidy_result(no_info)$.variable,
    no_info$identifiers$variable$data$id
  )
  expect_false("variable" %in% names(get_tidy_result(no_info)))
  expect_equal(get_tidy_result(no_info)$z, tidy$z[-bad])
})

test_that("sampled EBAM fits reproduce the direct backend and restore RNG", {
  skip_if_not_installed("siggenes")
  se <- ebam_test_se(n = 160L, n_ref = 5L, n_test = 7L)
  se <- se[, order(sequence(c(5L, 7L)))]
  withr::local_seed(23L)
  before <- .Random.seed
  res <- gly_ebam(se, n_perm = 30L, seed = 8L, a0 = 0.8)
  expect_identical(.Random.seed, before)
  direct <- withr::with_seed(
    8L,
    siggenes::ebam(
      SummarizedExperiment::assay(se),
      as.integer(se$group == "C"),
      gene.names = as.character(seq_len(nrow(se))),
      a0 = 0.8,
      B = 30L,
      B.more = 0,
      fast = FALSE
    )
  )
  expect_equal(get_raw_result(res), direct)
  repeat_fit <- gly_ebam(se, n_perm = 30L, seed = 8L, a0 = 0.8)
  expect_identical(get_raw_result(res), get_raw_result(repeat_fit))
  other <- gly_ebam(se, n_perm = 30L, seed = 9L, a0 = 0.8)
  expect_false(identical(get_raw_result(other)@mat.samp, direct@mat.samp))
  expect_identical(.Random.seed, before)
  rm(".Random.seed", envir = .GlobalEnv)
  gly_ebam(se, n_perm = 30L, seed = 8L, a0 = 0.8)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("EBAM handles extreme requests and large allocation spaces", {
  skip_if_not_installed("siggenes")
  small <- ebam_test_se(n = 160L)
  capped <- gly_ebam(small, n_perm = .Machine$integer.max)
  expect_identical(capped$settings$n_perm, .Machine$integer.max)
  expect_identical(capped$diagnostics$n_permutations, 20L)
  almost_all <- gly_ebam(small, n_perm = 19L)
  expect_identical(almost_all$diagnostics$n_distinct_permutations, 19L)
  # choose(1030, 515) overflows double precision. Only five rows are drawn.
  large <- ebam_test_se(n = 100L, n_ref = 515L, n_test = 515L)
  sampled <- gly_ebam(large, n_perm = 5L, seed = 4L)
  perm <- get_raw_result(sampled)@mat.samp
  expect_identical(dim(perm), c(5L, 1030L))
  expect_identical(
    sampled$diagnostics$n_distinct_permutations,
    nrow(unique(perm))
  )
  expect_true(all(rowSums(perm) == 515L))
  repeated <- gly_ebam(
    ebam_test_se(n = 100L, n_ref = 9L, n_test = 9L),
    n_perm = 600L,
    seed = 19L
  )
  perm <- get_raw_result(repeated)@mat.samp
  expect_identical(repeated$diagnostics$n_permutations, 600L)
  expect_lt(repeated$diagnostics$n_distinct_permutations, 600L)
  expect_identical(
    repeated$diagnostics$n_distinct_permutations,
    nrow(unique(perm))
  )
})

test_that("EBAM validates group labels, supported designs and tuning arguments", {
  skip_if_not_installed("siggenes")
  se <- ebam_test_se(n = 80L)
  se$three <- rep(letters[1:3], 2)
  se$bad <- c("C", "C", "H", "H", "H", NA)
  expect_error(gly_ebam(se, "three"), "exactly 2")
  expect_error(gly_ebam(se, "bad"), "group label")
  se$bad[6] <- " "
  expect_error(gly_ebam(se, "bad"), "group label")
  expect_error(gly_ebam(se, "absent"), "not found")
  expect_error(gly_ebam(se, ref_group = "missing"), "ref_group")
  expect_error(gly_ebam(se[, 1:3]), "exactly 2")
  expect_error(gly_ebam(se[, 1:4]), "at least two samples")
  expect_error(gly_ebam(se, paired = TRUE), "unused argument")
  expect_error(gly_ebam(se, subject_col = "id"), "unused argument")
  expect_error(gly_ebam(se, covariate_cols = "age"), "unused argument")
  for (p in list(
    NULL,
    numeric(),
    0,
    -1,
    1.01,
    NA_real_,
    Inf,
    c(1, 1),
    matrix(0.9)
  )) {
    expect_error(gly_ebam(se, posterior_threshold = p), "posterior_threshold")
  }
  for (b in c(0, 1, 1.5, NA, Inf)) {
    expect_error(gly_ebam(se, n_perm = b), "n_perm")
  }
  for (a0 in c(-1, 0, Inf)) {
    expect_error(gly_ebam(se, a0 = a0), "a0")
  }
  expect_error(gly_ebam(se, seed = -1L), "seed")
  expect_error(gly_ebam(se[1, ]), "at least two complete")
  mat <- SummarizedExperiment::assay(se)
  SummarizedExperiment::assay(se)[] <- NA_real_
  expect_error(gly_ebam(se), "Found 0")
  SummarizedExperiment::assay(se)[] <- 1
  expect_error(gly_ebam(se), "estimable within-group variation")
  mat[, 4:6] <- mat[, 1:3]
  SummarizedExperiment::assay(se) <- mat
  withr::local_seed(42L)
  before <- .Random.seed
  expect_error(gly_ebam(se), "could not fit the EBAM score density")
  expect_identical(.Random.seed, before)
})

test_that("EBAM uses real GlycomicSE identifiers, annotations and analysis scales", {
  skip_if_not_installed("siggenes")
  se <- glyexp::real_experiment2
  se <- se[, se$group %in% c("H", "C")]
  expect_s4_class(se, "GlycomicSE")
  expect_warning(
    res <- gly_ebam(se, ref_group = "H", seed = 2L),
    class = "glystats_ebam_exclusions"
  )
  tidy <- get_tidy_result(res)
  expect_identical(tidy$.variable, rownames(se))
  expect_identical(
    tidy$glycan_composition,
    SummarizedExperiment::rowData(se)$glycan_composition
  )
  expect_identical(res$meta_data, S4Vectors::metadata(se))
  expect_identical(res$preprocessing$analysis_scale, "log2")
  expect_equal(sum(tidy$called, na.rm = TRUE), res$thresholds$n_called)
  logged <- se
  SummarizedExperiment::assay(logged) <- log2(
    SummarizedExperiment::assay(se) + 1e-6
  )
  S4Vectors::metadata(logged)$scale <- "logx"
  expect_warning(
    native <- gly_ebam(logged, ref_group = "H", seed = 2L),
    class = "glystats_ebam_exclusions"
  )
  expect_identical(native$preprocessing$analysis_scale, "logx")
  expect_false(native$preprocessing$log_transform)
  expect_equal(get_tidy_result(native)$z, tidy$z)
  expect_equal(get_tidy_result(native)$posterior, tidy$posterior)
  expect_equal(native$thresholds, res$thresholds)
})
