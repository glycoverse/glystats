sam_test_se <- function(n_features = 60L, n_ref = 5L, n_test = 7L) {
  withr::local_seed(37)
  group <- factor(
    rep(c("H", "C"), c(n_ref, n_test)),
    levels = c("H", "C", "unused")
  )
  values <- matrix(rnorm(n_features * length(group)), n_features)
  values[1, group == "C"] <- values[1, group == "C"] + 8
  values[2, group == "C"] <- values[2, group == "C"] - 8
  dimnames(values) <- list(
    paste0("g", seq_len(n_features)),
    paste0("s", seq_along(group))
  )
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(
      group = group,
      sample = rev(colnames(values))
    ),
    rowData = S4Vectors::DataFrame(
      annotation = paste0("annotated ", rownames(values))
    ),
    metadata = list(scale = "log2")
  )
  exp[, order(sequence(c(n_ref, n_test)))]
}

test_that("SAM reproduces the direct fit, feature statistics and threshold grid", {
  skip_if_not_installed("siggenes")
  se <- sam_test_se()
  deltas <- c(0.2, 1, 100)
  res <- gly_sam(se, delta = deltas, n_perm = 80L, seed = 19L)
  direct <- withr::with_seed(
    19L,
    siggenes::sam(
      SummarizedExperiment::assay(se),
      as.integer(se$group == "C"),
      gene.names = as.character(seq_len(nrow(se))),
      control = siggenes::samControl(delta = deltas),
      B = 80,
      B.more = 0,
      R.fold = 0,
      na.replace = FALSE
    )
  )
  expect_s3_class(res, c("glystats_sam_res", "glystats_res"))
  expect_equal(get_raw_result(res), direct)
  tidy <- get_tidy_result(res)
  expect_identical(tidy$.variable, rownames(se))
  expect_equal(tidy$d, unname(direct@d))
  expect_equal(tidy$p_val, unname(direct@p.value))
  expect_equal(tidy$q_value, unname(direct@q.value))
  expect_false(any(c("p_adj", "fdr") %in% names(tidy)))
  expect_true(all(is.na(tidy$called)))
  expect_true(all(is.na(tidy$exclusion)))
  expect_equal(
    unname(as.matrix(res$thresholds)),
    unname(direct@mat.fdr[, c(
      "Delta",
      "cutlow",
      "cutup",
      "Called",
      "False",
      "p0",
      "FDR"
    )])
  )
  expect_equal(res$diagnostics$s0, unname(direct@s0))
  expect_identical(res$diagnostics$n_permutations, nrow(direct@mat.samp))
  expect_identical(
    res$diagnostics$n_distinct_permutations,
    nrow(unique(direct@mat.samp))
  )

  auto <- gly_sam(se, n_perm = 80L, seed = 19L)
  expect_equal(nrow(auto$thresholds), 10L)
  expect_equal(auto$thresholds$delta, get_raw_result(auto)@mat.fdr[, "Delta"])
  expect_true(all(is.na(get_tidy_result(auto)$called)))

  selected <- gly_sam(se, delta = 1, n_perm = 80L, seed = 19L)
  selected_tidy <- get_tidy_result(selected)
  direct_calls <- siggenes::summary(
    direct,
    delta = 1,
    what = "stats"
  )@row.sig.genes
  expect_identical(which(selected_tidy$called), as.integer(direct_calls))
  expect_true(all(selected_tidy$called[1:2]))
  expect_gt(selected_tidy$d[1], 0)
  expect_lt(selected_tidy$d[2], 0)
  expect_equal(sum(selected_tidy$called), selected$thresholds$n_called)
  empty <- gly_sam(se, delta = 100, n_perm = 80L, seed = 19L)
  expect_false(any(get_tidy_result(empty)$called))
  expect_equal(empty$thresholds$n_called, 0)
  expect_equal(unname(empty$thresholds$fdr), 0)
})

test_that("small panels retain long SE IDs, annotations and feature order", {
  skip_if_not_installed("siggenes")
  se <- sam_test_se(12, 3, 3)
  # Two structures from glyexp::real_experiment2 share their first 50 characters.
  ids <- c(
    "GlcNAc(??-?)Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-",
    "GlcNAc(??-?)Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)[Fuc(??-?)]GlcNAc(??-"
  )
  expect_identical(substr(ids[1], 1, 50), substr(ids[2], 1, 50))
  rownames(se)[1:2] <- ids
  SummarizedExperiment::rowData(se)$d <- seq_len(nrow(se))
  direct <- siggenes::sam(
    SummarizedExperiment::assay(se),
    as.integer(se$group == "C"),
    control = siggenes::samControl(delta = 1),
    s0 = 0.5,
    B = 20L,
    B.more = 0,
    R.fold = 0,
    na.replace = FALSE
  )
  res <- gly_sam(se, delta = 1, s0 = 0.5, n_perm = 100L)
  tidy <- get_tidy_result(res)
  expect_identical(tidy$.variable, rownames(se))
  expect_identical(tidy$annotation, paste0("annotated g", seq_len(nrow(se))))
  expect_identical(tidy$d_metadata, seq_len(nrow(se)))
  expect_equal(tidy$d, unname(direct@d))
  expect_equal(tidy$p_val, unname(direct@p.value))
  expect_equal(tidy$q_value, unname(direct@q.value))
  expect_identical(
    names(get_raw_result(res)@d),
    as.character(seq_len(nrow(se)))
  )
  expect_true(tidy$d[1] > 0 && tidy$d[2] < 0)
  expect_identical(res$settings$n_perm, 100L)
  expect_identical(res$diagnostics$n_permutations, 20L)
  expect_identical(res$diagnostics$n_distinct_permutations, 20L)
  expect_true(all(tidy$ref_group == "H" & tidy$test_group == "C"))
  reverse <- gly_sam(se, ref_group = "C", delta = 1, s0 = 0.5)
  expect_equal(get_tidy_result(reverse)$d, -tidy$d)
  expect_identical(get_tidy_result(reverse)$called, tidy$called)
  expect_true(all(get_tidy_result(reverse)$test_group == "H"))
  reordered <- gly_sam(se[c(2, 1, 3:12), 6:1], delta = 1, s0 = 0.5)
  original <- match(tidy$.variable, get_tidy_result(reordered)$.variable)
  expect_equal(get_tidy_result(reordered)$d[original], tidy$d)
  expect_identical(
    get_tidy_result(reordered)$annotation[original],
    tidy$annotation
  )
})

test_that("SAM reports the tested family and restores excluded positions", {
  skip_if_not_installed("siggenes")
  se <- sam_test_se()
  dimnames(se) <- list(rep("repeated", nrow(se)), NULL)
  mat <- SummarizedExperiment::assay(se)
  mat[2, ] <- NA_real_
  mat[4, 1:2] <- c(NA_real_, Inf)
  mat[6, ] <- 1
  mat[8, ] <- as.numeric(se$group == "C")
  mat[10, ] <- seq_len(ncol(se)) * 1e-12
  mat[12, 1] <- -93074815
  SummarizedExperiment::assay(se) <- mat
  bad <- c(2L, 4L, 6L, 8L, 10L, 12L)
  expect_warning(
    res <- gly_sam(se, delta = 1),
    class = "glystats_sam_exclusions"
  )
  tidy <- get_tidy_result(res)
  expect_identical(
    tidy$.variable,
    paste0(".glystats_variable_", seq_len(nrow(se)))
  )
  expect_identical(tidy$.variable, res$identifiers$variable$data$id)
  expect_true(res$identifiers$sample$generated)
  expect_identical(tidy$n_finite, as.integer(rowSums(is.finite(mat))))
  expect_identical(
    tidy$exclusion[bad],
    c(
      rep("non_finite", 2),
      rep("unestimable_variance", 3),
      "backend_missing_code"
    )
  )
  expect_true(all(is.na(tidy[bad, c("d", "p_val", "q_value", "called")])))
  expect_true(all(is.na(tidy$exclusion[-bad])))
  expect_identical(tidy$annotation, paste0("annotated g", seq_len(nrow(se))))
  expect_equal(res$diagnostics$n_tested, nrow(se) - length(bad))
  expect_identical(
    names(get_raw_result(res)@d),
    as.character(seq_len(nrow(se))[-bad])
  )
  complete <- gly_sam(se[-bad, ], delta = 1, add_info = FALSE)
  expect_equal(get_tidy_result(complete)$d, tidy$d[-bad])
  expect_equal(get_tidy_result(complete)$q_value, tidy$q_value[-bad])
  expect_equal(complete$thresholds, res$thresholds)
  expect_false("annotation" %in% names(get_tidy_result(complete)))
})

test_that("SAM preserves the backend's false-call count for a singleton tail", {
  skip_if_not_installed("siggenes")
  se <- sam_test_se(12, 3, 3)
  mat <- SummarizedExperiment::assay(se)
  test <- se$group == "C"
  difference <- rowMeans(mat[, test]) - rowMeans(mat[, !test])
  mat <- mat * sign(difference)
  mat[2, ] <- -mat[2, ]
  SummarizedExperiment::assay(se) <- mat
  direct <- siggenes::sam(
    mat,
    as.integer(test),
    control = siggenes::samControl(delta = c(0.05, 1)),
    R.fold = 0,
    na.replace = FALSE,
    s0 = 0.5,
    B = 20,
    B.more = 0,
    n.subset = 20
  )
  res <- gly_sam(se, delta = c(0.05, 1), s0 = 0.5)
  expect_equal(sum(get_tidy_result(res)$d < 0), 1L)
  expect_equal(res$thresholds$expected_false, direct@mat.fdr[, "False"])
  expect_equal(res$thresholds$fdr, direct@mat.fdr[, "FDR"])
})

test_that("SAM reproduces sampled permutations and restores the caller's RNG", {
  skip_if_not_installed("siggenes")
  se <- sam_test_se(60, 10, 10)
  withr::local_seed(27L)
  before <- .Random.seed
  a <- gly_sam(se, seed = 3L, n_perm = 32L)
  expect_identical(.Random.seed, before)
  b <- gly_sam(se, seed = 3L, n_perm = 32L)
  expect_identical(get_raw_result(a)@mat.samp, get_raw_result(b)@mat.samp)
  expect_identical(a$thresholds, b$thresholds)
  c <- gly_sam(se, seed = 4L, n_perm = 32L)
  expect_false(identical(
    get_raw_result(a)@mat.samp,
    get_raw_result(c)@mat.samp
  ))
  expect_identical(.Random.seed, before)
  expect_identical(
    a$diagnostics$n_distinct_permutations,
    nrow(unique(get_raw_result(a)@mat.samp))
  )
  rm(".Random.seed", envir = .GlobalEnv)
  gly_sam(se, seed = 3L, n_perm = 32L)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  # Above siggenes' enumeration limit it can sample an allocation twice.
  sampled <- gly_sam(sam_test_se(12, 9, 9), seed = 19L, n_perm = 600L, s0 = 1)
  perms <- get_raw_result(sampled)@mat.samp
  expect_lt(nrow(unique(perms)), nrow(perms))
  expect_identical(sampled$diagnostics$n_permutations, 600L)
  expect_identical(
    sampled$diagnostics$n_distinct_permutations,
    nrow(unique(perms))
  )
})

test_that("SAM validates independent groups, testable data and numerical settings", {
  skip_if_not_installed("siggenes")
  se <- sam_test_se(12, 3, 3)
  expect_error(gly_sam(se), "25 distinct")
  expect_error(gly_sam(se[1, ], s0 = 1), "at least two complete")
  missing <- se
  SummarizedExperiment::assay(missing)[] <- NA_real_
  expect_error(gly_sam(missing, s0 = 1), "Found 0")
  constant <- se
  SummarizedExperiment::assay(constant)[-1, ] <- 3
  expect_error(gly_sam(constant, s0 = 1), "within-group variation")
  expect_error(gly_sam(se[, se$group == "C"]), "exactly 2")
  expect_error(gly_sam(se[, c(1:2, 4, 6)]), "at least two samples")
  se$three <- c("A", "B", "C", "A", "B", "C")
  expect_error(gly_sam(se, "three"), "exactly 2")
  se$absent <- as.character(se$group)
  se$absent[2] <- NA_character_
  expect_error(gly_sam(se, "absent"), "group label")
  se$absent[se$absent == "C" | is.na(se$absent)] <- " "
  expect_error(gly_sam(se, "absent"), "group label")
  expect_error(gly_sam(se, "unknown"), "not found")
  expect_error(gly_sam(se, ref_group = "Y"), "ref_group")
  expect_error(gly_sam(se, paired = TRUE), "unused argument")
  expect_error(gly_sam(se, subject_col = "id"), "unused argument")
  expect_error(gly_sam(se, covariate_cols = "age"), "unused argument")
  for (delta in list(numeric(), 0, -1, NA_real_, Inf, c(1, 1), matrix(1))) {
    expect_error(gly_sam(se, delta = delta), "delta")
  }
  for (n_perm in c(0, 1, 2.5, NA, Inf)) {
    expect_error(gly_sam(se, n_perm = n_perm), "n_perm")
  }
  expect_error(gly_sam(se, seed = -1), "seed")
  expect_error(gly_sam(se, s0 = -1), "s0")
  expect_error(gly_sam(se, s0 = Inf), "s0")
})

test_that("an undefined automatic regularizer is reported and never silently replaced", {
  skip_if_not_installed("siggenes")
  se <- sam_test_se(40, 3, 3)
  # Identical means but different standard errors for every feature.
  SummarizedExperiment::assay(se)[] <- outer(1:40, rep(-1:1, each = 2))
  withr::local_seed(801)
  before <- .Random.seed
  expect_error(gly_sam(se), "prespecified.*s0")
  expect_identical(.Random.seed, before)
  tidy <- get_tidy_result(gly_sam(se, s0 = 1, delta = 1))
  expect_true(all(tidy$d == 0 & tidy$p_val == 1 & !tidy$called))
})

test_that("SAM accepts a real GlycomicSE and respects its abundance scale", {
  skip_if_not_installed("siggenes")
  se <- glyexp::real_experiment2
  se <- se[, se$group %in% c("H", "C")]
  expect_s4_class(se, "GlycomicSE")
  expect_warning(
    res <- gly_sam(se, ref_group = "H", delta = 1),
    class = "glystats_sam_exclusions"
  )
  tidy <- get_tidy_result(res)
  expect_identical(tidy$.variable, rownames(se))
  expect_identical(
    tidy$glycan_composition,
    SummarizedExperiment::rowData(se)$glycan_composition
  )
  expect_identical(res$meta_data, S4Vectors::metadata(se))
  expect_identical(res$preprocessing$analysis_scale, "log2")
  logged <- se
  SummarizedExperiment::assay(logged) <- log2(
    SummarizedExperiment::assay(se) + 1e-6
  )
  S4Vectors::metadata(logged)$scale <- "logx"
  expect_warning(
    native <- gly_sam(logged, ref_group = "H", delta = 1),
    class = "glystats_sam_exclusions"
  )
  expect_identical(native$preprocessing$analysis_scale, "logx")
  expect_false(native$preprocessing$log_transform)
  expect_equal(get_tidy_result(native)$d, tidy$d)
  expect_equal(native$thresholds, res$thresholds)
})

test_that("SAM aligns legacy sample metadata by its sample ID", {
  skip_if_not_installed("siggenes")
  withr::local_options(list(lifecycle_verbosity = "quiet"))
  se <- sam_test_se(12, 3, 3)
  mat <- SummarizedExperiment::assay(se)
  legacy <- suppressWarnings(glyexp::experiment(
    mat,
    tibble::tibble(sample = colnames(mat), label = se$group),
    tibble::tibble(variable = rownames(mat)),
    "others"
  ))
  legacy$sample_info <- legacy$sample_info[rev(seq_len(ncol(mat))), ]
  a <- gly_sam(legacy, group_col = "label", log_transform = FALSE, s0 = 0.5)
  b <- gly_sam(se, s0 = 0.5)
  expect_identical(get_tidy_result(a)$variable, get_tidy_result(b)$.variable)
  expect_equal(get_tidy_result(a)$d, get_tidy_result(b)$d)
  legacy$sample_info$sample[1] <- "unknown"
  expect_error(
    gly_sam(legacy, group_col = "label", s0 = 1),
    "match the assay columns"
  )
})
