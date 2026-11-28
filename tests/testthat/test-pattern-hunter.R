pattern_experiment <- function() {
  group <- factor(
    c("M", "A", "Z", "A", "A", "Z", "M", "A", "Z"),
    levels = c("unused", "Z", "A", "M")
  )
  score <- c(1, 0, 3, 0, 0, 3, 1, 0, 3)
  values <- rbind(
    up = 2 + score,
    down = 8 - score,
    noisy = c(5, 3, 4, 2, 6, 8, 7, 3, 6)
  )
  colnames(values) <- paste0("s", seq_along(group))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(group = group),
    metadata = list(scale = "log2")
  )
}

test_that("Pearson matches use individual samples in unequal-sized groups", {
  exp <- pattern_experiment()
  template <- c(Z = 3, M = 1, A = 0)
  result <- gly_pattern_hunter(exp, template)
  tidy <- get_tidy_result(result)
  target <- c(1, 0, 3, 0, 0, 3, 1, 0, 3)
  expected <- stats::cor.test(
    SummarizedExperiment::assay(exp)["noisy", ],
    target
  )

  expect_s3_class(result, "glystats_pattern_hunter_res")
  expect_s3_class(result, "glystats_res")
  expect_identical(tidy$.variable, rownames(exp))
  expect_equal(tidy$cor, c(1, -1, unname(expected$estimate)))
  expect_equal(tidy$p_val[3], expected$p.value)
  expect_identical(tidy$n_samples, rep(9L, 3))
  expect_identical(tidy$n_groups, rep(3L, 3))
  expect_identical(tidy$n_scores, rep(3L, 3))
  expect_identical(result$template, template)
  expect_named(get_raw_result(result), rownames(exp))
  expect_s3_class(get_raw_result(result, "noisy"), "htest")
  expect_equal(unname(get_raw_result(result, "noisy")$parameter), 7)
})

test_that("reversing a signed template preserves two-sided inference", {
  exp <- pattern_experiment()
  template <- c(A = -1, M = 0, Z = 2)
  for (method in c("pearson", "spearman", "kendall")) {
    forward <- get_tidy_result(gly_pattern_hunter(
      exp,
      template,
      method = method
    ))
    reverse <- get_tidy_result(gly_pattern_hunter(
      exp,
      -template,
      method = method
    ))
    expect_equal(reverse$cor, -forward$cor)
    expect_equal(reverse$p_val, forward$p_val)
    expect_equal(reverse$p_adj, forward$p_adj)
  }
})

test_that("sample order, template order and factor levels do not change the match", {
  exp <- pattern_experiment()
  template <- c(A = 0, M = 1, Z = 3)
  expected <- get_tidy_result(gly_pattern_hunter(exp, template))
  exp <- exp[, c(9, 1, 4, 3, 6, 2, 8, 5, 7)]
  # This ordinary sample annotation supplies the grouping variable.
  groups <- as.character(SummarizedExperiment::colData(exp)$group)
  labels <- c(A = "before dose", M = "dose-1", Z = "dose-2")
  SummarizedExperiment::colData(exp)$sample <- factor(
    unname(labels[groups]),
    levels = c("dose-1", "before dose", "dose-2", "absent")
  )
  result <- gly_pattern_hunter(
    exp,
    c(`dose-2` = 3, `before dose` = 0, `dose-1` = 1),
    group_col = "sample"
  )
  expect_equal(get_tidy_result(result), expected)
})

test_that("finite observations determine template coverage and the adjustment family", {
  values <- rbind(
    empty = c(NA, Inf, -Inf, NaN, NA, NA),
    constant = rep(7, 6),
    one_group = c(1, 2, 3, NA, NA, NA),
    one_score = c(1, 2, 3, 4, NA, NA),
    two_samples = c(1, NA, NA, NA, 5, NA),
    partial = c(1, NA, 2, Inf, 3, 4),
    full = c(4, 1, 3, 2, 5, 3)
  )
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(group = c("A", "A", "A", "B", "C", "C"))
  )
  template <- c(C = 1, B = 0, A = 0)
  result <- gly_pattern_hunter(exp, template, log_transform = FALSE)
  tidy <- get_tidy_result(result)
  expect_identical(tidy$.variable, rownames(values))
  expect_equal(tidy$n_samples, c(0, 6, 3, 4, 2, 4, 6))
  expect_equal(tidy$n_groups, c(0, 3, 1, 2, 2, 2, 3))
  expect_equal(tidy$n_scores, c(0, 2, 1, 1, 2, 2, 2))
  expect_true(all(is.na(tidy[1:5, c("cor", "p_val", "p_adj")])))
  expect_true(all(vapply(get_raw_result(result)[1:5], is.null, logical(1))))

  partial <- stats::cor.test(1:4, c(0, 0, 1, 1))
  full <- stats::cor.test(values["full", ], c(0, 0, 0, 0, 1, 1))
  p <- c(partial$p.value, full$p.value)
  expect_equal(tidy$cor[6:7], unname(c(partial$estimate, full$estimate)))
  expect_equal(tidy$p_val[6:7], p)
  expect_equal(tidy$p_adj[6:7], stats::p.adjust(p, "BH"))
  holm <- gly_pattern_hunter(
    exp,
    template,
    log_transform = FALSE,
    p_adj_method = "holm"
  )
  expect_equal(get_tidy_result(holm)$p_adj[6:7], stats::p.adjust(p, "holm"))
  no_adj <- gly_pattern_hunter(
    exp,
    template,
    log_transform = FALSE,
    p_adj_method = NULL
  )
  expect_false("p_adj" %in% names(get_tidy_result(no_adj)))
  empty <- gly_pattern_hunter(exp[integer(), ], template)
  expect_identical(get_tidy_result(empty), tidy[integer(), ])
})

test_that("rank tests handle ties in scores and abundances and untied small samples", {
  exp <- pattern_experiment()
  template <- c(M = 0, A = 0, Z = 1)
  target <- c(0, 0, 1, 0, 0, 1, 0, 0, 1)
  x <- SummarizedExperiment::assay(exp)["noisy", ]
  small <- exp["up", c(2, 1, 3)]
  for (method in c("spearman", "kendall")) {
    expect_no_warning(
      result <- gly_pattern_hunter(exp, template, method = method)
    )
    tidy <- get_tidy_result(result)
    expected <- stats::cor.test(x, target, method = method, exact = FALSE)
    expect_equal(tidy$cor[3], unname(expected$estimate))
    expect_equal(tidy$p_val[3], expected$p.value)
    expect_identical(tidy$n_scores, rep(2L, 3))
    expect_identical(tidy$n_groups, rep(3L, 3))

    untied <- gly_pattern_hunter(small, c(A = 0, M = 1, Z = 3), method = method)
    expect_equal(get_tidy_result(untied)$cor, 1)
    expect_equal(get_tidy_result(untied)$p_val, 2 / factorial(3))
  }
})

test_that("templates and sample group assignments are validated before testing", {
  exp <- pattern_experiment()
  valid <- c(A = 0, M = 1, Z = 3)
  for (bad in list(
    unname(valid),
    valid[-1],
    c(valid, unused = 2),
    c(a = 0, M = 1, Z = 3),
    c(valid, A = 2),
    setNames(valid, c("A", "", "Z")),
    setNames(valid, c("A", NA, "Z")),
    c(A = 2, M = 2, Z = 2),
    c(A = 0, M = Inf, Z = 1),
    c(A = -Inf, M = 0, Z = 1),
    c(A = NA_real_, M = 0, Z = 1),
    c(A = 0, M = NaN, Z = 1),
    c(A = "0", M = "1", Z = "3"),
    setNames(matrix(valid), names(valid))
  )) {
    expect_error(gly_pattern_hunter(exp, bad), "template")
  }
  expect_error(gly_pattern_hunter(exp, valid, group_col = "absent"), "absent")
  expect_error(gly_pattern_hunter(exp, valid, method = "cosine"), "method")
  expect_error(
    gly_pattern_hunter(exp, valid, p_adj_method = "bad"),
    "p_adj_method"
  )
  expect_error(
    gly_pattern_hunter(SummarizedExperiment::assay(exp), valid),
    "SummarizedExperiment"
  )

  groups <- as.character(SummarizedExperiment::colData(exp)$group)
  for (bad in c(NA_character_, "")) {
    SummarizedExperiment::colData(exp)$group <- replace(groups, 1, bad)
    expect_error(gly_pattern_hunter(exp, valid), "group labels")
  }
})

test_that("pattern matching follows the recorded analysis scale", {
  exp <- pattern_experiment()
  template <- c(A = 0, M = 1, Z = 3)
  target <- c(1, 0, 3, 0, 0, 3, 1, 0, 3)
  values <- SummarizedExperiment::assay(exp)["noisy", ]
  S4Vectors::metadata(exp)$scale <- "raw"
  logged <- gly_pattern_hunter(exp, template)
  expect_identical(logged$preprocessing$analysis_scale, "log2")
  expect_equal(
    get_tidy_result(logged)$cor[3],
    stats::cor(log2(values + 1e-6), target)
  )
  raw <- gly_pattern_hunter(exp, template, log_transform = FALSE)
  expect_equal(get_tidy_result(raw)$cor[3], stats::cor(values, target))
  S4Vectors::metadata(exp)$scale <- "logx"
  SummarizedExperiment::assay(exp) <- SummarizedExperiment::assay(exp) - 5
  expect_warning(
    native <- gly_pattern_hunter(exp, template, log_transform = TRUE),
    class = "glystats_already_logged"
  )
  expect_identical(native$meta_data, S4Vectors::metadata(exp))
  expect_identical(native$preprocessing$analysis_scale, "logx")
  expect_equal(get_tidy_result(native)$cor, get_tidy_result(raw)$cor)
})

test_that("SE position IDs and annotation collisions survive result extraction", {
  exp <- pattern_experiment()
  template <- c(Z = 3, M = 1, A = 0)
  rownames(exp) <- c("duplicate", "duplicate", "other")
  colnames(exp) <- NULL
  SummarizedExperiment::rowData(exp)$variable <- letters[1:3]
  SummarizedExperiment::rowData(exp)$cor <- 11:13
  SummarizedExperiment::rowData(exp)$p_adj <- 21:23
  result <- gly_pattern_hunter(exp, template, p_adj_method = NULL)
  tidy <- get_tidy_result(result)
  ids <- paste0(".glystats_variable_", 1:3)
  expect_identical(tidy$.variable, ids)
  expect_identical(tidy$variable, letters[1:3])
  expect_identical(tidy$cor_metadata, 11:13)
  expect_identical(tidy$p_adj_metadata, 21:23)
  expect_false("p_adj" %in% names(tidy))
  expect_named(get_raw_result(result), ids)
  expect_equal(tidy$cor[1:2], c(1, -1))
  expect_identical(result$identifiers$variable$data$name, rownames(exp))
  expect_true(result$identifiers$sample$generated)
  bare <- gly_pattern_hunter(exp, template, add_info = FALSE)
  expect_identical(get_tidy_result(bare)$.variable, ids)
  expect_false("variable" %in% names(get_tidy_result(bare)))
})

test_that("real glycomic and glycoproteomic SEs use the public result contract", {
  for (exp in list(glyexp::real_experiment2, glyexp::real_experiment)) {
    group <- SummarizedExperiment::colData(exp)$group
    exp <- exp[1:5, group %in% c("H", "C")]
    target <- as.integer(SummarizedExperiment::colData(exp)$group == "C")
    result <- glystats::gly_pattern_hunter(exp, c(C = 1, H = 0))
    tidy <- glystats::get_tidy_result(result)
    values <- SummarizedExperiment::assay(exp)
    expect_identical(tidy$.variable, rownames(exp))
    expect_identical(result$identifiers$sample$data$id, colnames(exp))
    expect_true("glycan_composition" %in% names(tidy))
    expect_equal(tidy$n_samples, unname(rowSums(is.finite(values))))
    i <- which.max(tidy$n_samples)
    expected <- stats::cor.test(log2(values[i, ] + 1e-6), target)
    expect_equal(tidy$cor[i], unname(expected$estimate))
    expect_equal(tidy$p_val[i], expected$p.value)
  }
})

test_that("legacy experiments retain matching statistics and undotted IDs", {
  withr::local_options(lifecycle_verbosity = "quiet")
  legacy <- glyexp::slice_head_var(exp_2groups(), n = 3)
  template <- c(H = 0, C = 1)
  original <- gly_pattern_hunter(legacy, template, add_info = FALSE)
  converted <- gly_pattern_hunter(
    as_test_se(legacy),
    template,
    add_info = FALSE
  )
  tidy <- get_tidy_result(original)
  se_tidy <- get_tidy_result(converted)
  expect_identical(tidy$variable, se_tidy$.variable)
  expect_equal(tidy[-1], se_tidy[-1])
  expect_true(original$preprocessing$scale_missing)
})
