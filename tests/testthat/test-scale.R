scale_test_experiment <- function(scale = "raw", groups = 2L) {
  withr::with_seed(823, {
    n <- 60L
    values <- matrix(
      stats::rnorm(8L * n, 3, 0.3),
      nrow = 8L,
      dimnames = list(paste0("v", 1:8), paste0("s", seq_len(n)))
    )
    group <- factor(rep(LETTERS[seq_len(groups)], each = n / groups))
    values[1:4, group == "B"] <- values[1:4, group == "B"] + 0.5
    abundance <- switch(
      scale,
      raw = 2^values,
      log = values * log(2),
      log2 = values,
      log10 = values * log10(2),
      values
    )
    SummarizedExperiment::SummarizedExperiment(
      assays = list(abundance = abundance),
      colData = S4Vectors::DataFrame(
        group = group,
        age = sample(seq_len(n)),
        time = stats::rexp(n),
        event = rep(0:1, length.out = n),
        subject = rep(seq_len(n / groups), groups),
        row.names = colnames(values)
      ),
      metadata = list(scale = scale, source = "scale test")
    )
  })
}

scale_test_methods <- function() {
  list(
    ttest = list(),
    wilcox = list(),
    anova = list(),
    ancova = list(covariate_cols = "age"),
    kruskal = list(),
    fold_change = list(),
    limma = list(),
    linear_model = list(formula = ~group),
    set_test = list(sets = list(pair = c("v1", "v2"))),
    pca = list(),
    tsne = list(perplexity = 3, max_iter = 250),
    umap = list(n_neighbors = 5, n_epochs = 10, n_threads = 1),
    plsda = list(
      ncomp = 1,
      crossvalI = 2,
      permI = 0,
      fig.pdfC = "none",
      info.txtC = "none"
    ),
    oplsda = list(
      pred_i = 1,
      ortho_i = 0,
      crossvalI = 2,
      permI = 0,
      fig.pdfC = "none",
      info.txtC = "none"
    ),
    cor = list(),
    hclust = list(k_values = 2),
    kmeans = list(centers = 2),
    cox = list(),
    roc = list()
  )
}

run_scale_test <- function(
  name,
  exp,
  option = NULL,
  args = scale_test_methods()[[name]]
) {
  dependency <- switch(
    name,
    limma = "limma",
    linear_model = "limma",
    ancova = "emmeans",
    kruskal = "FSA",
    cor = "Hmisc",
    plsda = "ropls",
    oplsda = "ropls",
    tsne = "Rtsne",
    umap = "uwot",
    roc = "pROC",
    NULL
  )
  if (!is.null(dependency)) {
    skip_if_not_installed(dependency)
  }
  withr::with_seed(
    523,
    suppressMessages(suppressWarnings(do.call(
      get(paste0("gly_", name), envir = asNamespace("glystats")),
      c(list(exp = exp, log_transform = option), args)
    )))
  )
}

test_that("scale resolution covers defaults, overrides, and invalid metadata", {
  for (scale in c("raw", "log", "log2", "log10")) {
    exp <- scale_test_experiment(scale)
    for (option in list(NULL, TRUE, FALSE)) {
      for (default in c(TRUE, FALSE)) {
        plan <- .resolve_transformation(exp, option, default)
        expected <- if (is.null(option)) scale == "raw" && default else option
        expect_identical(plan$log_transform, expected)
        expect_identical(plan$input_scale, scale)
        expect_identical(plan$analysis_scale, if (expected) "log2" else scale)
        expect_false(plan$scale_missing)
      }
    }
  }
  exp <- scale_test_experiment()
  S4Vectors::metadata(exp)$scale <- NULL
  expect_true(.resolve_transformation(exp, NULL)$scale_missing)
  expect_true(.resolve_transformation(exp, NULL)$log_transform)
  expect_false(.resolve_transformation(exp, NULL, FALSE)$log_transform)

  for (bad in list(
    NULL,
    NA_character_,
    "log3",
    "LOG2",
    c("raw", "log"),
    2,
    character()
  )) {
    S4Vectors::metadata(exp) <- list(scale = bad)
    for (option in list(NULL, TRUE, FALSE)) {
      expect_error(gly_fold_change(exp, log_transform = option), "scale")
    }
  }
  exp <- scale_test_experiment()
  for (bad in list(NA, logical(), c(TRUE, FALSE), 1, "auto")) {
    expect_error(gly_pca(exp, log_transform = bad), "log_transform")
  }
})

test_that("every public entry routes scales and options without forwarding them through dots", {
  # Intercept the fit boundary to exercise every wrapper across the full contract
  # without repeating expensive stochastic fits hundreds of times.
  for (name in names(scale_test_methods())) {
    capture <- function(expr_mat, ..., log_transform) {
      expect_identical(expr_mat, .get_expr_mat(exp))
      expect_identical(log_transform, expected)
      if (name == "fold_change") {
        tibble::tibble(variable = rownames(expr_mat), log2fc = 0)
      } else if (name == "linear_model") {
        list(
          tidy_result = tibble::tibble(variable = rownames(expr_mat)),
          raw_result = list()
        )
      } else {
        list(tidy_result = tibble::tibble(), raw_result = list())
      }
    }
    bindings <- stats::setNames(list(capture), paste0(".analyze_", name))
    for (scale in c("raw", "log", "log2", "log10")) {
      exp <- scale_test_experiment(scale)
      before <- exp
      for (option in list(NULL, TRUE, FALSE)) {
        expected <- if (is.null(option)) {
          scale == "raw" && !name %in% c("cox", "roc", "fold_change")
        } else {
          option
        }
        result <- testthat::with_mocked_bindings(
          run_scale_test(name, exp, option),
          !!!bindings
        )
        expect_identical(result$meta_data, S4Vectors::metadata(exp))
        expect_identical(result$preprocessing$log_transform, expected)
        expect_identical(exp, before)
      }
    }
  }
})

test_that("the log2 helper retains pseudo-counts and edge values", {
  values <- matrix(c(0, -2, NA_real_, 3), nrow = 1)
  expect_equal(
    suppressWarnings(.log_transform_expr_mat(values)),
    suppressWarnings(log2(values + 1e-6))
  )
  expect_identical(.log_transform_expr_mat(values, FALSE), values)
})

test_that("logged fold changes convert mean differences and retain group direction", {
  groups <- factor(rep(c("A", "B", "C"), each = 3))
  log_values <- matrix(
    c(-3, -1, 0, 0, 2, 3, 2, 4, 5),
    nrow = 1,
    dimnames = list("v", paste0("s", 1:9))
  )
  expected <- c(3, 5, 2)
  for (scale in c("log", "log2", "log10")) {
    values <- log_values / .log2_factor(scale)
    exp <- SummarizedExperiment::SummarizedExperiment(
      assays = list(abundance = values),
      colData = S4Vectors::DataFrame(group = groups),
      metadata = list(scale = scale)
    )
    expect_equal(
      unname(run_scale_test("fold_change", exp)$tidy_result$log2fc),
      expected
    )
    expect_equal(
      unname(run_scale_test("fold_change", exp, FALSE)$tidy_result$log2fc),
      expected
    )
    two <- exp[, groups != "C"]
    for (name in c("ttest", "wilcox")) {
      result <- run_scale_test(name, two)
      expect_equal(result$tidy_result$log2fc, 3)
      result <- run_scale_test(name, two, args = list(ref_group = "B"))
      expect_equal(result$tidy_result$log2fc, -3)
    }
  }
  exp <- scale_test_experiment()
  values <- .get_expr_mat(exp)
  groups <- .extract_groups(.get_sample_info(exp), "group")
  expected <- log2(
    rowMeans(values[, groups == "B"]) / rowMeans(values[, groups == "A"])
  )
  expect_equal(
    unname(run_scale_test("fold_change", exp)$tidy_result$log2fc),
    unname(expected)
  )
  logged <- log2(values + 1e-6)
  expected <- rowMeans(logged[, groups == "B"]) -
    rowMeans(logged[, groups == "A"])
  expect_equal(
    unname(run_scale_test("fold_change", exp, TRUE)$tidy_result$log2fc),
    unname(expected)
  )
})

test_that("two-group tests use supplied log units consistently for estimates and intervals", {
  results <- lapply(c("log", "log2", "log10"), function(scale) {
    exp <- scale_test_experiment(scale)
    result <- run_scale_test("ttest", exp)
    values <- .get_expr_mat(exp)[1, ]
    expected <- stats::t.test(values[31:60], values[1:30])
    row <- result$tidy_result[1, ]
    expect_equal(row$estimate, unname(diff(rev(expected$estimate))))
    expect_equal(row$conf_low, expected$conf.int[1])
    expect_equal(row$conf_high, expected$conf.int[2])
    expect_equal(row$p_val, expected$p.value)
    result
  })
  expect_equal(results[[1]]$tidy_result$log2fc, results[[2]]$tidy_result$log2fc)
  expect_equal(results[[3]]$tidy_result$log2fc, results[[2]]$tidy_result$log2fc)
})

test_that("paired and multi-group post-hoc fold changes honor all log bases", {
  for (scale in c("log", "log2", "log10")) {
    exp <- scale_test_experiment(scale, groups = 3L)
    # Missing pairs must be selected on the fitted variables, not other samples.
    SummarizedExperiment::assay(exp)[1, 1] <- NA_real_
    values <- .get_expr_mat(exp)
    values_log2 <- values * .log2_factor(scale)
    for (name in c("anova", "ancova", "kruskal")) {
      args <- scale_test_methods()[[name]]
      for (paired in if (name == "ancova") FALSE else c(FALSE, TRUE)) {
        if (paired) {
          args$subject_col <- "subject"
        }
        result <- run_scale_test(name, exp, args = args)
        post <- result$tidy_result$post_hoc_test
        expect_gt(nrow(post), 0)
        sample_groups <- .get_sample_info(exp)$group
        expected <- vapply(
          seq_len(nrow(post)),
          function(i) {
            values_i <- values_log2[post$variable[i], ]
            if (paired) {
              valid_subject <- apply(
                matrix(is.finite(values_i), ncol = 3),
                1,
                all
              )
              values_i[!rep(valid_subject, 3)] <- NA_real_
            }
            mean(
              values_i[sample_groups == post$test_group[i]],
              na.rm = paired
            ) -
              mean(values_i[sample_groups == post$ref_group[i]], na.rm = paired)
          },
          numeric(1)
        )
        expect_equal(unname(post$log2fc), unname(expected))
      }
    }
  }
})

test_that("limma contrasts convert log2fc while fits retain native units", {
  skip_if_not_installed("limma")
  results <- lapply(c("log", "log2", "log10"), function(scale) {
    exp <- scale_test_experiment(scale)
    result <- run_scale_test("limma", exp)
    coefficient <- unname(result$raw_result$coefficients[, 1])
    expect_equal(result$tidy_result$log2fc, coefficient * .log2_factor(scale))
    if (scale != "log2") {
      expect_equal(result$tidy_result$estimate, coefficient)
    }
    expect_equal(
      run_scale_test("limma", exp, FALSE)$tidy_result,
      result$tidy_result
    )
    result
  })
  expect_equal(results[[1]]$tidy_result$log2fc, results[[2]]$tidy_result$log2fc)
  expect_equal(results[[3]]$tidy_result$log2fc, results[[2]]$tidy_result$log2fc)
  exp <- scale_test_experiment()
  result <- run_scale_test("limma", exp, FALSE)
  expect_equal(
    result$tidy_result$estimate,
    unname(result$raw_result$coefficients[, 1])
  )
  expect_true(all(is.na(result$tidy_result$log2fc)))
  SummarizedExperiment::assay(exp)[,] <- NA_real_
  for (scale in c("raw", "log", "log2", "log10")) {
    S4Vectors::metadata(exp)$scale <- scale
    expect_no_error(run_scale_test("limma", exp, FALSE))
  }
})

for (name in names(scale_test_methods())) {
  local({
    method <- name
    test_that(
      paste(method, "fits logged data once and honors explicit log controls"),
      {
        exp <- scale_test_experiment("log2")
        original <- exp
        auto <- run_scale_test(method, exp)
        unchanged <- run_scale_test(method, exp, FALSE)
        expect_equal(auto$tidy_result, unchanged$tidy_result)
        expect_identical(auto$preprocessing$analysis_scale, "log2")
        expect_false(auto$preprocessing$log_transform)
        forced <- run_scale_test(method, exp, TRUE)
        expect_true(forced$preprocessing$log_transform)
        manual <- exp
        values <- .get_expr_mat(exp)
        natural <- method %in% c("pca", "tsne", "umap", "plsda", "oplsda")
        SummarizedExperiment::assay(manual) <- if (natural) {
          log(values + 1)
        } else {
          log2(values + 1e-6)
        }
        S4Vectors::metadata(manual)$scale <- if (natural) "log" else "log2"
        manual_result <- run_scale_test(method, manual, FALSE)
        # DEA fold changes use original inputs, independently of test transforms.
        if (method %in% c("ttest", "wilcox")) {
          forced$tidy_result$log2fc <- manual_result$tidy_result$log2fc <- NULL
        }
        if (method %in% c("anova", "ancova", "kruskal")) {
          forced$tidy_result$post_hoc_test$log2fc <- NULL
          manual_result$tidy_result$post_hoc_test$log2fc <- NULL
        }
        expect_equal(forced$tidy_result, manual_result$tidy_result)
        expect_identical(exp, original)
      }
    )
  })
}

test_that("legacy S3 metadata and missing scale retain historical behavior", {
  legacy <- exp_2groups()
  before <- legacy
  expect_s3_class(legacy, "glyexp_experiment")
  result <- run_scale_test("ttest", legacy)
  expect_true(result$preprocessing$scale_missing)
  expect_true(result$preprocessing$log_transform)
  expect_identical(legacy, before)
  for (method in c("cox", "roc", "fold_change")) {
    exp <- scale_test_experiment()
    S4Vectors::metadata(exp)$scale <- NULL
    auto <- run_scale_test(method, exp)
    unchanged <- run_scale_test(method, exp, FALSE)
    expect_false(auto$preprocessing$log_transform)
    expect_equal(auto$tidy_result, unchanged$tidy_result)
  }
})


test_that("fold changes retain raw zero and missing-value behavior", {
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      abundance = matrix(
        c(0, 0, NA, 1, 0, 1),
        nrow = 3,
        dimnames = list(c("zero_ref", "both_zero", "missing"), c("a", "b"))
      )
    ),
    colData = S4Vectors::DataFrame(group = c("A", "B")),
    metadata = list(scale = "raw")
  )
  result <- run_scale_test("fold_change", exp)
  expect_equal(unname(result$tidy_result$log2fc), c(Inf, NaN, NA_real_))
  expect_true(all(is.finite(run_scale_test(
    "fold_change",
    exp,
    TRUE
  )$tidy_result$log2fc[1:2])))
})

test_that("legacy S3 log bases use the same scale path as SE inputs", {
  legacy <- exp_2groups()
  for (scale in c("log", "log2", "log10")) {
    legacy$meta_data$scale <- scale
    legacy$expr_mat <- matrix(
      seq_len(length(legacy$expr_mat)),
      nrow = nrow(legacy$expr_mat),
      dimnames = dimnames(legacy$expr_mat)
    )
    legacy$expr_mat[1, 1] <- -1
    result <- run_scale_test("ttest", legacy)
    expect_identical(result$preprocessing$input_scale, scale)
    expect_false(result$preprocessing$log_transform)
    se <- suppressWarnings(glyexp::as_se(legacy))
    expect_equal(run_scale_test("ttest", se)$tidy_result, result$tidy_result)
  }
})

test_that("automatic set construction uses the recorded scale once", {
  exp <- scale_test_experiment("log10")
  result <- run_scale_test("set_test", exp, args = list())
  expected <- .construct_correlated_sets(
    .get_expr_mat(exp),
    strata = rep("all", nrow(exp)),
    threshold = 0.9,
    correlation = "pearson",
    clustering = "connected",
    log_transform = FALSE
  )
  expect_equal(
    result$raw_result$set_construction$correlation_matrix,
    expected$correlation_matrix
  )
  expect_false(result$preprocessing$log_transform)
})
