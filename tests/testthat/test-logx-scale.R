logx_test_experiment <- function(scale = "logx", groups = 2L) {
  withr::with_seed(823, {
    n <- 60L
    values <- matrix(
      stats::rnorm(8L * n, 0, 0.3),
      nrow = 8L,
      dimnames = list(paste0("v", 1:8), paste0("s", seq_len(n)))
    )
    values[8, 1] <- 0
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

logx_test_methods <- function() {
  list(
    ttest = list(),
    wilcox = list(),
    anova = list(),
    ancova = list(covariate_cols = "age"),
    kruskal = list(),
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

run_logx_test <- function(
  name,
  exp,
  option = NULL,
  args = logx_test_methods()[[name]]
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

test_that("logx resolves to unchanged native coordinates only", {
  exp <- logx_test_experiment()
  for (option in list(NULL, FALSE)) {
    for (default in c(TRUE, FALSE)) {
      plan <- .resolve_transformation(exp, option, default)
      expect_false(plan$log_transform)
      expect_identical(plan$input_scale, "logx")
      expect_identical(plan$analysis_scale, "logx")
      expect_false(plan$scale_missing)
    }
  }
  expect_warning(plan <- .resolve_transformation(exp, TRUE), "log_transform")
  expect_false(plan$log_transform)
  expect_true(plan$requested_log_transform)
  for (option in list(NULL, FALSE, TRUE)) {
    expect_error(
      suppressWarnings(gly_fold_change(exp, log_transform = option)),
      "logx"
    )
  }
})

for (name in names(logx_test_methods())) {
  local({
    method <- name
    test_that(paste(method, "uses native logx coordinates without mutation"), {
      exp <- logx_test_experiment()
      original <- exp
      auto <- run_logx_test(method, exp)
      explicit <- run_logx_test(method, exp, FALSE)
      expect_equal(auto$tidy_result, explicit$tidy_result)
      expect_identical(auto$preprocessing$analysis_scale, "logx")
      expect_false(auto$preprocessing$log_transform)
      expect_identical(auto$meta_data, S4Vectors::metadata(exp))
      expect_identical(exp, original)
      ignored <- run_logx_test(method, exp, TRUE)
      expect_equal(ignored$tidy_result, explicit$tidy_result)
      expect_false(ignored$preprocessing$log_transform)
      expect_true(ignored$preprocessing$requested_log_transform)
      # Compare with the established native-coordinate path, not just the new
      # automatic and explicit branches with each other.
      reference <- exp
      S4Vectors::metadata(reference)$scale <- "log2"
      expected <- run_logx_test(method, reference, FALSE)$tidy_result
      actual <- auto$tidy_result
      if (method %in% c("ttest", "wilcox", "limma")) {
        expect_type(actual$log2fc, "double")
        expect_true(all(is.na(actual$log2fc)))
        actual$log2fc <- expected$log2fc <- NULL
        if (method == "limma") {
          expect_equal(
            actual$estimate,
            unname(auto$raw_result$coefficients[, 1])
          )
          actual$estimate <- expected$estimate <- NULL
        }
      }
      if (method %in% c("anova", "ancova", "kruskal")) {
        expect_true(all(is.na(actual$post_hoc_test$log2fc)))
        actual$post_hoc_test$log2fc <- expected$post_hoc_test$log2fc <- NULL
      }
      expect_equal(actual, expected)
    })
  })
}

for (name in c(
  "ttest",
  "wilcox",
  "anova",
  "ancova",
  "kruskal",
  "limma",
  "linear_model"
)) {
  local({
    method <- name
    test_that(
      paste(method, "filters logx by significance without fold changes"),
      {
        exp <- logx_test_experiment()
        S4Vectors::metadata(exp)$exp_type <- "glycoproteomics"
        res <- run_logx_test(method, exp)
        expect_no_error(filtered <- filter_sig_vars(exp, res, p_adj_cutoff = 1))
        expect_gt(nrow(filtered), 0)
        expect_error(filter_sig_vars(exp, res, fc_cutoff = 2), "fc_cutoff")
        if (method %in% c("anova", "ancova", "kruskal")) {
          post <- res$tidy_result$post_hoc_test
          comparison <- paste0(post$ref_group[1], "_vs_", post$test_group[1])
          expect_no_error(filter_sig_vars(
            exp,
            res,
            p_adj_cutoff = 1,
            on = "post_hoc_test",
            comparison = comparison
          ))
          expect_error(
            filter_sig_vars(
              exp,
              res,
              fc_cutoff = 2,
              on = "post_hoc_test",
              comparison = comparison
            ),
            "fc_cutoff"
          )
        }
      }
    )
  })
}

test_that("paired logx tests retain native estimates with incomplete pairs", {
  exp <- logx_test_experiment()
  SummarizedExperiment::assay(exp)[1, 1] <- NA_real_
  for (method in c("ttest", "wilcox")) {
    res <- run_logx_test(method, exp, args = list(subject_col = "subject"))
    expect_true(all(is.na(res$tidy_result$log2fc)))
    reference <- exp
    S4Vectors::metadata(reference)$scale <- "log2"
    expected <- run_logx_test(
      method,
      reference,
      FALSE,
      args = list(subject_col = "subject")
    )$tidy_result
    actual <- res$tidy_result
    actual$log2fc <- expected$log2fc <- NULL
    expect_equal(actual, expected)
  }
})

for (name in c("anova", "ancova", "kruskal", "limma")) {
  local({
    method <- name
    test_that(paste(method, "handles multigroup and missing logx features"), {
      exp <- logx_test_experiment(groups = 3L)
      SummarizedExperiment::assay(exp)[1, 1] <- NA_real_
      SummarizedExperiment::assay(exp)[2, ] <- NA_real_
      res <- run_logx_test(method, exp)
      post <- if (method == "limma") {
        res$tidy_result
      } else {
        res$tidy_result$post_hoc_test
      }
      expect_type(post$log2fc, "double")
      expect_true(all(is.na(post$log2fc)))
      SummarizedExperiment::assay(exp)[,] <- NA_real_
      expect_no_error(run_logx_test(method, exp))
    })
  })
}

test_that("native glycomic and glycoproteomic logx classes preserve their inputs", {
  base <- logx_test_experiment()
  rd <- S4Vectors::DataFrame(
    protein = paste0("P", seq_len(nrow(base))),
    protein_site = seq_len(nrow(base)),
    glycan_composition = glyrepr::as_glycan_composition(rep(
      "H5N2",
      nrow(base)
    )),
    row.names = rownames(base)
  )
  for (constructor in list(glyexp::GlycomicSE, glyexp::GlycoproteomicSE)) {
    exp <- constructor(
      SummarizedExperiment::assay(base),
      rowData = rd,
      colData = SummarizedExperiment::colData(base),
      metadata = list(scale = "logx", glycan_type = "N")
    )
    before <- exp
    res <- run_logx_test("ttest", exp)
    expect_true(all(is.na(res$tidy_result$log2fc)))
    filtered <- filter_sig_vars(exp, res, p_adj_cutoff = 1)
    expect_identical(class(filtered), class(exp))
    expect_identical(rownames(filtered), rownames(exp))
    expect_identical(exp, before)
    expect_error(filter_sig_vars(exp, res, fc_cutoff = 2), "fc_cutoff")
  }
})

test_that("legacy S3 logx inputs retain native estimates and metadata", {
  withr::local_options(lifecycle_verbosity = "quiet")
  se <- logx_test_experiment()
  legacy <- glyexp::experiment(
    SummarizedExperiment::assay(se),
    sample_info = as.data.frame(SummarizedExperiment::colData(se)),
    var_info = data.frame(row.names = rownames(se)),
    exp_type = "others"
  )
  legacy$meta_data$scale <- "logx"
  original <- legacy
  expect_s3_class(legacy, "glyexp_experiment")
  for (method in c("pca", "ttest")) {
    res <- run_logx_test(method, legacy)
    expect_equal(res$tidy_result, run_logx_test(method, se)$tidy_result)
    expect_identical(res$preprocessing$analysis_scale, "logx")
    expect_false(res$preprocessing$log_transform)
    expect_identical(res$meta_data, legacy$meta_data)
  }
  expect_identical(legacy, original)
  expect_error(gly_fold_change(legacy), "logx")
})
