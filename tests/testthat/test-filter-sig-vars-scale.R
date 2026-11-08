raw_limma_filter_experiment <- function(n_groups = 2L) {
  groups <- factor(rep(LETTERS[seq_len(n_groups)], each = 6))
  values <- outer(c(100, 100, 200, 300), rep(1, length(groups))) +
    outer(c(1, 0.8, 1.5, 2.5), rep(c(-2, -1, 0, 0, 1, 2), n_groups))
  values[1, groups == "B"] <- values[1, groups == "B"] + 10
  dimnames(values) <- list(
    c("signal", "control1", "control2", "control3"),
    paste0("s", seq_along(groups))
  )
  glyexp::GlycoproteomicSE(
    values,
    rowData = S4Vectors::DataFrame(
      protein = paste0("P", seq_len(nrow(values))),
      protein_site = seq_len(nrow(values)),
      glycan_composition = glyrepr::as_glycan_composition(rep(
        "H5N2",
        nrow(values)
      )),
      row.names = rownames(values)
    ),
    colData = S4Vectors::DataFrame(
      group = groups,
      row.names = colnames(values)
    ),
    metadata = list(
      glycan_type = "N",
      exp_type = "glycoproteomics",
      scale = "raw"
    )
  )
}

test_that("raw-scale limma filtering keeps significant glycoproteomics variables", {
  skip_if_not_installed("limma")
  glyco <- raw_limma_filter_experiment()
  for (exp in list(glyco, methods::as(glyco, "SummarizedExperiment"))) {
    result <- suppressMessages(gly_limma(
      exp,
      log_transform = FALSE,
      add_info = FALSE
    ))
    original <- exp
    expect_identical(result$preprocessing$analysis_scale, "raw")
    expect_true(all(is.na(result$tidy_result$log2fc)))
    expect_identical(
      result$tidy_result$.variable[result$tidy_result$p_adj < 0.05],
      "signal"
    )
    expect_identical(rownames(filter_sig_vars(exp, result)), "signal")
    expect_identical(
      rownames(filter_sig_vars(exp, result, fc_cutoff = NULL)),
      "signal"
    )
    expect_identical(
      rownames(filter_sig_vars(
        exp,
        result,
        p_adj_cutoff = NULL,
        p_val_cutoff = 0.05
      )),
      "signal"
    )
    expect_error(
      filter_sig_vars(exp, result, fc_cutoff = 2),
      "raw-scale",
      class = "rlang_error"
    )
    expect_error(
      filter_sig_vars(exp, result, p_adj_cutoff = NULL, fc_cutoff = 2),
      "raw-scale",
      class = "rlang_error"
    )
    expect_identical(exp, original)
  }
})

test_that("raw-scale limma filtering respects comparison and p-value selection", {
  skip_if_not_installed("limma")
  exp <- raw_limma_filter_experiment(3L)
  result <- suppressMessages(gly_limma(
    exp,
    log_transform = FALSE,
    add_info = FALSE
  ))
  expect_identical(rownames(filter_sig_vars(exp, result)), "signal")
  expect_identical(
    rownames(filter_sig_vars(exp, result, comparison = "A_vs_B")),
    "signal"
  )
  expect_identical(
    rownames(filter_sig_vars(exp, result, comparison = "B_vs_C")),
    "signal"
  )
  expect_length(
    rownames(filter_sig_vars(exp, result, comparison = "A_vs_C")),
    0
  )
  unadjusted <- suppressMessages(gly_limma(
    exp,
    log_transform = FALSE,
    p_adj_method = NULL,
    add_info = FALSE
  ))
  expect_identical(
    rownames(filter_sig_vars(
      exp,
      unadjusted,
      p_adj_cutoff = NULL,
      p_val_cutoff = 0.05
    )),
    "signal"
  )
})

test_that("logged and legacy limma results retain fold-change filtering", {
  skip_if_not_installed("limma")
  exp <- raw_limma_filter_experiment()
  result <- suppressMessages(gly_limma(exp, add_info = FALSE))
  expect_identical(result$preprocessing$analysis_scale, "log2")
  expect_length(rownames(filter_sig_vars(exp, result)), 0)
  expect_identical(
    rownames(filter_sig_vars(exp, result, fc_cutoff = 1)),
    "signal"
  )
  result$preprocessing <- NULL
  expect_length(rownames(filter_sig_vars(exp, result)), 0)
  expect_identical(
    rownames(filter_sig_vars(exp, result, fc_cutoff = 1)),
    "signal"
  )
})

for (method in c("ttest", "wilcox", "limma", "anova", "ancova", "kruskal")) {
  local({
    name <- method
    test_that(
      paste(name, "distinguishes omitted FC defaults from explicit NULL"),
      {
        # Fixed results isolate filtering from fit variability: all variables are
        # significant but have a log2 fold change below the default cutoff.
        exp <- raw_limma_filter_experiment()
        tab <- tibble::tibble(
          variable = rownames(exp),
          p_val = 0.001,
          p_adj = 0.001,
          log2fc = 0.5,
          ref_group = "A",
          test_group = "B"
        )
        multi <- name %in% c("anova", "ancova", "kruskal")
        tidy <- if (multi) list(main_test = tab, post_hoc_test = tab) else tab
        if (name == "limma") {
          tidy$comparison <- "A_vs_B"
        }
        res <- structure(
          list(tidy_result = tidy),
          class = paste0("glystats_", name, "_res")
        )
        modes <- if (multi) {
          list(list(on = "post_hoc_test", comparison = "A_vs_B"))
        } else {
          list(list())
        }
        for (scale in c("raw", "log", "log2", "log10")) {
          S4Vectors::metadata(exp)$scale <- scale
          res$preprocessing <- list(
            analysis_scale = if (scale == "raw") "log2" else scale
          )
          for (mode in modes) {
            args <- c(list(exp = exp, res = res), mode)
            original <- exp
            expect_identical(nrow(do.call(filter_sig_vars, args)), 0L)
            explicit <- do.call(
              filter_sig_vars,
              c(args, list(fc_cutoff = NULL))
            )
            expect_identical(rownames(explicit), rownames(exp))
            expect_identical(
              S4Vectors::metadata(explicit),
              S4Vectors::metadata(exp)
            )
            expect_identical(exp, original)
            expect_error(
              do.call(
                filter_sig_vars,
                c(args, list(fc_cutoff = NULL, p_adj_cutoff = NULL))
              ),
              "All are NULL"
            )
            generic <- methods::as(exp, "SummarizedExperiment")
            S4Vectors::metadata(generic)$exp_type <- "others"
            args$exp <- generic
            expect_identical(
              do.call(filter_sig_vars, args),
              do.call(filter_sig_vars, c(args, list(fc_cutoff = NULL)))
            )
            expect_identical(nrow(do.call(filter_sig_vars, args)), nrow(exp))
          }
          if (multi) {
            expect_identical(
              rownames(filter_sig_vars(exp, res, fc_cutoff = NULL)),
              rownames(exp)
            )
          }
        }
      }
    )
  })
}

test_that("known-scale linear models still reject explicit fold-change cutoffs", {
  exp <- raw_limma_filter_experiment()
  res <- structure(
    list(
      tidy_result = tibble::tibble(
        variable = rownames(exp),
        term = "groupB",
        p_adj = 0.001,
        p_val = 0.001
      ),
      preprocessing = list(analysis_scale = "log2")
    ),
    class = "glystats_linear_model_res"
  )
  for (scale in c("raw", "log", "log2", "log10")) {
    S4Vectors::metadata(exp)$scale <- scale
    res$preprocessing$analysis_scale <- scale
    expect_error(filter_sig_vars(exp, res, fc_cutoff = 2), "fc_cutoff")
    expect_identical(
      filter_sig_vars(exp, res),
      filter_sig_vars(exp, res, fc_cutoff = NULL)
    )
    expect_identical(
      rownames(filter_sig_vars(exp, res, fc_cutoff = NULL)),
      rownames(exp)
    )
  }
})
