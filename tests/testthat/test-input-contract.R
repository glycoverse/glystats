test_that("analysis functions use a unified data-container input", {
  analysis_functions <- grep(
    "^gly_",
    getNamespaceExports("glystats"),
    value = TRUE
  )

  expect_false(any(grepl("_$", analysis_functions)))
  purrr::walk(analysis_functions, function(function_name) {
    analysis_function <- getExportedValue("glystats", function_name)
    expect_identical(names(formals(analysis_function))[1], "exp")
  })
})

test_that("analysis functions reject matrix inputs", {
  expr_mat <- matrix(1:4, nrow = 2)

  expect_error(gly_fold_change(expr_mat), "SummarizedExperiment")
  expect_error(gly_pca(expr_mat), "SummarizedExperiment")
  expect_error(gly_cox(expr_mat), "SummarizedExperiment")
})

test_that("analysis functions support SummarizedExperiment natively", {
  se <- glyexp::as_glycoproteomic_se(test_gp_exp)
  two_group_se <- se[, SummarizedExperiment::colData(se)$group %in% c("C", "H")]

  pca_result <- gly_pca(se)
  fold_change_result <- suppressMessages(gly_fold_change(two_group_se))

  expect_s4_class(se, "GlycoproteomicSE")
  expect_s3_class(pca_result, "glystats_pca_res")
  expect_s3_class(fold_change_result, "glystats_fc_res")
  expect_equal(pca_result$meta_data, S4Vectors::metadata(se))
  expect_true("protein" %in% names(fold_change_result$tidy_result))
})

test_that("plain SummarizedExperiment inputs do not require glyexp subclasses", {
  glyco_se <- glyexp::as_glycomic_se(glyexp::real_experiment2)
  abundance <- SummarizedExperiment::assay(glyco_se)
  complete_rows <- rowSums(!is.finite(abundance)) == 0
  varying_rows <- apply(abundance, 1, stats::var) > 0
  se <- methods::as(
    glyco_se[complete_rows & varying_rows, ],
    "SummarizedExperiment"
  )

  result <- gly_pca(se)

  expect_s4_class(se, "SummarizedExperiment")
  expect_false(methods::is(se, "GlycomicSE"))
  expect_s3_class(result, "glystats_pca_res")
  expect_equal(result$meta_data, S4Vectors::metadata(se))
})

test_that("SummarizedExperiment inputs preserve ordinary sample and variable metadata", {
  se <- methods::as(test_gp_se, "SummarizedExperiment")
  SummarizedExperiment::colData(se)$sample <- colnames(se)
  SummarizedExperiment::rowData(se)$variable <- rownames(se)

  result <- gly_pca(se)

  expect_setequal(result$tidy_result$samples$.sample, colnames(se))
  expect_setequal(result$tidy_result$variables$.variable, rownames(se))
  expect_setequal(
    result$tidy_result$samples$sample,
    SummarizedExperiment::colData(se)$sample
  )
  expect_setequal(
    result$tidy_result$variables$variable,
    SummarizedExperiment::rowData(se)$variable
  )
})

test_that("merged glyexp development constructors supply synchronized analysis identities", {
  values <- matrix(
    seq_len(24),
    nrow = 3,
    dimnames = list(c("v1", "v2", "v3"), paste0("s", seq_len(8)))
  )
  rd <- S4Vectors::DataFrame(
    glycan_composition = rep(glyrepr::glycan_composition(c(Hex = 1)), 3),
    variable = c("annotation 3", "annotation 1", "annotation 2"),
    row.names = c("v3", "v1", "v2")
  )
  cd <- S4Vectors::DataFrame(
    group = rep(c("A", "B"), each = 4),
    sample = rep("ordinary sample annotation", 8),
    row.names = colnames(values)
  )
  for (proteomic in c(FALSE, TRUE)) {
    info <- rd
    ctor <- glyexp::GlycomicSE
    if (proteomic) {
      info$protein <- c("p3", "p1", "p2")
      info$protein_site <- c(3L, 1L, 2L)
      ctor <- glyexp::GlycoproteomicSE
    }
    se <- ctor(
      values,
      rowData = info,
      colData = cd,
      metadata = list(glycan_type = "N")
    )
    result <- suppressMessages(gly_ttest(se, log_transform = FALSE))
    expect_identical(result$tidy_result$.variable, rownames(values))
    expect_identical(result$tidy_result$variable, paste("annotation", 1:3))
    expect_identical(
      colnames(SummarizedExperiment::assay(se, withDimnames = FALSE)),
      colnames(se)
    )
    expect_false(result$identifiers$variable$generated)
    expect_false(result$identifiers$sample$generated)
  }
})
