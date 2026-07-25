make_linear_model_exp <- function() {
  set.seed(2026)
  n_per_cell <- 8
  sample_info <- data.frame(
    treatment = factor(rep(c("A", "B"), each = n_per_cell * 2)),
    time = factor(rep(rep(c("T1", "T2"), each = n_per_cell), 2)),
    age = seq(30, 61),
    batch = factor(rep(c("X", "Y"), times = 16))
  )
  rownames(sample_info) <- paste0("S", seq_len(nrow(sample_info)))
  design <- stats::model.matrix(~ treatment * time + age + batch, sample_info)

  coefficients <- rbind(
    c(3, 1, 0.5, 0.01, 0, 1.5),
    c(3, 0, 0, 0, 0, 0),
    c(2, -0.5, 0.2, 0, 0.2, 0),
    c(4, 0.3, -0.2, -0.01, 0, 0.4),
    c(3, 0, 0, 0, 0, 0),
    c(2, 0.1, 0.1, 0, -0.1, 0),
    c(3, -0.2, 0, 0.01, 0, -0.3),
    c(4, 0, 0.4, 0, 0, 0)
  )
  log_expr <- coefficients %*%
    t(design) +
    matrix(
      stats::rnorm(8 * nrow(sample_info), sd = 0.08),
      nrow = 8
    )
  expr_mat <- 2^log_expr
  rownames(expr_mat) <- paste0("V", seq_len(nrow(expr_mat)))
  colnames(expr_mat) <- rownames(sample_info)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr_mat),
    rowData = S4Vectors::DataFrame(
      feature_type = rep(c("signal", "background"), each = 4),
      row.names = rownames(expr_mat)
    ),
    colData = S4Vectors::DataFrame(sample_info),
    metadata = list(exp_type = "others", study = "linear-model-test")
  )
}

test_that("gly_linear_model returns moderated coefficient results", {
  exp <- make_linear_model_exp()
  result <- gly_linear_model(
    exp,
    ~ treatment * time + age + batch,
    add_info = TRUE
  )

  expect_s3_class(
    result,
    c("glystats_linear_model_res", "glystats_res")
  )
  expect_named(result, c("tidy_result", "raw_result", "meta_data"))
  expect_named(
    result$raw_result,
    c("fit", "design", "contrast_matrix", "coefficient_mapping")
  )
  expect_identical(get_tidy_result(result), result$tidy_result)
  expect_identical(
    get_raw_result(result, "design"),
    result$raw_result$design
  )
  expect_null(result$raw_result$contrast_matrix)
  expect_identical(
    unique(result$tidy_result$term),
    c("treatmentB", "timeT2", "age", "batchY", "treatmentB:timeT2")
  )
  expect_identical(
    result$tidy_result$variable,
    rep(rownames(exp), times = 5)
  )
  expect_named(
    result$tidy_result,
    c(
      "variable",
      "feature_type",
      "term",
      "estimate",
      "AveExpr",
      "t",
      "p_val",
      "p_adj",
      "b"
    )
  )
  expect_equal(result$meta_data, S4Vectors::metadata(exp))
})

test_that("gly_linear_model agrees with a direct limma fit", {
  exp <- make_linear_model_exp()
  result <- gly_linear_model(
    exp,
    ~ treatment * time + age + batch,
    add_info = FALSE
  )

  expression <- log2(SummarizedExperiment::assay(exp) + 1e-6)
  design <- stats::model.matrix(
    ~ treatment * time + age + batch,
    data = as.data.frame(SummarizedExperiment::colData(exp))
  )
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  reference_fit <- limma::lmFit(expression, design) |>
    limma::eBayes(trend = TRUE)
  reference <- limma::topTable(
    reference_fit,
    coef = "treatmentB",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )
  observed <- dplyr::filter(result$tidy_result, .data$term == "treatmentB")

  expect_equal(observed$estimate, reference$logFC)
  expect_equal(observed$t, reference$t)
  expect_equal(observed$p_val, reference$P.Value)
  expect_equal(observed$p_adj, reference$adj.P.Val)
})

test_that("gly_linear_model supports named interaction contrasts", {
  exp <- make_linear_model_exp()
  coefficient_result <- gly_linear_model(
    exp,
    ~ treatment * time + age + batch,
    add_info = FALSE
  )
  contrast_result <- gly_linear_model(
    exp,
    ~ treatment * time + age + batch,
    contrasts = c(
      B_at_T2_vs_A_at_T2 = "treatmentB + `treatmentB:timeT2`"
    ),
    add_info = FALSE
  )

  treatment <- dplyr::filter(
    coefficient_result$tidy_result,
    .data$term == "treatmentB"
  )
  interaction <- dplyr::filter(
    coefficient_result$tidy_result,
    .data$term == "treatmentB:timeT2"
  )

  expect_identical(
    unique(contrast_result$tidy_result$term),
    "B_at_T2_vs_A_at_T2"
  )
  expect_equal(
    contrast_result$tidy_result$estimate,
    treatment$estimate + interaction$estimate
  )
  expect_identical(
    colnames(contrast_result$raw_result$contrast_matrix),
    "B_at_T2_vs_A_at_T2"
  )
})

test_that("gly_linear_model adjusts p-values within each term", {
  exp <- make_linear_model_exp()
  result <- gly_linear_model(
    exp,
    ~ treatment * time + age,
    add_info = FALSE
  )

  expected <- result$tidy_result |>
    dplyr::group_by(.data$term) |>
    dplyr::mutate(expected = stats::p.adjust(.data$p_val, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::pull("expected")
  expect_equal(result$tidy_result$p_adj, expected)

  unadjusted <- gly_linear_model(
    exp,
    ~ treatment * time + age,
    p_adj_method = NULL,
    add_info = FALSE
  )
  expect_false("p_adj" %in% colnames(unadjusted$tidy_result))
})

test_that("gly_linear_model preserves all-NA variables", {
  exp <- make_linear_model_exp()
  SummarizedExperiment::assay(exp)["V1", ] <- NA_real_

  result <- gly_linear_model(exp, ~ treatment + age, add_info = FALSE)
  failed <- dplyr::filter(result$tidy_result, .data$variable == "V1")

  expect_identical(unique(failed$term), c("treatmentB", "age"))
  expect_true(all(is.na(failed$estimate)))
  expect_true(all(is.na(failed$p_val)))
})

test_that("gly_linear_model supports legacy experiments", {
  se <- make_linear_model_exp()
  exp <- suppressWarnings(glyexp::experiment(
    expr_mat = SummarizedExperiment::assay(se),
    sample_info = as.data.frame(SummarizedExperiment::colData(se)),
    var_info = as.data.frame(SummarizedExperiment::rowData(se)),
    exp_type = "others"
  ))

  suppressWarnings(
    result <- suppressMessages(
      gly_linear_model(exp, ~ treatment + age, add_info = FALSE)
    )
  )
  expect_s3_class(result, "glystats_linear_model_res")
})

test_that("gly_linear_model validates formulas and designs", {
  exp <- make_linear_model_exp()

  expect_snapshot(
    error = TRUE,
    gly_linear_model(exp, expression ~ treatment)
  )
  expect_snapshot(
    error = TRUE,
    gly_linear_model(exp, ~missing_predictor)
  )
  expect_snapshot(
    error = TRUE,
    gly_linear_model(exp, ~1)
  )

  missing_exp <- exp
  SummarizedExperiment::colData(missing_exp)$age[1] <- NA_real_
  expect_snapshot(
    error = TRUE,
    gly_linear_model(missing_exp, ~age)
  )

  SummarizedExperiment::colData(exp)$duplicate <-
    SummarizedExperiment::colData(exp)$treatment
  expect_snapshot(
    error = TRUE,
    gly_linear_model(exp, ~ treatment + duplicate)
  )
})

test_that("gly_linear_model validates contrasts", {
  exp <- make_linear_model_exp()

  expect_snapshot(
    error = TRUE,
    gly_linear_model(exp, ~treatment, contrasts = "treatmentB")
  )
  expect_snapshot(
    error = TRUE,
    gly_linear_model(
      exp,
      ~treatment,
      contrasts = c(test = "treatmentB", test = "-treatmentB")
    )
  )
  expect_snapshot(
    error = TRUE,
    gly_linear_model(
      exp,
      ~ treatment * time,
      contrasts = c(test = "unknown_coefficient")
    )
  )
  expect_snapshot(
    error = TRUE,
    gly_linear_model(
      exp,
      ~ treatment * time,
      contrasts = c(test = "treatmentB + treatmentB:timeT2")
    )
  )
})

test_that("filter_sig_vars filters linear model results by term", {
  exp <- make_linear_model_exp()
  result <- gly_linear_model(
    exp,
    ~ treatment * time + age,
    add_info = FALSE
  )
  result$tidy_result$p_adj <- 1
  result$tidy_result$p_val <- 1
  result$tidy_result$p_adj[
    result$tidy_result$variable == "V1" &
      result$tidy_result$term == "treatmentB"
  ] <- 0.001
  result$tidy_result$p_adj[
    result$tidy_result$variable == "V2" &
      result$tidy_result$term == "timeT2"
  ] <- 0.001

  any_term <- filter_sig_vars(exp, result)
  treatment <- filter_sig_vars(exp, result, term = "treatmentB")

  expect_setequal(rownames(any_term), c("V1", "V2"))
  expect_identical(rownames(treatment), "V1")

  unadjusted <- gly_linear_model(
    exp,
    ~treatment,
    p_adj_method = NULL,
    add_info = FALSE
  )
  unadjusted$tidy_result$p_val <- 1
  unadjusted$tidy_result$p_val[1] <- 0.001
  raw_filtered <- filter_sig_vars(
    exp,
    unadjusted,
    p_adj_cutoff = NULL,
    p_val_cutoff = 0.01
  )
  expect_identical(rownames(raw_filtered), "V1")
})

test_that("filter_sig_vars validates linear model filtering", {
  exp <- make_linear_model_exp()
  result <- gly_linear_model(exp, ~treatment, add_info = FALSE)

  expect_snapshot(
    error = TRUE,
    filter_sig_vars(exp, result, term = "unknown")
  )
  expect_snapshot(
    error = TRUE,
    filter_sig_vars(exp, result, fc_cutoff = 2)
  )

  unadjusted <- gly_linear_model(
    exp,
    ~treatment,
    p_adj_method = NULL,
    add_info = FALSE
  )
  expect_snapshot(
    error = TRUE,
    filter_sig_vars(exp, unadjusted)
  )
})
