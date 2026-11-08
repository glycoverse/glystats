identity_se <- function() {
  values <- rbind(
    c(1, 2, 4, 3, 10, 11, 13, 12),
    c(100, 110, 90, 105, 2, 4, 1, 3),
    c(2, 3, 1, 5, 4, 2, 5, 3)
  )
  dimnames(values) <- list(c("g1", "g2", "g3"), paste0("s", 1:8))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    rowData = S4Vectors::DataFrame(annotation = c("first", "second", "third")),
    colData = S4Vectors::DataFrame(
      group = factor(rep(c("A", "B"), each = 4)),
      person = rep(letters[1:4], 2),
      age = c(20, 40, 30, 50, 22, 42, 32, 52)
    )
  )
}

quiet_identity <- function(x) suppressWarnings(suppressMessages(x))

test_that("ordinary metadata IDs never override SE axes", {
  se <- identity_se()
  baseline <- gly_pca(se, log_transform = FALSE)
  SummarizedExperiment::rowData(se)$variable <- c("g2", "g1", "g1")
  SummarizedExperiment::colData(se)$sample <- rep("s1", ncol(se))
  result <- gly_pca(se, log_transform = FALSE)
  vars <- result$tidy_result$variables
  samples <- result$tidy_result$samples
  expect_equal(vars$value, baseline$tidy_result$variables$value)
  expect_equal(samples$value, baseline$tidy_result$samples$value)
  expect_identical(vars$variable, baseline$tidy_result$variables$variable)
  expect_identical(samples$sample, baseline$tidy_result$samples$sample)
  expect_equal(nrow(vars), nrow(baseline$tidy_result$variables))
  expect_identical(vars$annotation, c("first", "second", "third")[match(vars$variable, rownames(se))])
  expect_identical(vars$variable_metadata, c("g2", "g1", "g1")[match(vars$variable, rownames(se))])
  expect_true(all(samples$sample_metadata == "s1"))
})

test_that("metadata sample and variable remain usable model columns", {
  se <- identity_se()
  SummarizedExperiment::colData(se)$sample <- SummarizedExperiment::colData(se)$group
  SummarizedExperiment::rowData(se)$variable <- c("x", "x", "y")
  by_group <- quiet_identity(gly_ttest(se, log_transform = FALSE))
  by_sample <- quiet_identity(gly_ttest(se, group_col = "sample", log_transform = FALSE))
  expect_equal(by_group$tidy_result, by_sample$tidy_result)
  skip_if_not_installed("limma")
  lm_group <- quiet_identity(gly_linear_model(se, ~ group, log_transform = FALSE))
  lm_sample <- quiet_identity(gly_linear_model(se, ~ sample, log_transform = FALSE))
  expect_equal(lm_group$tidy_result$p_val, lm_sample$tidy_result$p_val)
  expect_identical(.get_var_info(se)$variable, c("x", "x", "y"))
})

test_that("annotation collisions preserve data and cannot add spurious joins", {
  se <- identity_se()
  SummarizedExperiment::rowData(se)$variable <- c("x", "x", "y")
  SummarizedExperiment::rowData(se)$variable_metadata <- 1:3
  SummarizedExperiment::rowData(se)$sample <- c("not", "sample", "ids")
  SummarizedExperiment::rowData(se)$p_val <- c(99, 98, 97)
  result <- quiet_identity(gly_ttest(se, log_transform = FALSE))$tidy_result
  expect_equal(nrow(result), 3L)
  expect_true(all(result$p_val >= 0 & result$p_val <= 1))
  expect_identical(result$p_val_metadata, c(99, 98, 97)[match(result$variable, rownames(se))])
  expect_identical(result$variable_metadata, (1:3)[match(result$variable, rownames(se))])
  expect_identical(result$variable_metadata_metadata, c("x", "x", "y")[match(result$variable, rownames(se))])
  expect_identical(result$sample, c("not", "sample", "ids")[match(result$variable, rownames(se))])
})

test_that("NULL, duplicate and empty names use consistent position IDs", {
  for (kind in c("NULL", "duplicate", "empty")) {
    se <- identity_se()
    if (kind == "NULL") dimnames(se) <- NULL
    if (kind == "duplicate") dimnames(se) <- list(rep("g", 3), rep("s", 8))
    if (kind == "empty") dimnames(se) <- list(c("", "g", "g"), rep("", 8))
    before <- se
    result <- quiet_identity(gly_ttest(se, log_transform = FALSE))
    map <- result$identifiers$variable
    expect_true(map$generated)
    expect_identical(map$data$id, paste0(".glystats_variable_", 1:3))
    expect_equal(result$tidy_result$p_val, quiet_identity(gly_ttest(identity_se(), log_transform = FALSE))$tidy_result$p_val)
    expect_identical(result$tidy_result$variable, map$data$id)
    expect_identical(result$tidy_result$annotation, c("first", "second", "third"))
    expect_identical(map$names_missing, kind == "NULL")
    expect_identical(se, before)
    no_info <- quiet_identity(gly_ttest(se, log_transform = FALSE, add_info = FALSE))
    expect_identical(no_info$identifiers, result$identifiers)
  }
})

test_that("paired post-hoc effects use positions even for duplicate samples", {
  se <- identity_se()
  ref <- quiet_identity(gly_anova(se, subject_col = "person", log_transform = FALSE, add_info = FALSE))
  colnames(se) <- rep("same", ncol(se))
  got <- quiet_identity(gly_anova(se, subject_col = "person", log_transform = FALSE, add_info = FALSE))
  expect_equal(got$tidy_result, ref$tidy_result)
})

test_that("position-keyed filtering selects original rows and rejects stale input", {
  se <- identity_se()
  rownames(se) <- c("same", "same", "")
  res <- quiet_identity(gly_ttest(se, log_transform = FALSE))
  res$tidy_result$p_adj <- c(.01, .9, .9)
  filtered <- filter_sig_vars(se, res, fc_cutoff = NULL)
  expect_identical(filtered, se[1, , drop = FALSE])
  expect_error(filter_sig_vars(se[c(2, 1, 3), ], res), "same unmodified experiment")
  changed <- se
  SummarizedExperiment::assay(changed)[1, 1] <- 999
  expect_error(filter_sig_vars(changed, res), "same unmodified experiment")
  old_res <- res
  old_res$identifiers <- NULL
  expect_error(filter_sig_vars(se, old_res), "same unmodified experiment")
  rownames(se) <- NULL
  res <- quiet_identity(gly_ttest(se, log_transform = FALSE))
  res$tidy_result$p_adj <- c(.9, .01, .9)
  expect_identical(filter_sig_vars(se, res, fc_cutoff = NULL), se[2, , drop = FALSE])
})

test_that("named filtering remains compatible with reordered SE inputs", {
  se <- identity_se()
  res <- quiet_identity(gly_ttest(se, log_transform = FALSE))
  res$tidy_result$p_adj <- c(.01, .9, .9)
  expect_identical(filter_sig_vars(se[3:1, ], res, fc_cutoff = NULL), se[1, , drop = FALSE])
})

test_that("missing variables are restored using unique position IDs", {
  se <- identity_se()
  rownames(se) <- rep("same", nrow(se))
  SummarizedExperiment::assay(se)[2, ] <- NA_real_
  result <- quiet_identity(gly_pca(se, log_transform = FALSE))
  vars <- result$tidy_result$variables
  expect_setequal(vars$variable, result$identifiers$variable$data$id)
  expect_true(all(is.na(vars$value[vars$variable == ".glystats_variable_2"])))
  expect_true(all(vars$annotation[vars$variable == ".glystats_variable_2"] == "second"))
})

test_that("identity helpers support empty axes without inventing rows", {
  se <- identity_se()[FALSE, FALSE]
  dimnames(se) <- NULL
  expect_identical(.se_axis_ids(se, "variable"), character())
  expect_identical(.se_axis_ids(se, "sample"), character())
})

test_that("analysis families share position IDs and metadata alignment", {
  se <- identity_se()
  dimnames(se) <- list(rep("same", nrow(se)), rep("", ncol(se)))
  SummarizedExperiment::colData(se)$time <- c(8, 3, 5, 9, 4, 7, 6, 2)
  SummarizedExperiment::colData(se)$event <- c(1, 1, 0, 1, 0, 1, 1, 0)
  ref <- se
  dimnames(ref) <- list(.se_axis_ids(se, "variable"), .se_axis_ids(se, "sample"))
  methods <- list(
    pca = function(x) gly_pca(x, log_transform = FALSE),
    fold_change = function(x) gly_fold_change(x, log_transform = FALSE),
    cor = function(x) gly_cor(x, log_transform = FALSE),
    kmeans = function(x) gly_kmeans(x, centers = 2, log_transform = FALSE),
    cox = function(x) gly_cox(x, log_transform = FALSE),
    set_test = function(x) gly_set_test(x, log_transform = FALSE),
    custom_sets = function(x) gly_set_test(x, sets = list(a = c(".glystats_variable_1", ".glystats_variable_2")), log_transform = FALSE)
  )
  if (requireNamespace("limma", quietly = TRUE)) {
    methods$limma <- function(x) gly_limma(x, log_transform = FALSE)
    methods$linear_model <- function(x) gly_linear_model(x, ~ group, log_transform = FALSE)
  }
  if (requireNamespace("pROC", quietly = TRUE)) methods$roc <- function(x) gly_roc(x)
  if (requireNamespace("ggdendro", quietly = TRUE)) methods$hclust <- function(x) gly_hclust(x, k_values = 2, log_transform = FALSE)
  if (requireNamespace("Rtsne", quietly = TRUE)) methods$tsne <- function(x) gly_tsne(x, perplexity = 2, max_iter = 250, log_transform = FALSE)
  if (requireNamespace("uwot", quietly = TRUE)) methods$umap <- function(x) gly_umap(x, n_neighbors = 3, n_epochs = 20, n_threads = 1, log_transform = FALSE)
  for (method in methods) {
    set.seed(42)
    got <- quiet_identity(method(se))
    set.seed(42)
    expected <- quiet_identity(method(ref))
    expect_equal(got$tidy_result, expected$tidy_result)
    expect_true(got$identifiers$variable$generated)
    expect_true(got$identifiers$sample$generated)
    expect_identical(got$identifiers$variable$data$name, rownames(se))
    expect_identical(got$identifiers$sample$data$name, colnames(se))
  }
})

test_that("supervised models keep distinct unnamed samples and variables", {
  skip_if_not_installed("ropls")
  se <- exp_topliss_valid()
  dimnames(se) <- NULL
  ref <- se
  dimnames(ref) <- list(.se_axis_ids(se, "variable"), .se_axis_ids(se, "sample"))
  for (method in list(gly_plsda, gly_oplsda)) {
    set.seed(42)
    got <- quiet_identity(method(se, permI = 0, log_transform = FALSE))
    set.seed(42)
    expected <- quiet_identity(method(ref, permI = 0, log_transform = FALSE))
    expect_equal(got$tidy_result, expected$tidy_result)
    expect_true(got$identifiers$variable$names_missing)
    expect_true(got$identifiers$sample$names_missing)
  }
})

test_that("top-level names override stale names stored in the assay", {
  se <- identity_se()
  rownames(se) <- paste0("new", 1:3)
  colnames(se) <- paste0("new_sample", 1:8)
  expect_false(identical(rownames(SummarizedExperiment::assay(se, withDimnames = FALSE)), rownames(se)))
  result <- gly_pca(se, log_transform = FALSE)
  expect_setequal(result$tidy_result$variables$variable, rownames(se))
  expect_setequal(result$tidy_result$samples$sample, colnames(se))
})

test_that("linear model filters use the stored position mapping", {
  skip_if_not_installed("limma")
  se <- identity_se()
  rownames(se) <- rep("same", nrow(se))
  for (method in list(
    function(x) gly_limma(x, log_transform = FALSE),
    function(x) gly_linear_model(x, ~ group, log_transform = FALSE)
  )) {
    result <- quiet_identity(method(se))
    result$tidy_result$p_adj <- ifelse(result$tidy_result$variable == ".glystats_variable_2", .01, .9)
    expect_identical(filter_sig_vars(se, result, fc_cutoff = NULL), se[2, , drop = FALSE])
    expect_error(filter_sig_vars(se[3:1, ], result), "same unmodified experiment")
  }
})

test_that("generated IDs cannot accidentally match legacy experiment IDs", {
  se <- identity_se()
  rownames(se) <- NULL
  result <- quiet_identity(gly_ttest(se, log_transform = FALSE))
  mat <- .get_expr_mat(se)
  legacy <- suppressWarnings(glyexp::experiment(
    mat,
    sample_info = tibble::tibble(sample = colnames(mat), group = c(rep("A", 4), rep("B", 4))),
    var_info = tibble::tibble(variable = rownames(mat)),
    exp_type = "others"
  ))
  expect_error(filter_sig_vars(legacy, result), "same unmodified experiment")
})
