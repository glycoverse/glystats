identifier_test_matrix <- function() {
  ids <- c(
    "ordinary",
    "sample",
    "group",
    "subject",
    "variable",
    "log_value",
    "age"
  )
  values <- outer(seq_along(ids), c(1, 3, 2, 5, 4, 8, 6, 10), "+")
  dimnames(values) <- list(ids, paste0("s", seq_len(ncol(values))))
  values
}

test_that("two-group tests preserve identifiers matching analysis columns", {
  values <- identifier_test_matrix()
  groups <- factor(rep(c("A", "B"), each = 4))
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(
      group = groups,
      person = rep(letters[1:4], 2)
    )
  )

  for (paired in c(FALSE, TRUE)) {
    for (method in c("ttest", "wilcox")) {
      analyze <- if (method == "ttest") gly_ttest else gly_wilcox
      reference <- if (method == "ttest") stats::t.test else stats::wilcox.test
      result <- suppressWarnings(suppressMessages(analyze(
        exp,
        subject_col = if (paired) "person" else NULL,
        log_transform = FALSE,
        add_info = FALSE,
        p_adj_method = NULL
      )))

      expect_setequal(result$tidy_result$.variable, rownames(values))
      expect_setequal(names(result$raw_result), rownames(values))
      expect_equal(nrow(result$tidy_result), nrow(values))
      for (id in rownames(values)) {
        expected <- suppressWarnings(reference(
          values[id, groups == "B"],
          values[id, groups == "A"],
          paired = paired
        ))
        expect_equal(result$raw_result[[id]]$p.value, expected$p.value)
        expect_equal(
          result$tidy_result$p_val[result$tidy_result$.variable == id],
          expected$p.value
        )
      }
    }
  }
})

test_that("long analysis data preserve feature, subject and covariate alignment", {
  values <- identifier_test_matrix()
  groups <- factor(rep(c("A", "B"), each = 4))
  subjects <- factor(rep(letters[1:4], 2))
  ages <- c(20, 30, 40, 50, 21, 31, 41, 51)
  data <- .prepare_multi_group_data(
    values,
    groups,
    covariates = data.frame(age = ages),
    subjects = subjects,
    is_logged = TRUE
  )

  expect_equal(nrow(data), length(values))
  expect_setequal(data$variable, rownames(values))
  for (id in rownames(values)) {
    feature <- data[data$variable == id, ]
    expect_identical(feature$sample, colnames(values))
    expect_equal(feature$log_value, unname(values[id, ]))
    expect_identical(feature$group, groups)
    expect_identical(feature$subject, subjects)
    expect_identical(feature$age, ages)
  }
})

test_that("ANOVA retains reserved feature identifiers in paired and unpaired fits", {
  values <- identifier_test_matrix()
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(
      group = factor(rep(c("A", "B"), each = 4)),
      person = rep(letters[1:4], 2)
    )
  )
  for (paired in c(FALSE, TRUE)) {
    result <- suppressWarnings(suppressMessages(gly_anova(
      exp,
      subject_col = if (paired) "person" else NULL,
      log_transform = FALSE,
      add_info = FALSE
    )))
    expect_setequal(result$tidy_result$main_test$.variable, rownames(values))
    expect_setequal(names(result$raw_result$main_test), rownames(values))
    expect_true(all(is.finite(result$tidy_result$main_test$p_val)))
  }
})

test_that("ANCOVA preserves a feature sharing a covariate name", {
  skip_if_not_installed("emmeans")
  values <- identifier_test_matrix()
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(
      group = factor(rep(c("A", "B"), each = 4)),
      age = c(20, 40, 30, 50, 22, 42, 32, 52)
    )
  )
  result <- suppressWarnings(suppressMessages(gly_ancova(
    exp,
    covariate_cols = "age",
    log_transform = FALSE,
    add_info = FALSE
  )))
  expect_setequal(result$tidy_result$main_test$.variable, rownames(values))
  expect_setequal(names(result$raw_result$main_test), rownames(values))
  expect_true(all(is.finite(result$tidy_result$main_test$p_val)))
})


duplicate_identifier_experiment <- function() {
  values <- rbind(
    c(1, 2, 4, 3, 10, 11, 13, 12),
    c(100, 110, 90, 105, 2, 4, 1, 3)
  )
  dimnames(values) <- list(rep("duplicate", 2), paste0("s", seq_len(8)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    colData = S4Vectors::DataFrame(
      group = factor(rep(c("A", "B"), each = 4)),
      person = rep(letters[1:4], 2),
      age = c(20, 40, 30, 50, 22, 42, 32, 52)
    )
  )
}

test_that("the former wide-table preparation also rejected duplicate identifiers", {
  exp <- duplicate_identifier_experiment()
  wide <- as.data.frame(t(SummarizedExperiment::assay(exp)))
  expect_gt(anyDuplicated(names(wide)), 0L)
  expect_error(
    tibble::as_tibble(tibble::rownames_to_column(wide, "sample")),
    "unique"
  )
})

test_that("duplicate SE identifiers keep distinct differential fits", {
  exp <- duplicate_identifier_experiment()
  before <- exp
  reference <- exp
  rownames(reference) <- c("first", "second")
  for (paired in c(FALSE, TRUE)) {
    for (analyze in list(gly_ttest, gly_wilcox, gly_anova, gly_kruskal)) {
      if (identical(analyze, gly_kruskal) && !paired) {
        skip_if_not_installed("FSA")
      }
      run <- function(x) {
        suppressWarnings(suppressMessages(analyze(
          x,
          subject_col = if (paired) "person" else NULL,
          log_transform = FALSE,
          add_info = FALSE
        )))
      }
      result <- run(exp)
      expected <- run(reference)
      tbl <- result$tidy_result
      ref <- expected$tidy_result
      if (is.list(tbl) && !tibble::is_tibble(tbl)) {
        tbl <- tbl$main_test
        ref <- ref$main_test
      }
      expect_equal(nrow(tbl), 2L)
      expect_equal(tbl$p_val, ref$p_val)
      expect_identical(tbl$.variable, result$identifiers$variable$data$id)
      expect_identical(
        result$identifiers$variable$data$name,
        rep("duplicate", 2)
      )
    }
  }
  expect_identical(exp, before)
})

test_that("ANCOVA keeps distinct rows with duplicate SE identifiers", {
  skip_if_not_installed("emmeans")
  exp <- duplicate_identifier_experiment()
  result <- suppressWarnings(suppressMessages(gly_ancova(
    exp,
    covariate_cols = "age",
    log_transform = FALSE,
    add_info = FALSE
  )))
  expect_equal(nrow(result$tidy_result$main_test), 2L)
  expect_identical(
    result$tidy_result$main_test$.variable,
    result$identifiers$variable$data$id
  )
})
