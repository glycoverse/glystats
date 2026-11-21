skip_if_not_installed("emmeans")

ancova_test_experiment <- function(values, groups, ...) {
  colnames(values) <- paste0("S", seq_len(ncol(values)))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expression = values),
    colData = S4Vectors::DataFrame(group = groups, ...)
  )
}

test_that("ANCOVA excludes the variation explained by covariates from the group test", {
  group <- factor(rep(c("A", "B"), each = 4))
  z <- c(0:3, 3:6)
  e <- rep(c(1, -1, -1, 1), 2)
  y <- z + e
  d <- data.frame(y, group, z)
  first <- stats::anova(stats::lm(y ~ group + z, data = d))
  expect_equal(first["group", "F value"], 11.25)
  expect_lt(first["group", "Pr(>F)"], 0.05)

  exp <- ancova_test_experiment(rbind(null = y), group, z = z)
  result <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = "z",
    log_transform = FALSE
  ))
  main <- result$tidy_result$main_test
  expect_equal(main$.variable, "null")
  # The errors are orthogonal to both group and z, so the conditional
  # group effect is exactly zero. Subtracting nested-model residual sums of
  # squares here can instead produce tiny negative values on some builds.
  expect_equal(main$df, 1)
  expect_equal(main$sumsq, 0, tolerance = 1e-12)
  expect_equal(main$meansq, 0, tolerance = 1e-12)
  expect_equal(main$statistic, 0, tolerance = 1e-12)
  expect_equal(main$p_val, 1, tolerance = 1e-12)
  expect_equal(main$p_adj, 1, tolerance = 1e-12)
  expect_true(is.na(main$post_hoc))
  expect_equal(nrow(result$tidy_result$post_hoc_test), 0L)
  expect_length(result$raw_result$post_hoc_test, 0L)
  expect_s3_class(result$raw_result$main_test[["null"]], "aov")
})

test_that("ANCOVA uses the same complete weighted observations for both tests", {
  group <- factor(rep(c("A", "B"), each = 6))
  z <- c(0, 1, NA_real_, 3, 4, 5, 2, 3, 4, 5, 6, 7)
  e <- c(0.2, -0.1, 0.3, -0.3, -0.2, 0.1, 0.1, -0.3, 0.2, -0.1, 0.4, -0.2)
  weights <- c(1, 3, 1, 2, 1, 1, 2, 1, 2, 4, 1, 0)
  y <- 4 + ifelse(group == "B", 3, 0) + 0.8 * z + e
  y[7] <- NA_real_
  d <- data.frame(y, group, z, w = weights)
  retained <- seq_along(y) != 2L
  d <- d[retained & stats::complete.cases(d), ]
  full <- stats::aov(y ~ z + group, data = d, weights = w)
  reduced <- stats::lm(y ~ z, data = d, weights = w)
  expected <- stats::anova(reduced, full)
  exp <- ancova_test_experiment(rbind(signal = y), group, z = z)

  result <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = "z",
    weights = weights,
    subset = retained,
    log_transform = FALSE
  ))
  main <- result$tidy_result$main_test
  expect_equal(main$statistic, expected$F[2L])
  expect_equal(main$sumsq, expected$"Sum of Sq"[2L])
  expect_equal(main$p_val, expected$"Pr(>F)"[2L])
  expect_identical(main$post_hoc, "A_vs_B")
  raw_main <- result$raw_result$main_test$signal
  expect_equal(stats::nobs(raw_main), stats::nobs(full))
  expect_equal(unname(stats::coef(raw_main)), unname(stats::coef(full)))

  emm <- emmeans::emmeans(full, specs = "group")
  pair <- as.data.frame(summary(emmeans::contrast(
    emm,
    method = "pairwise",
    adjust = "tukey"
  )))
  expect_equal(result$raw_result$post_hoc_test$signal$p.value, pair$p.value)
  expect_equal(result$raw_result$post_hoc_test$signal$estimate, pair$estimate)
  expect_equal(result$tidy_result$post_hoc_test$p_adj, pair$p.value)
})

test_that("ANCOVA gates post-hoc comparisons using the same reported adjusted p-values", {
  group <- factor(rep(c("A", "B"), each = 4))
  z <- rep(0:3, 2)
  e <- rep(c(1, -1, -1, 1), 2)
  # Residual MSE is 8/5, with group design information of 2.
  d <- sqrt(stats::qf(0.96, 1, 5) * (8 / 5) / 2)
  values <- rbind(
    borderline = 8 + z + e + d * (group == "B"),
    null = 8 + z + e
  )
  exp <- ancova_test_experiment(values, group, z = z)
  with_bh <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = "z",
    log_transform = FALSE
  ))
  main <- with_bh$tidy_result$main_test
  expect_equal(main$p_val, c(.04, 1))
  expect_equal(main$p_adj, stats::p.adjust(main$p_val, method = "BH"))
  expect_equal(main$p_adj, c(.08, 1))
  expect_length(with_bh$raw_result$post_hoc_test, 0L)

  unadjusted <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = "z",
    p_adj_method = NULL,
    log_transform = FALSE
  ))
  expect_false("p_adj" %in% names(unadjusted$tidy_result$main_test))
  expect_equal(unadjusted$tidy_result$main_test$p_val, main$p_val)
  expect_named(unadjusted$raw_result$post_hoc_test, "borderline")
  expect_equal(
    unadjusted$tidy_result$main_test$post_hoc,
    c("A_vs_B", NA_character_)
  )
})

test_that("ANCOVA retains non-estimable rows and only rejects aliasing of the group term", {
  group <- factor(rep(c("A", "B", "C"), each = 4))
  z <- rep(0:3, 3)
  e <- rep(c(1, -1, -1, 1), 3)
  y <- 8 + z + e + as.integer(group)
  d <- data.frame(y, group, z)
  saturated <- y
  saturated[-c(1, 2, 5, 9)] <- NA_real_
  exp <- ancova_test_experiment(
    rbind(estimable = y, saturated = saturated),
    group,
    z = z,
    duplicate_z = 2 * z
  )
  out <- suppressMessages(gly_ancova(
    exp,
    covariate_cols = c("z", "duplicate_z"),
    log_transform = FALSE
  ))
  expected <- stats::anova(
    stats::lm(y ~ z, data = d),
    stats::lm(y ~ z + group, data = d)
  )
  main <- out$tidy_result$main_test
  expect_equal(main$.variable, c("estimable", "saturated"))
  expect_equal(main$df, c(2, NA_real_))
  expect_equal(main$p_val, c(expected$"Pr(>F)"[2L], NA_real_))
  expect_equal(main$statistic, c(expected$F[2L], NA_real_))
  expect_true(is.na(main$post_hoc[2L]))
  expect_s3_class(out$raw_result$main_test$saturated, "aov")

  # One of the two nominal group contrasts is fully confounded.
  alias <- ancova_test_experiment(
    rbind(aliased = y),
    group,
    z = z,
    group_b = as.integer(group == "B")
  )
  out <- suppressMessages(gly_ancova(
    alias,
    covariate_cols = c("z", "group_b"),
    log_transform = FALSE
  ))
  expect_equal(out$tidy_result$main_test$.variable, "aliased")
  expect_true(is.na(out$tidy_result$main_test$df))
  expect_true(is.na(out$tidy_result$main_test$p_val))
  expect_length(out$raw_result$post_hoc_test, 0L)
  expect_s3_class(out$raw_result$main_test$aliased, "aov")
})
