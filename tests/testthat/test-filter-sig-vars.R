# ===== ANOVA =====
small_exp <- function() {
  test_gp_exp |>
    glyexp::slice_head_var(n = 10) |>
    glyexp::mutate_obs(
      group = factor(.data$group, levels = c("H", "M", "Y", "C"))
    ) |>
    as_test_se()
}

# ----- ANOVA main test results -----
# > result$tidy_result$main_test |> dplyr::select(variable, p_val, p_adj)
# # A tibble: 10 × 3
#    variable       p_val      p_adj
#    <chr>          <dbl>      <dbl>
#  1 V1       0.0145      0.0207
#  2 V2       0.0000424   0.000141
#  3 V3       0.000105    0.000262
#  4 V4       0.000000285 0.00000285
#  5 V5       0.156       0.173
#  6 V6       0.000000600 0.00000300
#  7 V7       0.0127      0.0207
#  8 V8       0.00316     0.00632
#  9 V9       0.0463      0.0579
# 10 V10      0.432       0.432

test_that("filter_sig_vars works for ANOVA with default settings", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  result$tidy_result$main_test$p_adj <- c(rep(0.01, 5), rep(1, 5))
  sig_exp <- filter_sig_vars(exp, result)
  expect_equal(nrow(sig_exp), 5)
  expect_s4_class(sig_exp, "GlycoproteomicSE")
})

test_that("filter_sig_vars works for ANOVA with p_val_cutoff", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  result$tidy_result$main_test$p_val <- c(rep(0.01, 5), rep(1, 5))
  sig_exp <- filter_sig_vars(
    exp,
    result,
    p_val_cutoff = 0.05,
    p_adj_cutoff = NULL
  )
  expect_equal(nrow(sig_exp), 5)
})

test_that("filter_sig_vars fails for ANOVA with fc_cutoff on main test", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(filter_sig_vars(exp, result, fc_cutoff = 2), error = TRUE)
})

test_that("filter_sig_vars fails for ANOVA with comparison on main test", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(
    filter_sig_vars(exp, result, comparison = "C_vs_H"),
    error = TRUE
  )
})

test_that("filter_sig_vars fails for ANOVA with both p_val_cutoff and p_adj_cutoff", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(
    filter_sig_vars(exp, result, p_val_cutoff = 0.05, p_adj_cutoff = 0.05),
    error = TRUE
  )
})

test_that("filter_sig_vars fails for ANOVA with neither p_val_cutoff, p_adj_cutoff", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(
    filter_sig_vars(exp, result, p_val_cutoff = NULL, p_adj_cutoff = NULL),
    error = TRUE
  )
})

test_that("filter_sig_vars works for ANOVA on post-hoc results", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  post_hoc_res <- result$tidy_result$post_hoc_test |>
    dplyr::mutate(
      comparison = paste0(.data$ref_group, "_vs_", .data$test_group)
    ) |>
    dplyr::arrange(.data$comparison) |>
    dplyr::mutate(row_num = dplyr::row_number()) |>
    dplyr::mutate(
      p_adj = dplyr::if_else(
        .data$row_num %% 2 == 0 & .data$comparison == "H_vs_C",
        0.01,
        1
      )
    )
  result$tidy_result$post_hoc_test <- post_hoc_res
  sig_exp <- filter_sig_vars(
    exp,
    result,
    on = "post_hoc_test",
    comparison = "H_vs_C"
  )
  expect_equal(nrow(sig_exp), 3)
})

test_that("filter_sig_vars fails for ANOVA on post-hoc results without comparison", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(
    filter_sig_vars(exp, result, on = "post_hoc_test"),
    error = TRUE
  )
})

test_that("filter_sig_vars fails for ANOVA with invalid comparison format", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(
    filter_sig_vars(
      exp,
      result,
      on = "post_hoc_test",
      comparison = "H_vs_C_vs_M"
    ),
    error = TRUE
  )
})

test_that("filter_sig_vars fails for ANOVA with non-existent comparison", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(
    filter_sig_vars(exp, result, on = "post_hoc_test", comparison = "H_vs_D"),
    error = TRUE
  )
})

test_that("filter_sig_vars fails for ANOVA with comparison of wrong order", {
  exp <- small_exp()
  result <- suppressMessages(gly_anova(exp))
  expect_snapshot(
    filter_sig_vars(exp, result, on = "post_hoc_test", comparison = "C_vs_H"),
    error = TRUE
  )
})

test_that("filter_sig_vars works for Kruskal-Wallis test", {
  exp <- small_exp()
  result <- suppressMessages(gly_kruskal(exp))
  result$tidy_result$main_test$p_adj <- c(rep(0.01, 5), rep(1, 5))
  sig_exp <- filter_sig_vars(exp, result)
  expect_equal(nrow(sig_exp), 5)
})


# ===== t-test =====
small_exp2 <- function() {
  test_gp_exp |>
    glyexp::slice_head_var(n = 10) |>
    glyexp::filter_obs(.data$group %in% c("C", "H")) |>
    glyexp::mutate_obs(group = factor(.data$group, levels = c("H", "C"))) |>
    as_test_se()
}

# ----- t-test results -----
# > result$tidy_result[c("variable", "p_val", "p_adj", "log2fc")]
# # A tibble: 10 × 4
#    variable     p_val    p_adj log2fc
#    <chr>        <dbl>    <dbl>  <dbl>
#  1 V1       0.960     0.992    -0.111
#  2 V10      0.992     0.992    -0.930
#  3 V2       0.00308   0.0103   10.9
#  4 V3       0.0221    0.0368    7.28
#  5 V4       0.00566   0.0135    8.36
#  6 V5       0.0985    0.123    -1.94
#  7 V6       0.00677   0.0135   14.1
#  8 V7       0.0661    0.0944    4.74
#  9 V8       0.0000361 0.000361 10.7
# 10 V9       0.000696  0.00348   8.29

test_that("filter_sig_vars works for t-test with default settings", {
  exp <- small_exp2()
  result <- suppressMessages(gly_ttest(exp))
  result$tidy_result$p_adj <- c(rep(0.01, 5), rep(1, 5))
  sig_exp <- filter_sig_vars(exp, result)
  expect_equal(nrow(sig_exp), 4)
})

# ===== limma =====
test_that("filter_sig_vars works for limma with default settings for 2 groups", {
  exp <- small_exp2()
  result <- suppressMessages(gly_limma(exp))
  result$tidy_result$p_adj <- c(rep(0.01, 5), rep(1, 5))
  result$tidy_result$log2fc <- c(rep(c(1, 2), 5))
  sig_exp <- filter_sig_vars(exp, result)
  expect_equal(nrow(sig_exp), 2)
})

test_that("filter_sig_vars works for limma for multiple groups", {
  exp <- small_exp()
  result <- suppressMessages(gly_limma(exp))
  sig_exp <- filter_sig_vars(exp, result, comparison = "H_vs_C")
  expect_equal(nrow(sig_exp), 7)
})
