# Test Correlated Variable Sets

Jointly test correlated variable sets with Hotelling's \\T^2\\
statistic.

## Usage

``` r
gly_set_test(
  exp,
  sets,
  group_col = "group",
  subject_col = NULL,
  method = "hotelling",
  p_adj_method = "BH"
)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment` containing expression data
  and sample information.

- sets:

  A `glystats_correlated_sets` object returned by
  [`gly_correlated_sets()`](https://glycoverse.github.io/glystats/dev/reference/gly_correlated_sets.md),
  or a uniquely named list of character vectors.

- group_col:

  A character string naming the two-level grouping column in sample
  information. Default is `"group"`.

- subject_col:

  An optional character string naming subject identifiers. When
  supplied, one-sample Hotelling tests are applied to paired
  test-minus-reference difference vectors. Otherwise independent
  two-sample Hotelling tests are used.

- method:

  Statistical method. Currently only `"hotelling"` is supported.

- p_adj_method:

  A character string specifying the method used to adjust set-level
  p-values. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  is `"BH"`. If `NULL`, no adjustment is performed.

## Value

A list with classes `glystats_set_test_res` and `glystats_res`
containing:

- `tidy_result$sets`: One row per set with its estimate vector,
  Hotelling statistic, degrees of freedom, Mahalanobis effect size,
  p-value, adjusted p-value, and fit status.

- `tidy_result$members`: One row per set member with its marginal
  estimate and mean within-set correlation.

- `raw_result`: Set definitions, their correlation matrix, and
  individual Hotelling test objects.

- `meta_data`: Metadata copied from `exp`.

## Details

Expression values are transformed using `log2(x + 1e-6)`. The first
factor level of `group_col` is the reference group and the second is the
test group. The estimate is therefore `test - reference`.

Equivalent profiles are included in member-level output but represented
once in the covariance matrix. Sets that are too large for their
residual sample size or have a singular covariance matrix are retained
with `status = "failed"`, `NA` statistics, and an explicit
`failure_reason`.

## See also

[`gly_correlated_sets()`](https://glycoverse.github.io/glystats/dev/reference/gly_correlated_sets.md),
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)

## Examples

``` r
set.seed(1)
expression <- matrix(
  stats::rexp(36),
  nrow = 3,
  dimnames = list(c("A", "B", "C"), paste0("S", 1:12))
)
expression["B", ] <- expression["A", ] * exp(stats::rnorm(12, sd = 0.01))
exp <- SummarizedExperiment::SummarizedExperiment(
  assays = list(expression = expression),
  colData = S4Vectors::DataFrame(
    group = factor(
      rep(c("control", "case"), each = 6),
      levels = c("control", "case")
    )
  )
)
sets <- gly_correlated_sets(exp, threshold = 0.8)
result <- gly_set_test(exp, sets)
#> ℹ Ref Group: "control"
#> ℹ Test Group: "case"
result$tidy_result$sets
#> # A tibble: 1 × 18
#>   set_id n_variables test_dimension variables estimate  statistic   df1   df2
#>   <chr>        <int>          <int> <list>    <list>        <dbl> <int> <dbl>
#> 1 set_1            2              2 <chr [2]> <dbl [2]>      6.25     2     9
#> # ℹ 10 more variables: effect_size <dbl>, p_val <dbl>, n_ref <int>,
#> #   n_test <int>, ref_group <chr>, test_group <chr>, paired <lgl>,
#> #   status <chr>, failure_reason <chr>, p_adj <dbl>
```
