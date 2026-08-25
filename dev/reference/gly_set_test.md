# Construct and Test Correlated Variable Sets

Construct and jointly test correlated variable sets with Hotelling's
\\T^2\\ statistic.

## Usage

``` r
gly_set_test(
  exp,
  sets = NULL,
  group_col = "group",
  subject_col = NULL,
  method = "hotelling",
  p_adj_method = "BH",
  threshold = 0.9,
  correlation = "pearson",
  clustering = "connected",
  within = NULL
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

  A uniquely named list of character vectors defining custom sets. When
  `NULL` (the default), correlated sets are constructed automatically
  from `exp` using `threshold`, `correlation`, `clustering`, and
  `within`.

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

- threshold:

  A numeric correlation threshold between 0 and 1 used when `sets` is
  `NULL`. Variables are linked only when their correlation is strictly
  greater than this value. Default is `0.9`.

- correlation:

  Correlation coefficient used when `sets` is `NULL`. Either `"pearson"`
  (default) or `"spearman"`.

- clustering:

  How correlated variables are grouped when `sets` is `NULL`.
  `"connected"` (default) allows transitive chaining; `"complete"`
  requires every pair in a set to exceed `threshold`.

- within:

  Optional variable-information columns defining strata within which
  automatic sets are constructed, such as
  `c("protein", "protein_site")`.

## Value

A list with classes `glystats_set_test_res` and `glystats_res`
containing:

- `tidy_result$sets`: One row per set with estimates for every member,
  the effective covariance rank in `test_dimension` (`NA` when rejected
  before covariance estimation), Hotelling statistic, degrees of
  freedom, Mahalanobis effect size, p-value, adjusted p-value, and fit
  status.

- `tidy_result$members`: One row per set member with its marginal
  estimate and mean within-set correlation.

- `raw_result`: Set definitions, their correlation matrix, individual
  Hotelling test objects, and automatic set-construction details (when
  `sets` is `NULL`).

- `meta_data`: Metadata copied from `exp`.

## Details

Expression values are transformed using `log2(x + 1e-6)`. The first
factor level of `group_col` is the reference group and the second is the
test group. The estimate is therefore `test - reference`.

Every usable variable is assigned to a set, including a one-variable set
when it is not correlated above `threshold` with another variable.
Structurally rank-deficient covariance matrices are standardized to unit
member variance and tested in their nonredundant subspace using a
Moore-Penrose inverse when the original set dimension satisfies the
classical sample-size requirement and the observed mean contrast lies in
the estimable subspace. Other ineligible sets are retained with
`status = "failed"`, `NA` statistics, and an explicit `failure_reason`.

## See also

[`stats::cor()`](https://rdrr.io/r/stats/cor.html),
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html),
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
result <- gly_set_test(exp, threshold = 0.8)
#> ℹ Ref Group: "control"
#> ℹ Test Group: "case"
result$tidy_result$sets
#> # A tibble: 2 × 18
#>   set_id n_variables test_dimension variables estimate  statistic   df1   df2
#>   <chr>        <int>          <int> <list>    <list>        <dbl> <int> <dbl>
#> 1 set_1            2              2 <chr [2]> <dbl [2]>   6.25        2     9
#> 2 set_2            1              1 <chr [1]> <dbl [1]>   0.00412     1    10
#> # ℹ 10 more variables: effect_size <dbl>, p_val <dbl>, n_ref <int>,
#> #   n_test <int>, ref_group <chr>, test_group <chr>, paired <lgl>,
#> #   status <chr>, failure_reason <chr>, p_adj <dbl>
```
