# Wilcoxon Test for Differential Expression Analysis

Perform Wilcoxon rank-sum or signed-rank tests for glycomics or
glycoproteomics data. P-values are adjusted for multiple testing using
the method specified by `p_adj_method`.

## Usage

``` r
gly_wilcox(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  ref_group = NULL,
  add_info = TRUE,
  ...,
  subject_col = NULL,
  log_transform = TRUE
)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment` containing an expression
  matrix and sample information.

- group_col:

  A character string specifying the column name of the grouping variable
  in the sample information. Default is `"group"`.

- p_adj_method:

  A character string specifying the method to adjust p-values. See
  `p.adjust.methods` for available methods. Default is "BH". If NULL, no
  adjustment is performed.

- ref_group:

  A character string specifying the reference group. If NULL (default),
  the first level of the group factor is used as the reference.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned.

- ...:

  Additional arguments passed to
  [`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html).

- subject_col:

  An optional character string naming the subject identifier column in
  sample information. When supplied, samples are matched by subject and
  a paired Wilcoxon signed-rank test is performed.

- log_transform:

  A single non-missing logical value. If `TRUE` (default), transform
  expression values using `log2(x + 1e-6)` before analysis. If `FALSE`,
  use the supplied values unchanged, for example for already
  log-transformed data.

## Value

A list with three elements:

- `tidy_result`: A tibble with Wilcoxon test results containing the
  following columns:

  - `variable`: Variable name

  - `statistic`: Wilcoxon test statistic

  - `p_val`: Raw p-value from Wilcoxon test

  - `effect_size`: Rank-biserial correlation, using matched pairs in
    paired analyses

  - `method`: Statistical method used

  - `alternative`: Alternative hypothesis

  - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

  - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean))
    Additional columns from experiment metadata may be included if
    add_info = TRUE.

- `raw_result`: A list of `wilcox.test` model objects

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_wilcox_res` and `glystats_res`.

## Details

Tests and effect sizes use the selected analysis scale. Reported fold
changes continue to use ratios of group means from the original supplied
values.

By default, expression data are transformed using `log2(x + 1e-6)`
before statistical testing. Set `log_transform = FALSE` to skip this
transformation. Exactly 2 groups are required in the grouping variable.
In paired analyses, only subjects observed in both groups are used and
`effect_size` is the matched-pairs rank-biserial correlation.

## See also

[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
