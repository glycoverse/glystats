# Wilcoxon rank-sum test for Differential Expression Analysis

Perform Wilcoxon rank-sum test (Mann-Whitney U test) for glycomics or
glycoproteomics data. The function supports non-parametric comparison of
two groups. P-values are adjusted for multiple testing using the method
specified by `p_adj_method`.

## Usage

``` r
gly_wilcox(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  ref_group = NULL,
  add_info = TRUE,
  ...
)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  or `SummarizedExperiment` object containing an expression matrix and
  sample information.

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

## Value

A list with three elements:

- `tidy_result`: A tibble with Wilcoxon test results containing the
  following columns:

  - `variable`: Variable name

  - `statistic`: Wilcoxon test statistic

  - `p_val`: Raw p-value from Wilcoxon test

  - `effect_size`: Rank-biserial correlation

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

The function performs log2 transformation on the expression data
(log2(x + 1e-6)) before statistical testing. Exactly 2 groups are
required in the grouping variable.

## See also

[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
