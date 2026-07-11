# Two-sample t-test for Differential Expression Analysis

Perform two-sample t-test for glycomics or glycoproteomics data. The
function supports Student's t-test for comparing two groups. P-values
are adjusted for multiple testing using the method specified by
`p_adj_method`.

## Usage

``` r
gly_ttest(
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
  [`stats::t.test()`](https://rdrr.io/r/stats/t.test.html).

## Value

A list with three elements:

- `tidy_result`: A tibble with t-test results containing the following
  columns:

  - `variable`: Variable name

  - `estimate`: Difference in group means (group2 - group1)

  - `estimate1`: Mean of group 1

  - `estimate2`: Mean of group 2

  - `statistic`: t-statistic

  - `p_val`: Raw p-value from t-test

  - `parameter`: Degrees of freedom

  - `conf_low`: Lower bound of 95% confidence interval

  - `conf_high`: Upper bound of 95% confidence interval

  - `effect_size`: Cohen's d

  - `method`: Statistical method used

  - `alternative`: Alternative hypothesis

  - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

  - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean))

- `raw_result`: A list of `t.test` model objects

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_ttest_res` and `glystats_res`.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1e-6)) before statistical testing. Exactly 2 groups are
required in the grouping variable.

## See also

[`stats::t.test()`](https://rdrr.io/r/stats/t.test.html)
