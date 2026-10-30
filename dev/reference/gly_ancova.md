# Analysis of Covariance (ANCOVA) for Differential Expression Analysis

Perform ANCOVA for glycomics or glycoproteomics data with sample-level
covariates. The function supports parametric comparison of multiple
groups while adjusting for covariates. For significant results, emmeans
post-hoc comparisons (Tukey adjustment) are automatically performed.
P-values are adjusted for multiple testing using the method specified by
`p_adj_method`.

## Usage

``` r
gly_ancova(
  exp,
  group_col = "group",
  covariate_cols = NULL,
  p_adj_method = "BH",
  add_info = TRUE,
  ...,
  log_transform = NULL
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

- covariate_cols:

  A character vector specifying column names in sample information to
  include as covariates. At least one covariate must be provided.

- p_adj_method:

  A character string specifying the method to adjust p-values. See
  `p.adjust.methods` for available methods. Default is "BH". If NULL, no
  adjustment is performed.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned.

- ...:

  Additional arguments passed to
  [`stats::aov()`](https://rdrr.io/r/stats/aov.html).

- log_transform:

  `NULL` (default) uses `metadata(exp)$scale` to avoid logging data
  twice. `TRUE` forces the method's log transformation on the supplied
  values, even if they are already logged; `FALSE` uses the supplied
  values unchanged. This argument must be named when it follows `...`.

## Value

A list containing four elements:

- `tidy_result`: A list containing:

  - `main_test`: A tibble with ANCOVA results containing the following
    columns:

    - `variable`: Variable name

    - `term`: ANCOVA term (usually "group")

    - `df`: Degrees of freedom

    - `sumsq`: Sum of squares

    - `meansq`: Mean squares

    - `statistic`: F-statistic

    - `p_val`: Raw p-value from ANCOVA

    - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

    - `post_hoc`: Significant group pairs from post-hoc test, in the
      format of "ref_vs_test".

  - `post_hoc_test`: A tibble with pairwise comparison results
    containing the following columns:

    - `variable`: Variable name

    - `ref_group`: Reference group

    - `test_group`: Test/treatment/case group

    - `p_val`: Adjusted p-value from emmeans pairwise test

    - `p_adj`: Adjusted p-value from emmeans pairwise test

- `raw_result`: A list containing:

  - `main_test`: A list of raw `aov` model objects.

  - `post_hoc_test`: A list of raw emmeans results.

- `meta_data`: Metadata copied from the input experiment.

- `preprocessing`: Input and analysis scales and the applied
  transformation; see
  [gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md).

## Details

For raw data with automatic transformation, the function performs log2
transformation on the expression data (log2(x + 1e-6)) before
statistical testing. At least 2 groups and at least 1 covariate are
required.

For any variable failed to fit a
[`stats::aov()`](https://rdrr.io/r/stats/aov.html) model, NAs will be
assigned to the results in both main test and post-hoc test.

**Post-hoc Test:** emmeans pairwise comparisons with Tukey adjustment
([`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html))
are performed for variables with significant main effects (p_adj \<
0.05).

## Required packages

This function requires the `emmeans` package for post-hoc comparisons.

## See also

[`stats::aov()`](https://rdrr.io/r/stats/aov.html),
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html),
[`emmeans::contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.html)

[gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md)
for input scales, transformations, and fold changes.
