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
  ...
)

gly_ancova_(expr_mat, groups, covariates, p_adj_method = "BH", ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- group_col:

  (Only for `gly_ancova()`) A character string specifying the column
  name of the grouping variable in the sample information. Default is
  `"group"`.

- covariate_cols:

  (Only for `gly_ancova()`) A character vector specifying column names
  in sample information to include as covariates. At least one covariate
  must be provided.

- p_adj_method:

  A character string specifying the method to adjust p-values. See
  `p.adjust.methods` for available methods. Default is "BH". If NULL, no
  adjustment is performed.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned. Only applicable to `gly_ancova()`.

- ...:

  Additional arguments passed to
  [`stats::aov()`](https://rdrr.io/r/stats/aov.html).

- expr_mat:

  (Only for `gly_ancova_()`) A numeric matrix with variables as rows and
  samples as columns.

- groups:

  (Only for `gly_ancova_()`) A factor or character vector specifying
  group membership for each sample. Must have at least 2 levels.
  Character vectors will be automatically converted to factors.

- covariates:

  (Only for `gly_ancova_()`) A data frame, matrix, or vector of
  sample-level covariates. At least one covariate must be provided.

## Value

A list containing two elements:

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

## Details

The function performs log2 transformation on the expression data
(log2(x + 1)) before statistical testing. At least 2 groups and at least
1 covariate are required.

`gly_ancova()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_ancova_()` is the underlying API that works with matrices, factor
vectors, and covariate data directly, providing more flexibility for
users who don't use the glyexp package.

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
