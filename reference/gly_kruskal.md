# Kruskal-Wallis test for Differential Expression Analysis

Perform Kruskal-Wallis test for glycomics or glycoproteomics data. The
function supports non-parametric comparison of multiple groups. For
significant results, Dunn's post-hoc test is automatically performed.
P-values are adjusted for multiple testing using the method specified by
`p_adj_method`.

## Usage

``` r
gly_kruskal(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  add_info = TRUE,
  ...
)

gly_kruskal_(expr_mat, groups, p_adj_method = "BH", ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- group_col:

  (Only for `gly_kruskal()`) A character string specifying the column
  name of the grouping variable in the sample information. Default is
  `"group"`.

- p_adj_method:

  A character string specifying the method to adjust p-values. See
  `p.adjust.methods` for available methods. Default is "BH". If NULL, no
  adjustment is performed.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned. Only applicable to `gly_kruskal()`.

- ...:

  Additional arguments passed to
  [`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html).

- expr_mat:

  (Only for `gly_kruskal_()`) A numeric matrix with variables as rows
  and samples as columns.

- groups:

  (Only for `gly_kruskal_()`) A factor or character vector specifying
  group membership for each sample. Must have at least 2 levels.
  Character vectors will be automatically converted to factors.

## Value

A list containing three elements:

- `tidy_result`: A list containing:

  - `main_test`: A tibble with Kruskal-Wallis test results containing
    the following columns:

    - `variable`: Variable name

    - `statistic`: Kruskal-Wallis test statistic

    - `p_val`: Raw p-value from Kruskal-Wallis test

    - `parameter`: Degrees of freedom

    - `method`: Statistical method used

    - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

    - `post_hoc`: Significant group pairs from post-hoc test, in the
      format of "ref_vs_test".

  - `post_hoc_test`: A tibble with pairwise comparison results
    containing the following columns:

    - `variable`: Variable name

    - `ref_group`: Reference group

    - `test_group`: Test/treatment/case group

    - `p_val`: Raw p-value from Dunn's test

    - `p_adj`: Adjusted p-value from Dunn's test

- `raw_result`: A list containing:

  - `main_test`: A list of raw `kruskal.test` objects.

  - `post_hoc_test`: A list of raw `dunnTest` objects.

- `meta_data`: A list containing metadata from the input experiment.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1)) before statistical testing. At least 2 groups are required
in the grouping variable.

`gly_kruskal()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_kruskal_()` is the underlying API that works with matrices and
factor vectors directly, providing more flexibility for users who don't
use the glyexp package.

For any variable failed to fit a
[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)
model, NAs will be assigned to the results in both main test and
post-hoc test.

**Post-hoc Test:** Dunn's test with Holm correction for multiple
comparisons
([`FSA::dunnTest()`](https://fishr-core-team.github.io/FSA/reference/dunnTest.html))
is performed for variables with significant main effects (p_adj \<
0.05).

## Required packages

This function requires the `FSA` package for Dunn's post-hoc test.

## See also

[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html),
[`FSA::dunnTest()`](https://fishr-core-team.github.io/FSA/reference/dunnTest.html)
