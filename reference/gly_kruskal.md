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

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned.

- ...:

  Additional arguments passed to
  [`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html).

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

    - `effect_size`: Epsilon-squared

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

  - `post_hoc_test`: A list of raw `dunnTest` objects. Post-hoc
    comparison labels and direction-sensitive statistics follow the
    package direction, i.e. `test_group - ref_group`.

- `meta_data`: A list containing metadata from the input experiment.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1e-6)) before statistical testing. At least 2 groups are
required in the grouping variable.

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
