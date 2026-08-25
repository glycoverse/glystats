# Kruskal-Wallis or Friedman Test for Differential Expression Analysis

Perform Kruskal-Wallis or Friedman tests for glycomics or
glycoproteomics data. For significant results, Dunn's post-hoc test is
automatically performed. P-values are adjusted for multiple testing
using the method specified by `p_adj_method`.

## Usage

``` r
gly_kruskal(
  exp,
  group_col = "group",
  p_adj_method = "BH",
  add_info = TRUE,
  ...,
  subject_col = NULL
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

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned.

- ...:

  Additional arguments passed to
  [`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)
  or, for paired analyses,
  [`stats::friedman.test()`](https://rdrr.io/r/stats/friedman.test.html).

- subject_col:

  An optional character string naming the subject identifier column in
  sample information. When supplied, a Friedman test is performed.

## Value

A list containing three elements:

- `tidy_result`: A list containing:

  - `main_test`: A tibble with Kruskal-Wallis or Friedman test results
    containing the following columns:

    - `variable`: Variable name

    - `statistic`: Kruskal-Wallis test statistic

    - `p_val`: Raw p-value from Kruskal-Wallis test

    - `parameter`: Degrees of freedom

    - `method`: Statistical method used

    - `effect_size`: Epsilon-squared, or Kendall's W for paired analyses

    - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

    - `post_hoc`: Significant group pairs from post-hoc test, in the
      format of "ref_vs_test".

  - `post_hoc_test`: A tibble with pairwise comparison results
    containing the following columns:

    - `variable`: Variable name

    - `ref_group`: Reference group

    - `test_group`: Test/treatment/case group

    - `p_val`: Raw post-hoc p-value

    - `p_adj`: Holm-adjusted post-hoc p-value

- `raw_result`: A list containing:

  - `main_test`: A list of raw `kruskal.test` or `friedman.test`
    objects.

  - `post_hoc_test`: A list of raw Dunn or paired Wilcoxon summaries.

- `meta_data`: A list containing metadata from the input experiment.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1e-6)) before statistical testing. At least 2 groups are
required in the grouping variable. In paired analyses, only subjects
observed in every group are used, with at most one sample per subject
and group. The reported effect size is Kendall's \\W\\.

For any variable failed to fit a
[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html) or
[`stats::friedman.test()`](https://rdrr.io/r/stats/friedman.test.html)
model, NAs will be assigned to the results in both main test and
post-hoc test.

**Post-hoc Test:** Dunn's test with Holm correction for multiple
comparisons
([`FSA::dunnTest()`](https://fishr-core-team.github.io/FSA/reference/dunnTest.html))
is performed for variables with significant main effects (p_adj \<
0.05). Paired analyses instead use paired Wilcoxon signed-rank tests
with Holm correction.

## Required packages

Unpaired analyses require the `FSA` package for Dunn's post-hoc test.

## See also

[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html),
[`stats::friedman.test()`](https://rdrr.io/r/stats/friedman.test.html),
[`FSA::dunnTest()`](https://fishr-core-team.github.io/FSA/reference/dunnTest.html)
