# One-way or Repeated-Measures ANOVA for Differential Expression Analysis

Perform one-way or repeated-measures ANOVA for glycomics or
glycoproteomics data. For significant results, Tukey's HSD post-hoc test
is automatically performed. P-values are adjusted for multiple testing
using the method specified by `p_adj_method`.

## Usage

``` r
gly_anova(
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
  [`stats::aov()`](https://rdrr.io/r/stats/aov.html).

- subject_col:

  An optional character string naming the subject identifier column in
  sample information. When supplied, a subject-blocked repeated-measures
  ANOVA is performed.

## Value

A list containing three elements:

- `tidy_result`: A list containing:

  - `main_test`: A tibble with ANOVA results containing the following
    columns:

    - `variable`: Variable name

    - `term`: ANOVA term (usually "groups")

    - `df`: Degrees of freedom

    - `sumsq`: Sum of squares

    - `meansq`: Mean squares

    - `statistic`: F-statistic

    - `p_val`: Raw p-value from ANOVA

    - `effect_size`: Eta-squared, or partial eta-squared for paired
      analyses

    - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

    - `post_hoc`: Significant group pairs from post-hoc test, in the
      format of "ref_vs_test".

  - `post_hoc_test`: A tibble with pairwise comparison results
    containing the following columns:

    - `variable`: Variable name

    - `ref_group`: Reference group

    - `test_group`: Test/treatment/case group

    - `p_val`: Raw post-hoc p-value

    - `p_adj`: Adjusted post-hoc p-value

- `raw_result`: A list containing:

  - `main_test`: A list of raw `aov` model objects.

  - `post_hoc_test`: A list of raw `TukeyHSD` objects or paired t-test
    summaries.

- `meta_data`: A list containing metadata from the input experiment.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1e-6)) before statistical testing. At least 2 groups are
required in the grouping variable. In paired analyses, only subjects
observed in every group are used, with at most one sample per subject
and group. The reported effect size is partial eta-squared for the group
term.

For any variable failed to fit a
[`stats::aov()`](https://rdrr.io/r/stats/aov.html) model, NAs will be
assigned to the results in both main test and post-hoc test.

**Post-hoc Test:** Tukey's HSD test for pairwise comparisons
([`stats::TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html)) is
performed for variables with significant main effects (p_adj \< 0.05).
Paired analyses instead use paired t-tests with Holm correction.

## See also

[`stats::aov()`](https://rdrr.io/r/stats/aov.html),
[`stats::TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html)
