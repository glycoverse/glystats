# One-way ANOVA for Differential Expression Analysis

Perform one-way ANOVA for glycomics or glycoproteomics data. The
function supports parametric comparison of multiple groups. For
significant results, Tukey's HSD post-hoc test is automatically
performed. P-values are adjusted for multiple testing using the method
specified by `p_adj_method`.

## Usage

``` r
gly_anova(exp, group_col = "group", p_adj_method = "BH", add_info = TRUE, ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

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

    - `effect_size`: Eta-squared

    - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

    - `post_hoc`: Significant group pairs from post-hoc test, in the
      format of "ref_vs_test".

  - `post_hoc_test`: A tibble with pairwise comparison results
    containing the following columns:

    - `variable`: Variable name

    - `ref_group`: Reference group

    - `test_group`: Test/treatment/case group

    - `p_val`: Raw p-value from Tukey's HSD test

    - `p_adj`: Adjusted p-value from Tukey's HSD test

- `raw_result`: A list containing:

  - `main_test`: A list of raw `aov` model objects.

  - `post_hoc_test`: A list of raw `TukeyHSD` objects. Post-hoc
    comparison labels follow the package direction, i.e.
    `test_group - ref_group`.

- `meta_data`: A list containing metadata from the input experiment.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1e-6)) before statistical testing. At least 2 groups are
required in the grouping variable.

For any variable failed to fit a
[`stats::aov()`](https://rdrr.io/r/stats/aov.html) model, NAs will be
assigned to the results in both main test and post-hoc test.

**Post-hoc Test:** Tukey's HSD test for pairwise comparisons
([`stats::TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html)) is
performed for variables with significant main effects (p_adj \< 0.05).

## See also

[`stats::aov()`](https://rdrr.io/r/stats/aov.html),
[`stats::TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html)
