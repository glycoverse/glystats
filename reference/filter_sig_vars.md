# Filter Significant Variables

Filtering the experiment to keep only significant variables is a common
task. This function provides a convenient way to do this. It supports
results from all glystats DEA functions including
[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md),
[`gly_ancova()`](https://glycoverse.github.io/glystats/reference/gly_ancova.md),
[`gly_kruskal()`](https://glycoverse.github.io/glystats/reference/gly_kruskal.md),
[`gly_ttest()`](https://glycoverse.github.io/glystats/reference/gly_ttest.md),
[`gly_wilcox()`](https://glycoverse.github.io/glystats/reference/gly_wilcox.md),
[`gly_limma()`](https://glycoverse.github.io/glystats/reference/gly_limma.md),
and
[`gly_linear_model()`](https://glycoverse.github.io/glystats/reference/gly_linear_model.md).

## Usage

``` r
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  ...
)

# S3 method for class 'glystats_anova_res'
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  on = "main_test",
  comparison = NULL,
  ...
)

# S3 method for class 'glystats_ancova_res'
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  on = "main_test",
  comparison = NULL,
  ...
)

# S3 method for class 'glystats_kruskal_res'
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  on = "main_test",
  comparison = NULL,
  ...
)

# S3 method for class 'glystats_ttest_res'
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  ...
)

# S3 method for class 'glystats_wilcox_res'
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  ...
)

# S3 method for class 'glystats_limma_res'
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  comparison = NULL,
  ...
)

# S3 method for class 'glystats_linear_model_res'
filter_sig_vars(
  exp,
  res,
  p_adj_cutoff = 0.05,
  p_val_cutoff = NULL,
  fc_cutoff = NULL,
  term = NULL,
  ...
)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment`. Use the same object used to
  generate the DEA result.

- res:

  A glystats result object from a glystats DEA function.

- p_adj_cutoff:

  The threshold for p-adjusted values. Default is 0.05.

- p_val_cutoff:

  The threshold for p-values. We don't recommend using this. Default is
  NULL. If you insist to use it, please set `p_adj_cutoff` to NULL.

- fc_cutoff:

  The threshold for fold changes. Only positive value is needed. For
  example, `2` means fold change \> 2 or \< 1/2. Default is `2` for
  glycoproteomics data and `NULL` for others.

- ...:

  Additional arguments passed to methods. See the method-specific
  documentation for details.

- on:

  (For
  [`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/reference/gly_kruskal.md)
  results only) "main_test" or "post_hoc_test". Should the filter be
  applied on the main test results or the post-hoc test results? Default
  is "main_test". If "post_hoc_test", please set a `comparison` value.

- comparison:

  (For
  [`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md),
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/reference/gly_kruskal.md),
  and
  [`gly_limma()`](https://glycoverse.github.io/glystats/reference/gly_limma.md)
  results only) Specifies which comparison to filter on. A string with
  the format "group1_vs_group2". For
  [`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/reference/gly_kruskal.md)
  results, `comparison` is only used when `on` is "post_hoc_test". If
  not provided, filtering will be performed on the main test results for
  [`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/reference/gly_kruskal.md),
  and variables will be kept if they are significant in any comparison
  for
  [`gly_limma()`](https://glycoverse.github.io/glystats/reference/gly_limma.md).

- term:

  (For
  [`gly_linear_model()`](https://glycoverse.github.io/glystats/reference/gly_linear_model.md)
  results only) A coefficient or contrast name to filter on. If `NULL`,
  variables significant for any reported term are retained.

## Value

A filtered object of the same class as `exp`.

## Examples

``` r
library(glyexp)
library(glyclean)
#> 
#> Attaching package: ‘glyclean’
#> The following object is masked from ‘package:stats’:
#> 
#>     aggregate

exp <- auto_clean(real_experiment) |>
  glyexp::slice_head_row(n = 10)
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 24 (0.56%) variables.
#> ✔ Variable removal completed.
#> 
#> ── Normalizing data ──
#> 
#> ℹ Normalization method: `normalize_median()`
#> ℹ Reason: default for "glycoproteomics".
#> ✔ Normalization completed.
#> 
#> ── Imputing missing values ──
#> 
#> ℹ Imputation method: `impute_min_prob()`
#> ℹ Reason: default for "glycoproteomics" with n_samples < 30.
#> ✔ Imputation completed.
#> 
#> ── Aggregating data ──
#> 
#> ℹ Aggregating to "gfs" level
#> ✔ Aggregation completed.
#> 
#> ── Normalizing data again ──
#> 
#> ℹ Normalization method: `normalize_median()`
#> ℹ Reason: default for "glycoproteomics".
#> ✔ Normalization completed.
#> 
#> ── Correcting batch effects ──
#> 
#> ℹ Batch column batch not found in sample_info. Skipping batch correction.
#> ✔ Batch correction completed.
res <- gly_anova(exp)
#> ℹ Number of groups: 4
#> ℹ Groups: "H", "M", "Y", and "C"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
sig_exp <- filter_sig_vars(exp, res)
```
