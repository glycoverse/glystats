# Filter Significant Variables

Filtering the experiment to keep only significant variables is a common
task. This function provides a convenient way to do this. It supports
results from all glystats DEA functions including
[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md),
[`gly_kruskal()`](https://glycoverse.github.io/glystats/reference/gly_kruskal.md),
[`gly_ttest()`](https://glycoverse.github.io/glystats/reference/gly_ttest.md),
[`gly_wilcox()`](https://glycoverse.github.io/glystats/reference/gly_wilcox.md),
and
[`gly_limma()`](https://glycoverse.github.io/glystats/reference/gly_limma.md).

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
```

## Arguments

- exp:

  An
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html).

- res:

  A glystats result object from a glystats DEA function.

- p_adj_cutoff:

  The threshold for p-adjusted values. Default is 0.05.

- p_val_cutoff:

  The threshold for p-values. We don't recommend using this. Default is
  NULL. If you insist to use it, please set `p_adj_cutoff` to NULL.

- fc_cutoff:

  The threshold for fold changes. Default is NULL. Only positive value
  is needed. For example, `2` means fold change \> 2 or \< 1/2.

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

## Value

An new
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

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
  glyexp::slice_head_var(n = 10)
#> 
#> ── Normalizing data ──
#> 
#> ℹ No QC samples found. Using default normalization method based on experiment type.
#> ℹ Experiment type is "glycoproteomics". Using `normalize_median()`.
#> ✔ Normalization completed.
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ No QC samples found. Using all samples.
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 24 (0.56%) variables.
#> ✔ Variable removal completed.
#> 
#> ── Imputing missing values ──
#> 
#> ℹ No QC samples found. Using default imputation method based on sample size.
#> ℹ Sample size <= 30, using `impute_sample_min()`.
#> ✔ Imputation completed.
#> 
#> ── Aggregating data ──
#> 
#> ℹ Aggregating to "gfs" level
#> ✔ Aggregation completed.
#> 
#> ── Normalizing data again ──
#> 
#> ℹ No QC samples found. Using default normalization method based on experiment type.
#> ℹ Experiment type is "glycoproteomics". Using `normalize_median()`.
#> ✔ Normalization completed.
#> 
#> ── Correcting batch effects ──
#> 
#> ℹ Batch column  not found in sample_info. Skipping batch correction.
#> ✔ Batch correction completed.
res <- gly_anova(exp)
#> ℹ Number of groups: 4
#> ℹ Groups: "H", "M", "Y", and "C"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
sig_exp <- filter_sig_vars(exp, res)
```
