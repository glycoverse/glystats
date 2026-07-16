# Get the result of a glystats analysis

Syntax sugar for `$tidy_result` and `$raw_result` elements of a glystats
result object. It's useful to be used in pipes.

## Usage

``` r
get_tidy_result(res, which = NULL)

get_raw_result(res, which = NULL)
```

## Arguments

- res:

  A glystats result object.

- which:

  Used to specify which element to get, when the result is a list. For
  glystats results with only one tidy result (e.g.
  [`gly_ttest()`](https://glycoverse.github.io/glystats/reference/gly_ttest.md)),
  this argument is not needed. For others (e.g.
  [`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)),
  this argument is required to specify which tidy result to get.

## Value

A tibble.

## Examples

``` r
library(glyexp)
library(glyclean)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

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

# Using a pipe
sig_res <- exp |>
  gly_anova() |>
  get_tidy_result("main_test") |>
  filter(p_adj < 0.05)
#> ℹ Number of groups: 4
#> ℹ Groups: "H", "M", "Y", and "C"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.

# Equivalent to
anova_res <- gly_anova(exp)
#> ℹ Number of groups: 4
#> ℹ Groups: "H", "M", "Y", and "C"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
sig_res <- anova_res$tidy_result$main_test |>
  filter(p_adj < 0.05)
```
