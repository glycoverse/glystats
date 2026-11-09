# Getting Started with glystats

When working with omics data, it is often useful to begin with a small
set of reliable analyses: differential expression, PCA, survival
analysis, and more. `glystats` brings these methods together through one
consistent, user-friendly interface.

**A useful feature:** Every `glystats` analysis accepts
[`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
and
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
containers, the unified data interfaces used across the `glycoverse`
ecosystem.

If these containers are new to you, the [glyexp
introduction](https://glycoverse.github.io/glyexp/articles/glyexp.html)
provides a useful overview.

``` r

library(glystats)
library(glyexp)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
library(glyread)
library(glyclean)
library(dplyr)
```

## Abundance scale and transformation

Record the supplied abundance scale in `S4Vectors::metadata(exp)$scale`:
`"raw"`, `"log"` for natural logarithms, `"log2"`, `"log10"`, or
`"logx"`. Setting metadata describes the data; it does not transform
them. For older development objects labeled `"log"`, set the actual base
explicitly if it was 2 or 10.

Every analysis accepts `log_transform = NULL` by default. Known logged
data are used unchanged. Raw data and old objects without scale metadata
retain each method’s original transformation and pseudo-count. ROC, Cox,
and fold change continue to use unlogged data by default. Use named
`log_transform = TRUE` to request a log transformation of raw values, or
`log_transform = FALSE` to use the supplied values unchanged. On an
already logged scale, explicit `TRUE` produces a warning and is ignored,
so the analysis continues without another logarithm. See
[`?gly_scale`](https://glycoverse.github.io/glystats/reference/gly_scale.md)
for the formulas used by each method.

For example, after supplying base-10 logged abundances:

``` r

S4Vectors::metadata(exp)$scale <- "log10"
result <- gly_ttest(exp)
result$preprocessing # input_scale = "log10", analysis_scale = "log10"
```

The input is unchanged. Results retain its metadata and record the
actual transformation separately in `preprocessing`. Estimates, standard
errors, intervals, and raw models use the analysis scale. `log2fc` from
ordinary logged data is the difference of group means converted to base
2 (a geometric mean ratio), whereas raw inputs retain the log2 ratio of
arithmetic means. Natural-log differences are divided by `log(2)`;
base-10 differences are multiplied by `log2(10)`. Limma on unlogged raw
data reports `estimate` and `NA` `log2fc`.

### Transformed coordinates: logx

The `"logx"` scale covers CLR, ALR, VSN, and other transformed
logarithmic or log-like coordinates that cannot be converted back to
original abundances from the label alone. It permits negative and zero
values and does not imply a known logarithm base. All analyses except
[`gly_fold_change()`](https://glycoverse.github.io/glystats/reference/gly_fold_change.md)
can use `logx` with any `log_transform` setting, preserving the supplied
coordinates. Explicit `TRUE` produces a warning and is ignored. Check
the analysis assumptions for the particular transformation, especially
for compositional data.

For example, a CLR experiment from `glyclean` can be analyzed directly:

``` r

clr_exp <- glyclean::transform_clr(exp)
pca <- gly_pca(clr_exp)
de <- gly_ttest(clr_exp)
de$preprocessing # input_scale = "logx", analysis_scale = "logx"
de$tidy_result   # native estimates and effect sizes; log2fc is NA
sig <- filter_sig_vars(clr_exp, de) # significance only
```

[`gly_fold_change()`](https://glycoverse.github.io/glystats/reference/gly_fold_change.md)
rejects `logx` rather than reporting coordinate differences as abundance
fold changes. Differential tests retain their native estimates and
effect sizes but report `NA` for `log2fc`, including post-hoc
comparisons. Limma reports its native contrast in `estimate`.
[`filter_sig_vars()`](https://glycoverse.github.io/glystats/reference/filter_sig_vars.md)
defaults to significance-only filtering for `logx`, including for
glycoproteomics; a non-`NULL` `fc_cutoff` warns and is ignored. Explicit
`fc_cutoff = NULL` disables fold-change filtering on every scale.
Omitting the argument retains the default cutoff of 2 for
glycoproteomics comparisons with available fold changes. Supply an
experiment on a known abundance scale to compute fold changes; do not
relabel transformed coordinates as ordinary log values.

## Explore the data

We will start by exploring the demo dataset:

``` r

# Here we use `glyread` to read in the data and use `glyclean` to perform preprocessing
exp <- read_pglyco3_pglycoquant("glycopeptides.list", sample_info = "sample_info.csv") |>
  auto_clean()
```

``` r

exp
#> 
#> ── GlycoproteomicSE ────────────────────────────────────────────────────────────
#> ℹ Abundance assay: 12 samples, 225 variables
#> ℹ Glycan type: N
#> ℹ Row data fields: protein <chr>, glycan_composition <comp>, protein_site <int>, gene <chr>, protein_group <chr>
#> ℹ Column data fields: group <fct>
#> ℹ Metadata fields: exp_type <chr>, glycan_type <chr>, quant_method <chr>, scale <chr>, pipeline <glycln_p>
```

The result is a
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
containing 12 samples and 263 glycoforms, which is enough to illustrate
the main workflow.

``` r

tibble::as_tibble(rowData(exp), rownames = "variable")
#> # A tibble: 225 × 6
#>    variable          protein glycan_composition protein_site gene  protein_group
#>    <chr>             <chr>   <comp>                    <int> <chr> <chr>        
#>  1 P08185-176-Hex(5… P08185  Hex(5)HexNAc(4)Ne…          176 SERP… P08185       
#>  2 P04196-344-Hex(5… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   P04196       
#>  3 P04196-344-Hex(5… P04196  Hex(5)HexNAc(4)             344 HRG   P04196       
#>  4 P10909-291-Hex(6… P10909  Hex(6)HexNAc(5)             291 CLU   P10909       
#>  5 P04196-344-Hex(5… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   P04196       
#>  6 P04196-345-Hex(5… P04196  Hex(5)HexNAc(4)             345 HRG   P04196       
#>  7 P04196-344-Hex(5… P04196  Hex(5)HexNAc(4)dH…          344 HRG   P04196       
#>  8 P04196-344-Hex(4… P04196  Hex(4)HexNAc(3)             344 HRG   P04196       
#>  9 P04196-344-Hex(4… P04196  Hex(4)HexNAc(4)Ne…          344 HRG   P04196       
#> 10 P10909-291-Hex(5… P10909  Hex(5)HexNAc(4)             291 CLU   P10909       
#> # ℹ 215 more rows
```

The variable information tibble contains descriptive information for
each glycoform, including the protein, glycosylation site, and glycan
structures.

``` r

tibble::as_tibble(colData(exp), rownames = "sample")
#> # A tibble: 12 × 2
#>    sample                  group
#>    <chr>                   <fct>
#>  1 20241224-LXJ-Nglyco-H_1 H    
#>  2 20241224-LXJ-Nglyco-H_2 H    
#>  3 20241224-LXJ-Nglyco-H_3 H    
#>  4 20241224-LXJ-Nglyco-M_1 M    
#>  5 20241224-LXJ-Nglyco-M_2 M    
#>  6 20241224-LXJ-Nglyco-M_3 M    
#>  7 20241224-LXJ-Nglyco-Y_1 Y    
#>  8 20241224-LXJ-Nglyco-Y_2 Y    
#>  9 20241224-LXJ-Nglyco-Y_3 Y    
#> 10 20241224-LXJ-Nglyco-C_1 C    
#> 11 20241224-LXJ-Nglyco-C_2 C    
#> 12 20241224-LXJ-Nglyco-C_3 C
```

The sample information tibble includes a `group` column. In this
dataset, “H” represents healthy samples, “M” hepatitis, “Y” cirrhosis,
and “C” hepatocellular carcinoma. We will use this column for the
comparisons below.

## Use a consistent interface

`glystats` uses a consistent naming pattern: `gly_` followed by the
analysis name. For example,
[`gly_ttest()`](https://glycoverse.github.io/glystats/reference/gly_ttest.md)
runs a t-test and
[`gly_pca()`](https://glycoverse.github.io/glystats/reference/gly_pca.md)
performs PCA.

RStudio’s auto-completion can help you discover the available functions.

We can use ANOVA to identify glycoforms that differ between groups:

``` r

anova_res <- gly_anova(exp)
#> ℹ Number of groups: 4
#> ℹ Groups: "C", "H", "M", and "Y"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
```

[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
follows the `glycoverse` column conventions, detects the `group` column
in the sample information, and fits an ANOVA model for each glycoform.

## Work with the results

All `glystats` functions return a consistent list with two main
components:

- `tidy_result`: Analysis-ready tibbles in tidy format
- `raw_result`: The original result objects from the underlying
  statistical functions

For
[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md),
the `tidy_result` contains two informative tibbles: `main_test` and
`post_hoc_test`.

You can use
[`get_tidy_result()`](https://glycoverse.github.io/glystats/reference/get_tidy_result.md)
to get the tidy result tibble:

``` r

get_tidy_result(anova_res, "main_test")
#> # A tibble: 225 × 15
#>    .variable   protein glycan_composition protein_site gene  protein_group term 
#>    <chr>       <chr>   <comp>                    <int> <chr> <chr>         <chr>
#>  1 P08185-176… P08185  Hex(5)HexNAc(4)Ne…          176 SERP… P08185        group
#>  2 P04196-344… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   P04196        group
#>  3 P04196-344… P04196  Hex(5)HexNAc(4)             344 HRG   P04196        group
#>  4 P10909-291… P10909  Hex(6)HexNAc(5)             291 CLU   P10909        group
#>  5 P04196-344… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   P04196        group
#>  6 P04196-345… P04196  Hex(5)HexNAc(4)             345 HRG   P04196        group
#>  7 P04196-344… P04196  Hex(5)HexNAc(4)dH…          344 HRG   P04196        group
#>  8 P04196-344… P04196  Hex(4)HexNAc(3)             344 HRG   P04196        group
#>  9 P04196-344… P04196  Hex(4)HexNAc(4)Ne…          344 HRG   P04196        group
#> 10 P10909-291… P10909  Hex(5)HexNAc(4)             291 CLU   P10909        group
#> # ℹ 215 more rows
#> # ℹ 8 more variables: df <dbl>, sumsq <dbl>, meansq <dbl>, statistic <dbl>,
#> #   p_val <dbl>, p_adj <dbl>, effect_size <dbl>, post_hoc <chr>
```

[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
also adds the descriptive columns from the variable tibble to the tidy
results. Use the `add_info` parameter when you want to control this
behavior.

The `raw_result` houses two lists of models - one for the main test, one
for post hoc comparisons:

``` r

names(get_raw_result(anova_res))
#> [1] "main_test"     "post_hoc_test"
```

[`get_tidy_result()`](https://glycoverse.github.io/glystats/reference/get_tidy_result.md)
and
[`get_raw_result()`](https://glycoverse.github.io/glystats/reference/get_tidy_result.md)
can also be used in a pipe:

``` r

exp |>
  gly_anova() |>
  get_tidy_result("main_test") |>
  filter(p_adj < 0.05)
#> ℹ Number of groups: 4
#> ℹ Groups: "C", "H", "M", and "Y"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
#> # A tibble: 58 × 15
#>    .variable   protein glycan_composition protein_site gene  protein_group term 
#>    <chr>       <chr>   <comp>                    <int> <chr> <chr>         <chr>
#>  1 P04196-344… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   P04196        group
#>  2 P04196-344… P04196  Hex(5)HexNAc(4)             344 HRG   P04196        group
#>  3 P10909-291… P10909  Hex(6)HexNAc(5)             291 CLU   P10909        group
#>  4 P04196-344… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   P04196        group
#>  5 P04196-345… P04196  Hex(5)HexNAc(4)             345 HRG   P04196        group
#>  6 P04196-344… P04196  Hex(5)HexNAc(4)dH…          344 HRG   P04196        group
#>  7 P04196-344… P04196  Hex(4)HexNAc(3)             344 HRG   P04196        group
#>  8 P10909-291… P10909  Hex(5)HexNAc(4)             291 CLU   P10909        group
#>  9 P04196-344… P04196  Hex(5)HexNAc(4)dH…          344 HRG   P04196        group
#> 10 P13671-855… P13671  Hex(5)HexNAc(4)dH…          855 C6    P13671        group
#> # ℹ 48 more rows
#> # ℹ 8 more variables: df <dbl>, sumsq <dbl>, meansq <dbl>, statistic <dbl>,
#> #   p_val <dbl>, p_adj <dbl>, effect_size <dbl>, post_hoc <chr>
```

## Available analyses

The package includes the following analyses for glycomics and
glycoproteomics data:

- **Differential expression analysis:**
  - [`gly_ttest()`](https://glycoverse.github.io/glystats/reference/gly_ttest.md):
    Two-sample t-test
  - [`gly_wilcox()`](https://glycoverse.github.io/glystats/reference/gly_wilcox.md):
    Wilcoxon rank sum test
  - [`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md):
    One-way ANOVA
  - [`gly_kruskal()`](https://glycoverse.github.io/glystats/reference/gly_kruskal.md):
    Kruskal-Wallis rank sum test
  - [`gly_limma()`](https://glycoverse.github.io/glystats/reference/gly_limma.md):
    Linear models for microarray data (limma)
  - [`gly_ancova()`](https://glycoverse.github.io/glystats/reference/gly_ancova.md):
    ANCOVA (Analysis of Covariance)
  - [`gly_linear_model()`](https://glycoverse.github.io/glystats/reference/gly_linear_model.md):
    Formula-based moderated linear models
  - [`gly_set_test()`](https://glycoverse.github.io/glystats/reference/gly_set_test.md):
    correlated variable-set testing
  - [`gly_fold_change()`](https://glycoverse.github.io/glystats/reference/gly_fold_change.md):
    Calculate fold change
- **Dimensionality reduction:**
  - [`gly_pca()`](https://glycoverse.github.io/glystats/reference/gly_pca.md):
    Principal component analysis
  - [`gly_tsne()`](https://glycoverse.github.io/glystats/reference/gly_tsne.md):
    t-distributed stochastic neighbor embedding (t-SNE)
  - [`gly_umap()`](https://glycoverse.github.io/glystats/reference/gly_umap.md):
    Uniform manifold approximation and projection (UMAP)
  - [`gly_oplsda()`](https://glycoverse.github.io/glystats/reference/gly_oplsda.md):
    Orthogonal partial least squares discriminant analysis (OPLS-DA)
  - [`gly_plsda()`](https://glycoverse.github.io/glystats/reference/gly_plsda.md):
    Partial least squares discriminant analysis (PLS-DA)
- **Clustering:**
  - [`gly_kmeans()`](https://glycoverse.github.io/glystats/reference/gly_kmeans.md):
    K-means clustering
  - [`gly_hclust()`](https://glycoverse.github.io/glystats/reference/gly_hclust.md):
    Hierarchical clustering
- **Survival analysis:**
  - [`gly_cox()`](https://glycoverse.github.io/glystats/reference/gly_cox.md):
    Cox proportional hazards model
- **Additional tools:**
  - [`gly_cor()`](https://glycoverse.github.io/glystats/reference/gly_cor.md):
    Correlation analysis
  - [`gly_roc()`](https://glycoverse.github.io/glystats/reference/gly_roc.md):
    Receiver operating characteristic (ROC) analysis
