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
#> ℹ Row data fields: protein <chr>, glycan_composition <comp>, protein_site <int>, gene <chr>
#> ℹ Column data fields: group <fct>
#> ℹ Metadata fields: exp_type <chr>, glycan_type <chr>, quant_method <chr>
```

The result is a
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
containing 12 samples and 263 glycoforms, which is enough to illustrate
the main workflow.

``` r

tibble::as_tibble(rowData(exp), rownames = "variable")
#> # A tibble: 225 × 5
#>    variable                        protein glycan_composition protein_site gene 
#>    <chr>                           <chr>   <comp>                    <int> <chr>
#>  1 P08185-176-Hex(5)HexNAc(4)NeuA… P08185  Hex(5)HexNAc(4)Ne…          176 SERP…
#>  2 P04196-344-Hex(5)HexNAc(4)NeuA… P04196  Hex(5)HexNAc(4)Ne…          344 HRG  
#>  3 P04196-344-Hex(5)HexNAc(4)      P04196  Hex(5)HexNAc(4)             344 HRG  
#>  4 P10909-291-Hex(6)HexNAc(5)      P10909  Hex(6)HexNAc(5)             291 CLU  
#>  5 P04196-344-Hex(5)HexNAc(4)NeuA… P04196  Hex(5)HexNAc(4)Ne…          344 HRG  
#>  6 P04196-345-Hex(5)HexNAc(4)      P04196  Hex(5)HexNAc(4)             345 HRG  
#>  7 P04196-344-Hex(5)HexNAc(4)dHex… P04196  Hex(5)HexNAc(4)dH…          344 HRG  
#>  8 P04196-344-Hex(4)HexNAc(3)      P04196  Hex(4)HexNAc(3)             344 HRG  
#>  9 P04196-344-Hex(4)HexNAc(4)NeuA… P04196  Hex(4)HexNAc(4)Ne…          344 HRG  
#> 10 P10909-291-Hex(5)HexNAc(4)      P10909  Hex(5)HexNAc(4)             291 CLU  
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
#> # A tibble: 225 × 14
#>    variable     protein glycan_composition protein_site gene  term     df  sumsq
#>    <chr>        <chr>   <comp>                    <int> <chr> <chr> <dbl>  <dbl>
#>  1 P08185-176-… P08185  Hex(5)HexNAc(4)Ne…          176 SERP… group     3  55.7 
#>  2 P04196-344-… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 158.  
#>  3 P04196-344-… P04196  Hex(5)HexNAc(4)             344 HRG   group     3 138.  
#>  4 P10909-291-… P10909  Hex(6)HexNAc(5)             291 CLU   group     3  22.5 
#>  5 P04196-344-… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 496.  
#>  6 P04196-345-… P04196  Hex(5)HexNAc(4)             345 HRG   group     3  60.2 
#>  7 P04196-344-… P04196  Hex(5)HexNAc(4)dH…          344 HRG   group     3 173.  
#>  8 P04196-344-… P04196  Hex(4)HexNAc(3)             344 HRG   group     3  72.1 
#>  9 P04196-344-… P04196  Hex(4)HexNAc(4)Ne…          344 HRG   group     3   9.65
#> 10 P10909-291-… P10909  Hex(5)HexNAc(4)             291 CLU   group     3  75.0 
#> # ℹ 215 more rows
#> # ℹ 6 more variables: meansq <dbl>, statistic <dbl>, p_val <dbl>, p_adj <dbl>,
#> #   effect_size <dbl>, post_hoc <chr>
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
#> # A tibble: 54 × 14
#>    variable      protein glycan_composition protein_site gene  term     df sumsq
#>    <chr>         <chr>   <comp>                    <int> <chr> <chr> <dbl> <dbl>
#>  1 P04196-344-H… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 158. 
#>  2 P04196-344-H… P04196  Hex(5)HexNAc(4)             344 HRG   group     3 138. 
#>  3 P04196-344-H… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 496. 
#>  4 P10909-291-H… P10909  Hex(5)HexNAc(4)             291 CLU   group     3  75.0
#>  5 P04196-344-H… P04196  Hex(5)HexNAc(4)dH…          344 HRG   group     3 127. 
#>  6 P13671-855-H… P13671  Hex(5)HexNAc(4)dH…          855 C6    group     3  81.1
#>  7 P04196-344-H… P04196  Hex(4)HexNAc(3)dH…          344 HRG   group     3 195. 
#>  8 P04196-344-H… P04196  Hex(5)HexNAc(4)dH…          344 HRG   group     3 116. 
#>  9 P01019-161-H… P01019  Hex(5)HexNAc(4)Ne…          161 AGT   group     3  46.4
#> 10 P01019-161-H… P01019  Hex(4)HexNAc(3)Ne…          161 AGT   group     3  61.4
#> # ℹ 44 more rows
#> # ℹ 6 more variables: meansq <dbl>, statistic <dbl>, p_val <dbl>, p_adj <dbl>,
#> #   effect_size <dbl>, post_hoc <chr>
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
