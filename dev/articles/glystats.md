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
[`?gly_scale`](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md)
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
[`gly_fold_change()`](https://glycoverse.github.io/glystats/dev/reference/gly_fold_change.md)
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

[`gly_fold_change()`](https://glycoverse.github.io/glystats/dev/reference/gly_fold_change.md)
rejects `logx` rather than reporting coordinate differences as abundance
fold changes. Differential tests retain their native estimates and
effect sizes but report `NA` for `log2fc`, including post-hoc
comparisons. Limma reports its native contrast in `estimate`.
[`filter_sig_vars()`](https://glycoverse.github.io/glystats/dev/reference/filter_sig_vars.md)
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
[`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md)
runs a t-test and
[`gly_pca()`](https://glycoverse.github.io/glystats/dev/reference/gly_pca.md)
performs PCA.

RStudio’s auto-completion can help you discover the available functions.

We can use ANOVA to identify glycoforms that differ between groups:

``` r

anova_res <- gly_anova(exp)
#> ℹ Number of groups: 4
#> ℹ Groups: "C", "H", "M", and "Y"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
```

[`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
follows the `glycoverse` column conventions, detects the `group` column
in the sample information, and fits an ANOVA model for each glycoform.

## Work with the results

All `glystats` functions return a consistent list with two main
components:

- `tidy_result`: Analysis-ready tibbles in tidy format
- `raw_result`: The original result objects from the underlying
  statistical functions

For
[`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md),
the `tidy_result` contains two informative tibbles: `main_test` and
`post_hoc_test`.

You can use
[`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
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

[`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
also adds the descriptive columns from the variable tibble to the tidy
results. Use the `add_info` parameter when you want to control this
behavior.

The `raw_result` houses two lists of models - one for the main test, one
for post hoc comparisons:

``` r

names(get_raw_result(anova_res))
#> [1] "main_test"     "post_hoc_test"
```

[`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
and
[`get_raw_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
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

## Compare two independent groups with SAM

With the Bioconductor package `siggenes` installed, we can compare group
C with healthy samples using Significance Analysis of Microarrays (SAM).
The sign of `d` gives the direction on the analysis scale. Positive
values indicate higher mean abundance in C. For relative abundances, the
comparison is relative to the normalization reference.

``` r

two_groups <- exp[, exp$group %in% c("H", "C")]
sam <- gly_sam(two_groups, ref_group = "H", delta = 1, seed = 42, add_info = FALSE)
#> 
#> We're doing 20 complete permutations
sam$thresholds
#> # A tibble: 1 × 7
#>   delta cutoff_lower cutoff_upper n_called expected_false    p0    fdr
#>   <dbl>        <dbl>        <dbl>    <int>          <dbl> <dbl>  <dbl>
#> 1     1         -Inf         1.89        2            0.1 0.644 0.0322
get_tidy_result(sam) |>
  filter(called) |>
  select(.variable, d, p_val, q_value) |>
  head()
#> # A tibble: 2 × 4
#>   .variable                              d    p_val q_value
#>   <chr>                              <dbl>    <dbl>   <dbl>
#> 1 P04196-344-Hex(5)HexNAc(4)NeuAc(2)  2.27 0.000444  0.0644
#> 2 P00734-121-Hex(5)HexNAc(4)          1.89 0.000889  0.0644
```

Here `delta = 1` calls features whose scores reach the corresponding
lower or upper cutoff. The threshold table gives the number called and
their estimated permutation FDR. Each feature also has a pooled
permutation p-value and a separate q-value. Omit `delta` to inspect the
backend’s grid of thresholds before choosing one. Only complete features
with estimable within-group variation contribute to inference. The
`exclusion` column identifies the others, and `sam$diagnostics` reports
the tested family size and how many distinct label permutations were
used.

## Estimate posterior probabilities with EBAM

Empirical Bayes analysis of microarrays (EBAM) estimates the probability
that each glycoform differs between the two groups. With H as the
reference, a positive `z` indicates a higher mean in C on the analysis
scale.

``` r

ebam <- gly_ebam(two_groups, ref_group = "H", posterior_threshold = 0.9, seed = 2)
#> 
#> We're doing 20 complete permutations
ebam$thresholds
#> # A tibble: 1 × 6
#>   posterior_threshold cutoff_lower cutoff_upper n_called    p0    fdr
#>                 <dbl>        <dbl>        <dbl>    <int> <dbl>  <dbl>
#> 1                 0.9         -Inf         2.12       15 0.622 0.0518
get_tidy_result(ebam) |>
  filter(called) |>
  select(.variable, z, posterior, local_fdr) |>
  head()
#> # A tibble: 6 × 4
#>   .variable                                     z posterior local_fdr
#>   <chr>                                     <dbl>     <dbl>     <dbl>
#> 1 P04196-344-Hex(5)HexNAc(4)NeuAc(1)         4.32     0.965    0.0347
#> 2 P04196-344-Hex(5)HexNAc(4)                 2.12     0.906    0.0936
#> 3 P04196-344-Hex(5)HexNAc(4)NeuAc(2)         5.66     0.949    0.0514
#> 4 P04196-344-Hex(5)HexNAc(4)dHex(2)          3.46     0.964    0.0360
#> 5 P04196-344-Hex(4)HexNAc(3)                 3.22     0.961    0.0389
#> 6 P04196-344-Hex(5)HexNAc(4)dHex(1)NeuAc(2)  3.59     0.965    0.0350
```

`posterior` estimates the probability of differential abundance for one
feature. Its complement, `local_fdr`, estimates the probability of the
null for that feature. At the chosen threshold, siggenes identifies the
lower and upper score tails. Use `called` to recover the set of features
whose permutation FDR appears in `ebam$thresholds`. Supplying a vector
such as `c(0.8, 0.9, 0.95)` returns a row for each threshold and leaves
feature calls unspecified.

EBAM fits complete features with estimable within-group variation.
Excluded features keep their IDs and reasons in the tidy table. With
three samples per group, requesting 100 permutations uses all 20
possible label allocations. `ebam$settings$n_perm` and
`ebam$diagnostics$n_permutations` record the requested and actual
counts.

## Match a group template

To find glycoforms with higher abundance in group C relative to the
other three groups pooled, assign C a score of 1 and the others a score
of 0. Each sample takes its own group’s score. The signed correlation
measures the direction and strength of the match on the analysis scale.

``` r

exp |>
  gly_pattern_hunter(template = c(H = 0, M = 0, Y = 0, C = 1), add_info = FALSE) |>
  get_tidy_result() |>
  filter(cor > 0) |>
  arrange(p_adj, desc(cor)) |>
  head()
#> # A tibble: 6 × 7
#>   .variable                         cor n_samples n_groups n_scores  p_val p_adj
#>   <chr>                           <dbl>     <int>    <int>    <int>  <dbl> <dbl>
#> 1 P01877-92-Hex(3)HexNAc(5)dHex(… 0.604        12        4        2 0.0376 0.790
#> 2 P43652-109-Hex(5)HexNAc(4)NeuA… 0.595        12        4        2 0.0415 0.790
#> 3 P10909-291-Hex(3)HexNAc(8)      0.575        12        4        2 0.0503 0.790
#> 4 P10909-291-Hex(6)HexNAc(6)NeuA… 0.567        12        4        2 0.0543 0.790
#> 5 P04196-344-Hex(6)HexNAc(5)dHex… 0.553        12        4        2 0.0620 0.790
#> 6 O75882-1082-Hex(6)HexNAc(4)dHe… 0.543        12        4        2 0.0680 0.790
```

The sample and group counts show how much of the template each glycoform
covers. The default BH correction covers the variables tested for this
template.

## Map glycoform profiles across samples

With the optional `kohonen` package installed,
[`gly_som()`](https://glycoverse.github.io/glystats/dev/reference/gly_som.md)
places each glycoform on a two-dimensional grid using its abundance in
each sample. Glycoforms assigned to the same unit share a fitted
prototype. Use `on = "sample"` to map samples across glycoforms. To
explore group-mean curves, calculate those means first and supply one
assay column per group.

``` r

som <- gly_som(exp, grid = c(3, 2), seed = 42, add_info = FALSE)
head(get_tidy_result(som))
#> # A tibble: 6 × 7
#>   .variable                         unit     x     y distance n_finite exclusion
#>   <chr>                            <int> <int> <int>    <dbl>    <int> <chr>    
#> 1 P08185-176-Hex(5)HexNAc(4)NeuAc…     6     3     2    2.28        12 NA       
#> 2 P04196-344-Hex(5)HexNAc(4)NeuAc…     1     1     1    1.23        12 NA       
#> 3 P04196-344-Hex(5)HexNAc(4)           3     3     1    2.00        12 NA       
#> 4 P10909-291-Hex(6)HexNAc(5)           3     3     1    0.843       12 NA       
#> 5 P04196-344-Hex(5)HexNAc(4)NeuAc…     1     1     1    2.92        12 NA       
#> 6 P04196-345-Hex(5)HexNAc(4)           3     3     1    1.77        12 NA
som$diagnostics$units
#> # A tibble: 6 × 4
#>    unit     x     y n_profiles
#>   <int> <int> <int>      <int>
#> 1     1     1     1         22
#> 2     2     2     1         14
#> 3     3     3     1         57
#> 4     4     1     2         43
#> 5     5     2     2         63
#> 6     6     3     2         26
```

The `distance` column measures how closely each profile fits its
assigned prototype. An excluded profile keeps its ID and an `exclusion`
reason. The unit table counts the profiles assigned to each cell,
including empty cells.

By default each sample coordinate is standardized across the training
glycoforms after the abundance transformation. With `scale = FALSE`,
coordinates with a wider spread carry more weight. For relative
abundances, similarity describes variation relative to the normalization
reference. These assignments are exploratory. The seed reproduces a fit
for the same input and settings, while a different seed or grid size can
change membership.

## Available analyses

The package includes the following analyses for glycomics and
glycoproteomics data:

- **Differential expression analysis:**
  - [`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md):
    Two-sample t-test
  - [`gly_wilcox()`](https://glycoverse.github.io/glystats/dev/reference/gly_wilcox.md):
    Wilcoxon rank sum test
  - [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md):
    One-way ANOVA
  - [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md):
    Kruskal-Wallis rank sum test
  - [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md):
    Linear models for microarray data (limma)
  - [`gly_sam()`](https://glycoverse.github.io/glystats/dev/reference/gly_sam.md)
    compares two independent groups using SAM.
  - [`gly_ebam()`](https://glycoverse.github.io/glystats/dev/reference/gly_ebam.md)
    estimates posterior probabilities of differential abundance between
    two independent groups.
  - [`gly_ancova()`](https://glycoverse.github.io/glystats/dev/reference/gly_ancova.md):
    ANCOVA (Analysis of Covariance)
  - [`gly_linear_model()`](https://glycoverse.github.io/glystats/dev/reference/gly_linear_model.md):
    Formula-based moderated linear models
  - [`gly_set_test()`](https://glycoverse.github.io/glystats/dev/reference/gly_set_test.md):
    correlated variable-set testing
  - [`gly_fold_change()`](https://glycoverse.github.io/glystats/dev/reference/gly_fold_change.md):
    Calculate fold change
- **Dimensionality reduction:**
  - [`gly_pca()`](https://glycoverse.github.io/glystats/dev/reference/gly_pca.md):
    Principal component analysis
  - [`gly_tsne()`](https://glycoverse.github.io/glystats/dev/reference/gly_tsne.md):
    t-distributed stochastic neighbor embedding (t-SNE)
  - [`gly_umap()`](https://glycoverse.github.io/glystats/dev/reference/gly_umap.md):
    Uniform manifold approximation and projection (UMAP)
  - [`gly_oplsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_oplsda.md):
    Orthogonal partial least squares discriminant analysis (OPLS-DA)
  - [`gly_plsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_plsda.md):
    Partial least squares discriminant analysis (PLS-DA)
- **Clustering:**
  - [`gly_kmeans()`](https://glycoverse.github.io/glystats/dev/reference/gly_kmeans.md):
    K-means clustering
  - [`gly_hclust()`](https://glycoverse.github.io/glystats/dev/reference/gly_hclust.md):
    Hierarchical clustering
  - [`gly_som()`](https://glycoverse.github.io/glystats/dev/reference/gly_som.md)
    maps glycan or sample profiles to a self-organizing grid.
- **Survival analysis:**
  - [`gly_cox()`](https://glycoverse.github.io/glystats/dev/reference/gly_cox.md):
    Cox proportional hazards model
- **Additional tools:**
  - [`gly_cor()`](https://glycoverse.github.io/glystats/dev/reference/gly_cor.md):
    Correlation analysis
  - [`gly_pattern_hunter()`](https://glycoverse.github.io/glystats/dev/reference/gly_pattern_hunter.md)
    matches variables to a named group template.
  - [`gly_roc()`](https://glycoverse.github.io/glystats/dev/reference/gly_roc.md):
    Receiver operating characteristic (ROC) analysis
