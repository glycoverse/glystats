# Getting Started with glystats

When we talk about omics data analysis, we’re usually diving into
exciting territories like differential expression analysis, PCA,
survival analysis, and much more! 🧬 `glystats` brings all these
powerful analyses together under one unified, user-friendly interface.

**🎯 Key Feature:** Every `glystats` analysis accepts
\[glyexp::GlycomicSE()\] and \[glyexp::GlycoproteomicSE()\] containers,
the unified data interfaces at the heart of the `glycoverse` ecosystem
💓

New to these containers? No worries! Check out the [glyexp
introduction](https://glycoverse.github.io/glyexp/articles/glyexp.html)
first to get up to speed.

``` r

library(glystats)
library(glyexp)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> The following object is masked from 'package:glyexp':
#> 
#>     samples
library(glyread)
library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:S4Vectors':
#> 
#>     aggregate
#> The following object is masked from 'package:stats':
#> 
#>     aggregate
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following objects are masked from 'package:GenomicRanges':
#> 
#>     intersect, setdiff, union
#> The following object is masked from 'package:Seqinfo':
#> 
#>     intersect
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from 'package:generics':
#> 
#>     explain
#> The following object is masked from 'package:matrixStats':
#> 
#>     count
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## 🔍 Meet Your Data

Let’s start by exploring our demo dataset:

``` r

# Here we use `glyread` to read in the data and use `glyclean` to perform preprocessing
exp <- read_pglyco3_pglycoquant("glycopeptides.list", sample_info = "sample_info.csv") |> auto_clean()
#> ℹ Reading data
#> ℹ Finding leader proteins
#> ✔ Finding leader proteins [80ms]
#> 
#> ℹ Reading dataℹ Parsing glycan compositions and structures
#> ✔ Parsing glycan compositions and structures [325ms]
#> 
#> ℹ Reading data✔ Reading data [758ms]
#> 
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 2 (0.67%) variables.
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
#> ℹ Aggregating to "gf" level
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
exp
#> 
#> ── GlycoproteomicSE ────────────────────────────────────────────────────────────
#> ℹ Abundance assay: 12 samples, 225 variables
#> ℹ Glycan type: N
#> ℹ Row data fields: protein <chr>, glycan_composition <comp>, protein_site <int>, gene <chr>
#> ℹ Column data fields: group <fct>
#> ℹ Metadata fields: exp_type <chr>, glycan_type <chr>, quant_method <chr>
```

Look at that! 🎉 We’ve got a \[glyexp::GlycoproteomicSE()\] packed with
12 samples and 263 glycoforms. That’s plenty of data to work with!

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

The variable information tibble is like a detailed ID card for each
glycoform 🆔 - it contains everything you need to know: the protein,
glycosylation site, and glycan structures.

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

Our sample information tibble features a crucial “group” column 🏷️.
Here’s the key: “H” = healthy, “M” = hepatitis, “Y” = cirrhosis, and “C”
= hepatocellular carcinoma. This gives us a perfect setup for
comparative analysis!

## 🚀 One Interface to Rule Them All

Here’s where `glystats` really shines! ✨ Every function follows a
simple, intuitive naming pattern: `gly_` + analysis name. Think
[`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md)
for t-tests,
[`gly_pca()`](https://glycoverse.github.io/glystats/dev/reference/gly_pca.md)
for PCA, and so on.

**Pro tip:** leverage RStudio’s auto-completion to discover all
available functions! 💡

Let’s dive into action with an ANOVA analysis to identify differentially
expressed glycoforms:

``` r

anova_res <- gly_anova(exp)
#> ℹ Number of groups: 4
#> ℹ Groups: "C", "H", "M", and "Y"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
```

Boom! 💥 Analysis complete in just one line!
[`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
intelligently follows the `glycoverse` naming conventions, automatically
detecting the “group” column in your sample info and fitting an ANOVA
model for each glycoform.

## 📋 Understanding Your Results

All `glystats` functions return a consistent, well-structured list with
two key components:

- `tidy_result`: Clean, analysis-ready tibbles in tidy format 📊
- `raw_result`: The original result objects from underlying statistical
  functions 🔧

For
[`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md),
the `tidy_result` contains two informative tibbles: `main_test` and
`post_hoc_test`.

You can use
[`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
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

Notice something cool? 😎
[`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
thoughtfully adds back all the descriptive columns from your variable
tibble. Want to control this behavior? Just use the `add_info`
parameter!

The `raw_result` houses two lists of models - one for the main test, one
for post hoc comparisons:

``` r

names(get_raw_result(anova_res))
#> [1] "main_test"     "post_hoc_test"
```

[`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
and
[`get_raw_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
are useful to be used in pipes:

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

## 🎪 The Complete Analytical Arsenal

Ready to explore the full power of `glystats`? Here’s your complete
toolkit for glycomics and glycoproteomics data analysis:

- **🔬 Differential Expression Analysis:**
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
  - [`gly_linear_model()`](https://glycoverse.github.io/glystats/dev/reference/gly_linear_model.md):
    Formula-based moderated linear models
  - [`gly_fold_change()`](https://glycoverse.github.io/glystats/dev/reference/gly_fold_change.md):
    Calculate fold change

For multifactor experiments, use a one-sided formula to describe the
design. Named contrasts can combine model coefficients; wrap interaction
coefficient names in backticks:

``` r

model_res <- gly_linear_model(
  exp,
  ~ treatment * time + batch,
  contrasts = c(
    treatment_at_late = "treatmentB + `treatmentB:timelate`"
  )
)
```

- **📐 Dimensionality Reduction:**
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
- **🧩 Clustering:**
  - [`gly_kmeans()`](https://glycoverse.github.io/glystats/dev/reference/gly_kmeans.md):
    K-means clustering
  - [`gly_hclust()`](https://glycoverse.github.io/glystats/dev/reference/gly_hclust.md):
    Hierarchical clustering
- **⏱️ Survival Analysis:**
  - [`gly_cox()`](https://glycoverse.github.io/glystats/dev/reference/gly_cox.md):
    Cox proportional hazards model
- **🔧 Additional Tools:**
  - [`gly_cor()`](https://glycoverse.github.io/glystats/dev/reference/gly_cor.md):
    Correlation analysis
  - [`gly_roc()`](https://glycoverse.github.io/glystats/dev/reference/gly_roc.md):
    Receiver operating characteristic (ROC) analysis

## 🚀 What’s Next on Your Journey?

Ready to dive deeper into the `glycoverse`? Here’s your roadmap to
success:

1.  **📥 Data Import:** Start with
    [glyread](https://glycoverse.github.io/glyread/articles/glyread.html)
    to import your data into
    [`glyexp::GlycomicSE`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
    or
    [`glyexp::GlycoproteomicSE`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
    objects

2.  **🧹 Data Preprocessing:** Use
    [glyclean](https://glycoverse.github.io/glyclean/articles/glyclean.html)
    to polish and prepare your data for analysis

3.  **📊 Statistical Analysis:** You’re here! Use `glystats` to unlock
    powerful insights from your glycomics data

4.  **🎨 Visualization:** Stay tuned! We’re crafting an amazing `glyvis`
    package for stunning data visualizations

Happy analyzing! 🎉✨
