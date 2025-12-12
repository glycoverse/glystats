# Getting Started with glystats

When we talk about omics data analysis, we’re usually diving into
exciting territories like differential expression analysis, PCA,
survival analysis, and much more! 🧬 `glystats` brings all these
powerful analyses together under one unified, user-friendly interface.

**🎯 Key Feature:** `glystats` offers two levels of interfaces to fit
your workflow:

- `gly_xxx()`: seamlessly works with \[glyexp::experiment()\], the
  beating heart of the `glycoverse` ecosystem 💓
- `gly_xxx_()`: flexible enough to work with general inputs like
  matrices or data frames

This vignette focuses on the `gly_xxx()` interface. New to
\[glyexp::experiment()\]? No worries! Check out its
[introduction](https://glycoverse.github.io/glyexp/articles/glyexp.html)
first to get up to speed.

``` r
library(glystats)
library(glyexp)
library(glyread)
library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:stats':
#> 
#>     aggregate
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:glyexp':
#> 
#>     select_var
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
#> ✔ Finding leader proteins [84ms]
#> 
#> ℹ Reading dataColumn group converted to <factor>.ℹ Parsing glycan compositions and structures
#> Column group converted to <factor>.✔ Parsing glycan compositions and structures [270ms]
#> 
#> ℹ Reading data✔ Reading data [764ms]
#> 
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
#> ℹ Total removed: 2 (0.67%) variables.
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
#> ℹ Aggregating to "gf" level
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
exp
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 225 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, glycan_composition <comp>, protein_site <int>, gene <chr>
```

Look at that! 🎉 We’ve got a \[glyexp::experiment()\] packed with 12
samples and 263 glycoforms. That’s plenty of data to work with!

``` r
get_var_info(exp)
#> # A tibble: 225 × 5
#>    variable protein glycan_composition      protein_site gene    
#>    <chr>    <chr>   <comp>                         <int> <chr>   
#>  1 V1       P08185  Hex(5)HexNAc(4)NeuAc(2)          176 SERPINA6
#>  2 V2       P04196  Hex(5)HexNAc(4)NeuAc(1)          344 HRG     
#>  3 V3       P04196  Hex(5)HexNAc(4)                  344 HRG     
#>  4 V4       P10909  Hex(6)HexNAc(5)                  291 CLU     
#>  5 V5       P04196  Hex(5)HexNAc(4)NeuAc(2)          344 HRG     
#>  6 V6       P04196  Hex(5)HexNAc(4)                  345 HRG     
#>  7 V7       P04196  Hex(5)HexNAc(4)dHex(2)           344 HRG     
#>  8 V8       P04196  Hex(4)HexNAc(3)                  344 HRG     
#>  9 V9       P04196  Hex(4)HexNAc(4)NeuAc(1)          344 HRG     
#> 10 V10      P10909  Hex(5)HexNAc(4)                  291 CLU     
#> # ℹ 215 more rows
```

The variable information tibble is like a detailed ID card for each
glycoform 🆔 - it contains everything you need to know: the protein,
glycosylation site, and glycan structures.

``` r
get_sample_info(exp)
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
[`gly_ttest()`](https://glycoverse.github.io/glystats/reference/gly_ttest.md)
for t-tests,
[`gly_pca()`](https://glycoverse.github.io/glystats/reference/gly_pca.md)
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
[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
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
[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md),
the `tidy_result` contains two informative tibbles: `main_test` and
`post_hoc_test`.

You can use
[`get_tidy_result()`](https://glycoverse.github.io/glystats/reference/get_tidy_result.md)
to get the tidy result tibble:

``` r
get_tidy_result(anova_res, "main_test")
#> # A tibble: 225 × 13
#>    variable term     df sumsq meansq statistic   p_val  p_adj post_hoc   protein
#>    <chr>    <chr> <dbl> <dbl>  <dbl>     <dbl>   <dbl>  <dbl> <chr>      <chr>  
#>  1 V1       group     3 67.7   22.6      6.68  0.0143  0.0528 NA         P08185 
#>  2 V10      group     3 87.0   29.0      7.87  0.00903 0.0391 C_vs_H;H_… P10909 
#>  3 V100     group     3  6.57   2.19     2.59  0.125   0.228  NA         P01871 
#>  4 V101     group     3 15.4    5.12     6.13  0.0180  0.0580 NA         P01871 
#>  5 V102     group     3 43.5   14.5      0.520 0.680   0.758  NA         P01871 
#>  6 V103     group     3  8.61   2.87     3.93  0.0539  0.124  NA         P01871 
#>  7 V104     group     3  5.30   1.77     0.347 0.793   0.829  NA         P01871 
#>  8 V105     group     3  3.28   1.09     6.08  0.0185  0.0581 NA         P01871 
#>  9 V106     group     3 14.0    4.66     6.59  0.0148  0.0534 NA         P01871 
#> 10 V107     group     3  7.43   2.48     1.51  0.286   0.409  NA         P01871 
#> # ℹ 215 more rows
#> # ℹ 3 more variables: glycan_composition <comp>, protein_site <int>, gene <chr>
```

Notice something cool? 😎
[`gly_anova()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
thoughtfully adds back all the descriptive columns from your variable
tibble. Want to control this behavior? Just use the `add_info`
parameter!

The `raw_result` houses two lists of models - one for the main test, one
for post hoc comparisons:

``` r
names(get_raw_result(anova_res))
#> [1] "main_test"     "post_hoc_test"
```

[`get_tidy_result()`](https://glycoverse.github.io/glystats/reference/get_tidy_result.md)
and
[`get_raw_result()`](https://glycoverse.github.io/glystats/reference/get_tidy_result.md)
are useful to be used in pipes:

``` r
exp |>
  gly_anova() |>
  get_tidy_result("main_test") |>
  filter(p_adj < 0.05)
#> ℹ Number of groups: 4
#> ℹ Groups: "C", "H", "M", and "Y"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
#> # A tibble: 60 × 13
#>    variable term     df  sumsq meansq statistic   p_val   p_adj post_hoc protein
#>    <chr>    <chr> <dbl>  <dbl>  <dbl>     <dbl>   <dbl>   <dbl> <chr>    <chr>  
#>  1 V10      group     3  87.0   29.0       7.87 9.03e-3 3.91e-2 C_vs_H;… P10909 
#>  2 V11      group     3 176.    58.8      20.8  3.94e-4 5.90e-3 C_vs_H;… P04196 
#>  3 V113     group     3 226.    75.2      46.6  2.07e-5 5.17e-4 C_vs_H;… P10909 
#>  4 V116     group     3   9.72   3.24      8.08 8.36e-3 3.69e-2 C_vs_H;… P01871 
#>  5 V12      group     3  81.1   27.0       9.70 4.85e-3 2.66e-2 C_vs_H;… P13671 
#>  6 V125     group     3  11.0    3.68      9.72 4.81e-3 2.66e-2 C_vs_H;… P10909 
#>  7 V126     group     3   9.36   3.12     11.1  3.22e-3 2.41e-2 C_vs_H;… P02790 
#>  8 V13      group     3 159.    52.9     138.   3.12e-7 1.75e-5 C_vs_H;… P04196 
#>  9 V130     group     3 308.   103.      144.   2.66e-7 1.75e-5 C_vs_H;… P01860 
#> 10 V131     group     3 272.    90.7      11.2  3.12e-3 2.41e-2 H_vs_M;… P01860 
#> # ℹ 50 more rows
#> # ℹ 3 more variables: glycan_composition <comp>, protein_site <int>, gene <chr>
```

## 🔧 Maximum Flexibility Mode

Feeling constrained by \[glyexp::experiment()\]? Fear not! 🦸‍♀️ Every
`gly_xxx()` function comes with a flexible `gly_xxx_()` sibling that
works with standard R objects. For instance,
[`gly_anova_()`](https://glycoverse.github.io/glystats/reference/gly_anova.md)
happily accepts plain matrices:

``` r
expr_mat <- get_expr_mat(exp)
groups <- factor(get_sample_info(exp)$group)
anova_res2 <- gly_anova_(expr_mat, groups)
```

This adaptability makes `glystats` a perfect team player 🤝 - it
seamlessly integrates into existing analysis pipelines and workflows, no
matter what data structures you’re already using!

## 🎪 The Complete Analytical Arsenal

Ready to explore the full power of `glystats`? Here’s your complete
toolkit for glycomics and glycoproteomics data analysis:

- **🔬 Differential Expression Analysis:**
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
  - [`gly_fold_change()`](https://glycoverse.github.io/glystats/reference/gly_fold_change.md):
    Calculate fold change
- **📐 Dimensionality Reduction:**
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
- **🧩 Clustering:**
  - [`gly_kmeans()`](https://glycoverse.github.io/glystats/reference/gly_kmeans.md):
    K-means clustering
  - [`gly_hclust()`](https://glycoverse.github.io/glystats/reference/gly_hclust.md):
    Hierarchical clustering
  - [`gly_consensus_clustering()`](https://glycoverse.github.io/glystats/reference/gly_consensus_clustering.md):
    Consensus clustering
  - [`gly_wgcna()`](https://glycoverse.github.io/glystats/reference/gly_wgcna.md):
    Weighted gene co-expression network analysis (WGCNA)
- **⏱️ Survival Analysis:**
  - [`gly_cox()`](https://glycoverse.github.io/glystats/reference/gly_cox.md):
    Cox proportional hazards model
- **🎯 Enrichment Analysis:**
  - [`gly_enrich_go()`](https://glycoverse.github.io/glystats/reference/gly_enrich_go.md):
    Gene Ontology enrichment analysis
  - [`gly_enrich_kegg()`](https://glycoverse.github.io/glystats/reference/gly_enrich_go.md):
    KEGG enrichment analysis
  - [`gly_enrich_reactome()`](https://glycoverse.github.io/glystats/reference/gly_enrich_go.md):
    Reactome enrichment analysis
- **🔧 Additional Tools:**
  - [`gly_cor()`](https://glycoverse.github.io/glystats/reference/gly_cor.md):
    Correlation analysis
  - [`gly_roc()`](https://glycoverse.github.io/glystats/reference/gly_roc.md):
    Receiver operating characteristic (ROC) analysis

## 🚀 What’s Next on Your Journey?

Ready to dive deeper into the `glycoverse`? Here’s your roadmap to
success:

1.  **📥 Data Import:** Start with
    [glyread](https://glycoverse.github.io/glyread/articles/glyread.html)
    to seamlessly import your data into
    [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
    objects

2.  **🧹 Data Preprocessing:** Use
    [glyclean](https://glycoverse.github.io/glyclean/articles/glyclean.html)
    to polish and prepare your data for analysis

3.  **📊 Statistical Analysis:** You’re here! Use `glystats` to unlock
    powerful insights from your glycomics data

4.  **🎨 Visualization:** Stay tuned! We’re crafting an amazing `glyvis`
    package for stunning data visualizations

Happy analyzing! 🎉✨
