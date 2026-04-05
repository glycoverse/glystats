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
#> Column group converted to <factor>.✔ Parsing glycan compositions and structures [353ms]
#> 
#> ℹ Reading data✔ Reading data [833ms]
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
#>    variable                        protein glycan_composition protein_site gene 
#>    <glue>                          <chr>   <comp>                    <int> <chr>
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
#>    variable     protein glycan_composition protein_site gene  term     df  sumsq
#>    <glue>       <chr>   <comp>                    <int> <chr> <chr> <dbl>  <dbl>
#>  1 P08185-176-… P08185  Hex(5)HexNAc(4)Ne…          176 SERP… group     3  67.7 
#>  2 P04196-344-… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 161.  
#>  3 P04196-344-… P04196  Hex(5)HexNAc(4)             344 HRG   group     3 126.  
#>  4 P10909-291-… P10909  Hex(6)HexNAc(5)             291 CLU   group     3  23.1 
#>  5 P04196-344-… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 472.  
#>  6 P04196-345-… P04196  Hex(5)HexNAc(4)             345 HRG   group     3  60.0 
#>  7 P04196-344-… P04196  Hex(5)HexNAc(4)dH…          344 HRG   group     3 208.  
#>  8 P04196-344-… P04196  Hex(4)HexNAc(3)             344 HRG   group     3 109.  
#>  9 P04196-344-… P04196  Hex(4)HexNAc(4)Ne…          344 HRG   group     3   9.66
#> 10 P10909-291-… P10909  Hex(5)HexNAc(4)             291 CLU   group     3  87.0 
#> # ℹ 215 more rows
#> # ℹ 5 more variables: meansq <dbl>, statistic <dbl>, p_val <dbl>, p_adj <dbl>,
#> #   post_hoc <chr>
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
#>    variable      protein glycan_composition protein_site gene  term     df sumsq
#>    <glue>        <chr>   <comp>                    <int> <chr> <chr> <dbl> <dbl>
#>  1 P04196-344-H… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 161. 
#>  2 P04196-344-H… P04196  Hex(5)HexNAc(4)             344 HRG   group     3 126. 
#>  3 P04196-344-H… P04196  Hex(5)HexNAc(4)Ne…          344 HRG   group     3 472. 
#>  4 P04196-344-H… P04196  Hex(5)HexNAc(4)dH…          344 HRG   group     3 208. 
#>  5 P10909-291-H… P10909  Hex(5)HexNAc(4)             291 CLU   group     3  87.0
#>  6 P04196-344-H… P04196  Hex(5)HexNAc(4)dH…          344 HRG   group     3 176. 
#>  7 P13671-855-H… P13671  Hex(5)HexNAc(4)dH…          855 C6    group     3  81.1
#>  8 P04196-344-H… P04196  Hex(4)HexNAc(3)dH…          344 HRG   group     3 159. 
#>  9 P04196-344-H… P04196  Hex(5)HexNAc(4)dH…          344 HRG   group     3 150. 
#> 10 P01019-161-H… P01019  Hex(5)HexNAc(4)Ne…          161 AGT   group     3  46.4
#> # ℹ 50 more rows
#> # ℹ 5 more variables: meansq <dbl>, statistic <dbl>, p_val <dbl>, p_adj <dbl>,
#> #   post_hoc <chr>
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
- **⏱️ Survival Analysis:**
  - [`gly_cox()`](https://glycoverse.github.io/glystats/reference/gly_cox.md):
    Cox proportional hazards model
- **🎯 Enrichment Analysis:**
  - [`gly_enrich_go()`](https://glycoverse.github.io/glystats/reference/gly_enrich_go.md):
    Gene Ontology enrichment analysis
  - [`gly_enrich_kegg()`](https://glycoverse.github.io/glystats/reference/gly_enrich_kegg.md):
    KEGG enrichment analysis
  - [`gly_enrich_reactome()`](https://glycoverse.github.io/glystats/reference/gly_enrich_reactome.md):
    Reactome enrichment analysis
  - [`gly_enrich_wikipathways()`](https://glycoverse.github.io/glystats/reference/gly_enrich_wikipathways.md):
    WikiPathways enrichment analysis
  - [`gly_enrich_do()`](https://glycoverse.github.io/glystats/reference/gly_enrich_do.md):
    Disease Ontology enrichment analysis
  - [`gly_enrich_ncg()`](https://glycoverse.github.io/glystats/reference/gly_enrich_ncg.md):
    Network of Cancer Genes enrichment analysis
  - [`gly_enrich_dgn()`](https://glycoverse.github.io/glystats/reference/gly_enrich_dgn.md):
    DisGeNET enrichment analysis
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
