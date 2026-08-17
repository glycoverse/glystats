# Construct Correlated Variable Sets

Group highly correlated variables into sets for covariance-aware
differential testing with
[`gly_set_test()`](https://glycoverse.github.io/glystats/dev/reference/gly_set_test.md).

## Usage

``` r
gly_correlated_sets(
  exp,
  threshold = 0.9,
  correlation = "pearson",
  clustering = "complete",
  min_size = 2L,
  within = NULL
)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment` containing expression data
  and variable information.

- threshold:

  A numeric correlation threshold between 0 and 1. Variables are linked
  only when their correlation is strictly greater than this value.
  Default is `0.9`.

- correlation:

  Correlation coefficient to use. Either `"pearson"` (default) or
  `"spearman"`.

- clustering:

  How correlated variables are grouped. `"complete"` (default) uses
  complete-linkage clustering so every pair in a set exceeds
  `threshold`. `"connected"` uses connected components and therefore
  allows transitive chaining.

- min_size:

  Minimum number of distinct profiles in a set. Default is 2.

- within:

  Optional variable-information columns defining strata within which
  sets are constructed, such as `c("protein", "protein_site")`.

## Value

An object of class `glystats_correlated_sets` containing:

- `sets`: A named list mapping set identifiers to all member variables.

- `membership`: A tibble containing set, variable, representative, and
  alias information.

- `representatives`: The distinct profiles used to define each set.

- `correlation_matrix`: The pooled-sample correlation matrix.
  Correlations across different `within` strata are `NA`.

- `excluded_variables`: A tibble of excluded variables and reasons.

- `aliases`: A named list mapping representative variables to aliases.

- Construction settings and metadata copied from `exp`.

## Details

Correlations are calculated on `log2(x + 1e-6)` values, matching the
scale used by
[`gly_set_test()`](https://glycoverse.github.io/glystats/dev/reference/gly_set_test.md)
and other differential-testing functions in this package. All-missing,
insufficiently observed, and zero-variance variables are excluded.
Numerically equivalent profiles are represented once during clustering
and testing, while their aliases remain in the returned set membership.

`clustering = "connected"` reproduces the transitive set-construction
rule used by glycowork. `clustering = "complete"` avoids chaining by
requiring all within-set correlations to exceed `threshold`.

## See also

[`gly_set_test()`](https://glycoverse.github.io/glystats/dev/reference/gly_set_test.md),
[`stats::cor()`](https://rdrr.io/r/stats/cor.html),
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)

## Examples

``` r
set.seed(1)
expression <- matrix(
  stats::rexp(36),
  nrow = 3,
  dimnames = list(c("A", "B", "C"), paste0("S", 1:12))
)
expression["B", ] <- expression["A", ] * exp(stats::rnorm(12, sd = 0.01))
exp <- SummarizedExperiment::SummarizedExperiment(
  assays = list(expression = expression),
  colData = S4Vectors::DataFrame(
    group = factor(rep(c("control", "case"), each = 6))
  )
)
sets <- gly_correlated_sets(exp, threshold = 0.8)
sets$membership
#> # A tibble: 2 × 4
#>   set_id variable representative is_alias
#>   <chr>  <chr>    <chr>          <lgl>   
#> 1 set_1  A        A              FALSE   
#> 2 set_1  B        B              FALSE   
```
