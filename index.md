# glystats

The goal of glystats is to perform statistical analysis on
glycoproteomics and glycomics data. It works seamlessly with the
[glyexp](https://github.com/glycoverse/glyexp) package.

## Installation

You can install the latest release of glystats from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("glycoverse/glystats@*release")
```

Or install the development version:

``` r
remotes::install_github("glycoverse/glystats")
```

## Documentation

- 🚀 Get started:
  [Here](https://glycoverse.github.io/glystats/articles/glystats.html)
- 📚 Reference:
  [Here](https://glycoverse.github.io/glystats/reference/index.html)

## Role in `glycoverse`

`glystats` is the downstream analysis package in the `glycoverse`
ecosystem. It provides statistical analysis functions for
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects. A common workflow is to use
[glyread](https://github.com/glycoverse/glyread) to import data,
[glyclean](https://github.com/glycoverse/glyclean) to preprocess data,
and then `glystats` to perform statistical analysis.

## Example

Say we already have a preprocessed experiment object called `exp`:

``` r
# Two-sample t-test
ttest_res <- gly_ttest(exp)

# PCA analysis
pca_res <- gly_pca(exp)

# ROC analysis
roc_res <- gly_roc(exp)
```

That’s it! These functions use `glycoverse` column conventions to load
needed data and perform analysis. All functions start with `gly_` to
leverage the auto-completion in RStudio. They accept an
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object, and return analysis result as a tibble or a list of tibbles. See
documentation for each function for more details.
