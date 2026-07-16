# glystats

The goal of glystats is to perform statistical analysis on
glycoproteomics and glycomics data. It works seamlessly with the
[glyexp](https://github.com/glycoverse/glyexp) package.

## Installation

### Install glycoverse

We recommend installing the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse), which includes
this package and other core glycoverse packages.

### Install glystats alone

If you don’t want to install all glycoverse packages, you can only
install glystats.

You can install the latest release of glystats from
[r-universe](https://glycoverse.r-universe.dev/glystats)
(**recommended**):

``` r

# install.packages("pak")
pak::repo_add(glycoverse = "https://glycoverse.r-universe.dev")
pak::pkg_install("glystats")
```

Or from [GitHub](https://github.com/glycoverse/glystats):

``` r

pak::pkg_install("glycoverse/glystats@*release")
```

Or install the development version (NOT recommended):

``` r

pak::pkg_install("glycoverse/glystats")
```

**Note:** Tips and troubleshooting for the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse) are also
applicable here: [Installation of
glycoverse](https://github.com/glycoverse/glycoverse#installation).

## Documentation

- 🚀 Get started:
  [Here](https://glycoverse.github.io/glystats/articles/glystats.html)
- 📚 Reference:
  [Here](https://glycoverse.github.io/glystats/reference/index.html)

## Role in `glycoverse`

`glystats` is the downstream analysis package in the `glycoverse`
ecosystem. It provides statistical analysis functions for
[`glyexp::GlycomicSE`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
and
[`glyexp::GlycoproteomicSE`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
objects. A common workflow is to use
[glyread](https://github.com/glycoverse/glyread) to import data,
[glyclean](https://github.com/glycoverse/glyclean) to preprocess data,
and then `glystats` to perform statistical analysis.

## Example

Say we already have a preprocessed glycomics or glycoproteomics
container called `exp`:

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
leverage the auto-completion in RStudio. They accept a
[`glyexp::GlycomicSE`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
or
[`glyexp::GlycoproteomicSE`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
object, and return analysis result as a tibble or a list of tibbles. See
documentation for each function for more details.
