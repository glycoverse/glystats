# Principal Component Analysis (PCA)

Perform principal component analysis on the expression data. The
function uses [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) to
perform PCA and
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) to
tidy the results. If `scale = TRUE`, constant variables (zero variance)
will be removed before PCA.

## Usage

``` r
gly_pca(exp, center = TRUE, scale = TRUE, add_info = TRUE, ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- center:

  A logical indicating whether to center the data. Default is TRUE.

- scale:

  A logical indicating whether to scale the data. Default is TRUE.

- add_info:

  A logical value. If TRUE (default), sample and variable information
  from the experiment will be added to the result tibbles. If FALSE,
  only the PCA results are returned.

- ...:

  Additional arguments passed to
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html).

## Value

A list containing:

- `tidy_result`: A list of tibbles with PCA results:

  - `samples`: PCA scores for each sample containing the following
    columns:

    - `sample`: Sample name

    - `PC`: Principal component name (PC1, PC2, etc.)

    - `value`: Score value for the principal component

  - `variables`: PCA loadings for each variable containing the following
    columns:

    - `variable`: Variable name

    - `PC`: Principal component name (PC1, PC2, etc.)

    - `value`: Loading value for the principal component

  - `eigenvalues`: PCA eigenvalues containing the following columns:

    - `PC`: Principal component name (PC1, PC2, etc.)

    - `std.dev`: Standard deviation

    - `percent`: Percentage of variance explained

    - `cumulative`: Cumulative percentage of variance explained

- `raw_result`: The raw prcomp object from
  [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)

- `meta_data`: A list containing metadata from the input experiment

## Details

The function performs log transformation on the expression data (log(x +
1)) before PCA analysis.

## Required packages

This function only uses base R packages and does not require additional
dependencies.

## See also

[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)
