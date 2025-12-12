# K-means Clustering for Glycomics and Glycoproteomics Data

Perform k-means clustering on the expression data. The function uses
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) to perform
clustering and provides tidy results with cluster assignments.

## Usage

``` r
gly_kmeans(
  exp,
  on = "variable",
  centers = 3,
  scale = TRUE,
  add_info = TRUE,
  ...
)

gly_kmeans_(expr_mat, on = "variable", centers = 3, scale = TRUE, ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- on:

  A character string specifying what to cluster. Either "variable"
  (default) to cluster variables/features, or "sample" to cluster
  samples/observations.

- centers:

  Either the number of clusters (integer) or a set of initial cluster
  centers. Default is 3.

- scale:

  A logical indicating whether to scale the data before clustering.
  Default is TRUE.

- add_info:

  A logical value. If TRUE (default), sample information from the
  experiment will be added to the result tibbles. If FALSE, only the
  clustering results are returned. Only applicable to `gly_kmeans()`.

- ...:

  Additional arguments passed to
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).

- expr_mat:

  A numeric matrix with variables as rows and samples as columns.

## Value

A list with two elements:

- `tidy_result`: A tibble with cluster assignments containing the
  following columns:

  - `variable` or `sample`: Variable or sample name (depending on `on`
    parameter)

  - `cluster`: Cluster assignment

- `raw_result`: The raw kmeans object from
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).

## Details

The function performs log2 transformation on the expression data
(log2(x + 1)) before clustering. When `on = "variable"` (default),
variables are clustered based on their expression patterns across
samples. When `on = "sample"`, samples are clustered based on their
expression profiles across variables.

`gly_kmeans()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_kmeans_()` is the underlying API that works with matrices directly,
providing more flexibility for users who don't use the glyexp package.

**Data Preparation:** Data is log2-transformed and optionally scaled
before clustering.

**Clustering Method:** K-means clustering is performed using
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) with the
specified parameters.

## Required packages

This function only uses base R packages and does not require additional
dependencies.

## See also

[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html)
