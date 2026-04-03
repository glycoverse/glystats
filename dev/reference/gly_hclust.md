# Hierarchical Clustering for Glycomics and Glycoproteomics Data

Perform hierarchical clustering on the expression data. The function
uses [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) to perform
clustering and provides tidy results including cluster assignments,
dendrogram data for plotting, and merge heights.

## Usage

``` r
gly_hclust(
  exp,
  on = "variable",
  k_values = c(2, 3, 4, 5),
  scale = TRUE,
  add_info = TRUE,
  ...
)

gly_hclust_(
  expr_mat,
  on = "variable",
  k_values = c(2, 3, 4, 5),
  scale = TRUE,
  ...
)
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

- k_values:

  A numeric vector specifying the number of clusters to cut the tree
  into. Default is c(2, 3, 4, 5). If NULL, no cluster assignments are
  returned.

- scale:

  A logical indicating whether to scale the data before clustering.
  Default is TRUE.

- add_info:

  A logical value. If TRUE (default), sample information from the
  experiment will be added to the result tibbles. If FALSE, only the
  clustering results are returned. Only applicable to `gly_hclust()`.

- ...:

  Additional arguments passed to
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html) and
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Note: if
  both functions need a `method` parameter, use `dist.method` for
  distance and `hclust.method` for clustering method.

- expr_mat:

  (Only for `gly_hclust_()`) A numeric matrix with variables as rows and
  samples as columns.

## Value

A list containing:

- `tidy_result`: A list of tibbles with clustering results:

  - `clusters`: Cluster assignments containing the following columns:

    - `variable` or `sample`: Variable or sample name (depending on `on`
      parameter)

    - `cluster_k2`, `cluster_k3`, etc.: Cluster assignments for
      different k values

  - `dendrogram`: Dendrogram segment data for plotting (if ggdendro is
    available) containing:

    - `x`, `y`, `xend`, `yend`: Segment coordinates for plotting

  - `heights`: Merge heights and steps containing the following columns:

    - `merge_step`: Step number in the clustering process

    - `height`: Height at which clusters are merged

    - `n_clusters`: Number of clusters remaining after this merge

  - `labels`: Labels and their positions (if ggdendro is available)
    containing:

    - `x`, `y`: Position coordinates

    - `label`: Label text

- `raw_result`: The raw hclust object from
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)

- `meta_data`: A list containing metadata from the input experiment

## Details

The function performs log2 transformation on the expression data
(log2(x + 1)) before clustering. When `on = "variable"` (default),
variables are clustered based on their expression patterns across
samples. When `on = "sample"`, samples are clustered based on their
expression profiles across variables.

`gly_hclust()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_hclust_()` is the underlying API that works with matrices directly,
providing more flexibility for users who don't use the glyexp package.

**Distance Calculation:** Distance is calculated using
[`stats::dist()`](https://rdrr.io/r/stats/dist.html) with the specified
method.

**Clustering Method:** Hierarchical clustering is performed using
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) with the
specified method.

**Cluster Assignment:** The dendrogram is cut at different heights to
produce cluster assignments for the specified k values using
[`stats::cutree()`](https://rdrr.io/r/stats/cutree.html).

## Required packages

This function only uses base R packages and does not require additional
dependencies. For enhanced dendrogram plotting capabilities, the
`ggdendro` package is recommended but not required.

## See also

[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html),
[`stats::dist()`](https://rdrr.io/r/stats/dist.html),
[`stats::cutree()`](https://rdrr.io/r/stats/cutree.html)
