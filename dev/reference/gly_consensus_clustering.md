# Consensus Clustering for Glycomics and Glycoproteomics Data

**\[deprecated\]**

## Usage

``` r
gly_consensus_clustering(
  exp,
  on = "sample",
  max_k = 9,
  reps = 1000,
  p_item = 0.8,
  cluster_alg = "hc",
  scale = TRUE,
  output_file = NULL,
  add_info = TRUE,
  ...
)

gly_consensus_clustering_(
  expr_mat,
  on = "sample",
  max_k = 9,
  reps = 1000,
  p_item = 0.8,
  cluster_alg = "hc",
  scale = TRUE,
  output_file = NULL,
  ...
)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- on:

  A character string specifying what to cluster. Either "sample"
  (default) to cluster samples/observations, or "variable" to cluster
  variables/features.

- max_k:

  Maximum number of clusters to test. Default is 9. Passed to the `maxK`
  argument of
  [`ConsensusClusterPlus::ConsensusClusterPlus()`](https://rdrr.io/pkg/ConsensusClusterPlus/man/ConsensusClusterPlus.html).

- reps:

  Number of resampling iterations. Default is 1000. Passed to the `reps`
  argument of
  [`ConsensusClusterPlus::ConsensusClusterPlus()`](https://rdrr.io/pkg/ConsensusClusterPlus/man/ConsensusClusterPlus.html).

- p_item:

  Proportion of items to sample. Default is 0.8. Passed to the `pItem`
  argument of
  [`ConsensusClusterPlus::ConsensusClusterPlus()`](https://rdrr.io/pkg/ConsensusClusterPlus/man/ConsensusClusterPlus.html).

- cluster_alg:

  Clustering algorithm to use. Default is "hc" (hierarchical
  clustering). Other options include "km" (k-means), "pam" (partitioning
  around medoids). Passed to the `clusterAlg` argument of
  [`ConsensusClusterPlus::ConsensusClusterPlus()`](https://rdrr.io/pkg/ConsensusClusterPlus/man/ConsensusClusterPlus.html).

- scale:

  A logical indicating whether to scale the data. Default is TRUE.

- output_file:

  A character string specifying the output file path for plots. If NULL
  (default), no plots are saved and no files are generated. If provided,
  the consensus clustering plot will be saved directly to this file
  path. The function will create any necessary parent directories.

- add_info:

  A logical value. If TRUE (default), sample and variable information
  from the experiment will be added to the result tibbles. If FALSE,
  only the consensus clustering results are returned. Only applicable to
  `gly_consensus_clustering()`.

- ...:

  Additional arguments passed to
  [`ConsensusClusterPlus::ConsensusClusterPlus()`](https://rdrr.io/pkg/ConsensusClusterPlus/man/ConsensusClusterPlus.html).

- expr_mat:

  (Only for `gly_consensus_clustering_()`) A numeric matrix with
  variables as rows and samples as columns.

## Value

A list with three elements:

- `tidy_result`: A tibble with consensus clustering results containing
  the following columns:

  - `variable` or `sample`: Variable or sample name (depending on `on`
    parameter)

  - `k`: Number of clusters

  - `cluster`: Cluster assignment for the corresponding k

- `raw_result`: The raw ConsensusClusterPlus object

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_cc_res` and `glystats_res`.

## Details

Perform consensus clustering on the expression data using
ConsensusClusterPlus. The function returns cluster assignments for all k
values from 2 to max_k in long format. Setting `output_file` to
visualize consensus CDF curves and consensus matrices, which helps you
decide the optimal number of clusters. See [this
tutorial](https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf)
for more information.

The function performs log2 transformation on the expression data
(log2(x + 1)) before consensus clustering. When `on = "sample"`
(default), samples are clustered based on their expression profiles
across variables. When `on = "variable"`, variables are clustered based
on their expression patterns across samples.

`gly_consensus_clustering()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_consensus_clustering_()` is the underlying API that works with
matrices directly, providing more flexibility for users who don't use
the glyexp package.

**Consensus Clustering Process:**

1.  Perform consensus clustering for k = 2 to max_k using
    ConsensusClusterPlus

2.  Return cluster assignments for all k values in long format

## Required packages

This function requires the following packages to be installed:

- `ConsensusClusterPlus` for consensus clustering

## See also

[`ConsensusClusterPlus::ConsensusClusterPlus()`](https://rdrr.io/pkg/ConsensusClusterPlus/man/ConsensusClusterPlus.html)
