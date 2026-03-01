# Uniform Manifold Approximation and Projection (UMAP)

Perform UMAP dimensionality reduction on the expression data. The
function uses
[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
to perform UMAP analysis.

## Usage

``` r
gly_umap(exp, n_neighbors = 15, n_components = 2, add_info = TRUE, ...)

gly_umap_(expr_mat, n_neighbors = 15, n_components = 2, ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- n_neighbors:

  Number of neighbors to consider for each point. Default is 15.

- n_components:

  Number of output dimensions. Default is 2.

- add_info:

  A logical value. If TRUE (default), sample information from the
  experiment will be added to the result tibble. If FALSE, only the UMAP
  coordinates are returned. Only applicable to `gly_umap()`.

- ...:

  Additional arguments passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

- expr_mat:

  (Only for `gly_umap_()`) A numeric matrix with variables as rows and
  samples as columns.

## Value

A list with two elements:

- `tidy_result`: A tibble with UMAP coordinates containing the following
  columns:

  - `sample`: Sample name

  - `umap1`: First UMAP dimension

  - `umap2`: Second UMAP dimension

  - `umap3`, `umap4`, etc.: Additional UMAP dimensions (if n_components
    \> 2)

- `raw_result`: The raw UMAP result matrix The list has classes
  `glystats_umap_res` and `glystats_res`.

## Details

`gly_umap()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_umap_()` is the underlying API that works with matrices directly,
providing more flexibility for users who don't use the glyexp package.

## Required packages

This function requires the `uwot` package to be installed for UMAP
analysis.

## See also

[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
