# Uniform Manifold Approximation and Projection (UMAP)

Perform UMAP dimensionality reduction on the expression data. The
function uses
[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
to perform UMAP analysis.

## Usage

``` r
gly_umap(exp, n_neighbors = 15, n_components = 2, add_info = TRUE, ...)
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
  coordinates are returned.

- ...:

  Additional arguments passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

A list with three elements:

- `tidy_result`: A tibble with UMAP coordinates containing the following
  columns:

  - `sample`: Sample name

  - `umap1`: First UMAP dimension

  - `umap2`: Second UMAP dimension

  - `umap3`, `umap4`, etc.: Additional UMAP dimensions (if n_components
    \> 2)

- `raw_result`: The raw UMAP result matrix

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_umap_res` and `glystats_res`.

## Required packages

This function requires the `uwot` package to be installed for UMAP
analysis.

## See also

[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
