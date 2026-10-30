# Uniform Manifold Approximation and Projection (UMAP)

Perform UMAP dimensionality reduction on the expression data. The
function uses
[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
to perform UMAP analysis.

## Usage

``` r
gly_umap(
  exp,
  n_neighbors = 15,
  n_components = 2,
  add_info = TRUE,
  ...,
  log_transform = NULL
)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment` containing an expression
  matrix and sample information.

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

- log_transform:

  `NULL` (default) uses `metadata(exp)$scale` to avoid logging data
  twice. `TRUE` forces the method's log transformation on the supplied
  values, even if they are already logged; `FALSE` uses the supplied
  values unchanged. This argument must be named when it follows `...`.

## Value

A list with four elements:

- `tidy_result`: A tibble with UMAP coordinates containing the following
  columns:

  - `sample`: Sample name

  - `umap1`: First UMAP dimension

  - `umap2`: Second UMAP dimension

  - `umap3`, `umap4`, etc.: Additional UMAP dimensions (if n_components
    \> 2)

- `raw_result`: The raw UMAP result matrix

- `meta_data`: A list containing metadata from the input experiment

- `preprocessing`: Input and analysis scales and the applied
  transformation; see
  [gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md).
  The list has classes `glystats_umap_res` and `glystats_res`.

## Required packages

This function requires the `uwot` package to be installed for UMAP
analysis.

## See also

[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)

[gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md)
for input scales, transformations, and fold changes.
