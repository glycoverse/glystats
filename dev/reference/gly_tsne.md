# t-Distributed Stochastic Neighbor Embedding (t-SNE)

Perform t-SNE dimensionality reduction on the expression data. The
function uses
[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html) to perform
t-SNE analysis.

## Usage

``` r
gly_tsne(exp, dims = 2, perplexity = 30, add_info = TRUE, ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- dims:

  Number of output dimensions. Default is 2.

- perplexity:

  Perplexity parameter for t-SNE. Default is 30.

- add_info:

  A logical value. If TRUE (default), sample information from the
  experiment will be added to the result tibble. If FALSE, only the
  t-SNE coordinates are returned.

- ...:

  Additional arguments passed to
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html).

## Value

A list with three elements:

- `tidy_result`: A tibble with t-SNE coordinates containing the
  following columns:

  - `sample`: Sample name

  - `tsne1`: First t-SNE dimension

  - `tsne2`: Second t-SNE dimension

- `raw_result`: The raw Rtsne object

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_tsne_res` and `glystats_res`.

## Required packages

This function requires the `Rtsne` package to be installed for t-SNE
analysis.

## See also

[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html)
