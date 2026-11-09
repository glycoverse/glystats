# Principal Component Analysis (PCA)

Perform principal component analysis on the expression data. The
function uses [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) to
perform PCA and
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) to
tidy the results. If `scale = TRUE`, constant variables (zero variance)
will be removed before PCA.

## Usage

``` r
gly_pca(
  exp,
  center = TRUE,
  scale = TRUE,
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

- log_transform:

  `NULL` (default) uses `metadata(exp)$scale` to avoid logging data
  twice. `TRUE` requests the method's log transformation for raw data;
  if the input is already on `"log"`, `"log2"`, `"log10"`, or `"logx"`,
  it warns and uses the supplied values unchanged. `FALSE` always uses
  the supplied values unchanged. This argument must be named when it
  follows `...`.

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

- `preprocessing`: Input and analysis scales and the applied
  transformation; see
  [gly_scale](https://glycoverse.github.io/glystats/reference/gly_scale.md).

## Details

For raw data with automatic transformation, the function performs log
transformation on the expression data (log(x + 1)) before PCA analysis.

## Required packages

This function only uses base R packages and does not require additional
dependencies.

## SE identifiers

For SE inputs, the result also contains an `identifiers` mapping of
result IDs to input positions and original names; see
[gly_identity](https://glycoverse.github.io/glystats/reference/gly_identity.md).
In SE tidy results, identity columns are `.variable`/`.sample` (dotted
endpoints for correlations); undotted metadata columns remain
annotations.

## See also

[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)

[gly_scale](https://glycoverse.github.io/glystats/reference/gly_scale.md)
for input scales, transformations, and fold changes.
[gly_identity](https://glycoverse.github.io/glystats/reference/gly_identity.md)
for SE identifiers, annotation alignment, and unnamed inputs.
