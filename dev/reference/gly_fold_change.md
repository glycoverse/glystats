# Calculate fold change

Calculate fold change for a given expression matrix and group
information. Pairwise comparisons are returned when more than two groups
are supplied. When you run this function, you will see message about
"Ref Group" and "Test Group". "Ref Group" is the reference group, and
"Test Group" is the test/treatment/case group.

## Usage

``` r
gly_fold_change(
  exp,
  group_col = "group",
  add_info = TRUE,
  log_transform = NULL
)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment`.

- group_col:

  The column name of the group information in the sample information.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the fold
  change results are returned.

- log_transform:

  `NULL` (default) uses `metadata(exp)$scale` to avoid logging data
  twice. `TRUE` forces the method's log transformation on the supplied
  values, even if they are already logged; `FALSE` uses the supplied
  values unchanged. This argument must be named when it follows `...`.

## Value

A list with four elements:

- `tidy_result`: A tibble with fold change results containing the
  following columns:

  - `variable`: Variable name

  - `log2fc`: Log2 fold change; see
    [gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md)
    for the scale-dependent definition If more than two groups, two
    additional columns `ref_group` and `test_group` will be added.

- `raw_result`: The raw result (same as tidy_result for this function)

- `meta_data`: A list containing metadata from the input experiment

- `preprocessing`: Input and analysis scales and the applied
  transformation; see
  [gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md).

## See also

[gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md)
for input scales, transformations, and fold changes.
