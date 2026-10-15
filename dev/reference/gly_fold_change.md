# Calculate fold change

Calculate fold change for a given expression matrix and group
information. It could only be used for 2-group analysis. When you run
this function, you will see message about "Ref Group" and "Test Group".
"Ref Group" is the reference group, and "Test Group" is the
test/treatment/case group.

## Usage

``` r
gly_fold_change(
  exp,
  group_col = "group",
  add_info = TRUE,
  log_transform = FALSE
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

  A single non-missing logical. If `FALSE` (default), return the log2
  ratio of arithmetic group means as before. If `TRUE`, apply
  `log2(x + 1e-6)` and return the difference of group means on that
  scale, equivalent to the log2 ratio of geometric means of `x + 1e-6`.
  This option transforms input expression; it does not toggle the log2
  scale of the reported fold change.

## Value

A list with three elements:

- `tidy_result`: A tibble with fold change results containing the
  following columns:

  - `variable`: Variable name

  - `log2fc`: Log2 ratio of arithmetic group means by default, or
    difference of mean log2-transformed values when
    `log_transform = TRUE`. If more than two groups, two additional
    columns `ref_group` and `test_group` will be added.

- `raw_result`: The raw result (same as tidy_result for this function)

- `meta_data`: A list containing metadata from the input experiment
