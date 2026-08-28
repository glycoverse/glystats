# Calculate fold change

Calculate fold change for a given expression matrix and group
information. It could only be used for 2-group analysis. When you run
this function, you will see message about "Ref Group" and "Test Group".
"Ref Group" is the reference group, and "Test Group" is the
test/treatment/case group.

## Usage

``` r
gly_fold_change(exp, group_col = "group", add_info = TRUE)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://rdrr.io/pkg/glyexp/man/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://rdrr.io/pkg/glyexp/man/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment`.

- group_col:

  The column name of the group information in the sample information.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the fold
  change results are returned.

## Value

A list with three elements:

- `tidy_result`: A tibble with fold change results containing the
  following columns:

  - `variable`: Variable name

  - `log2fc`: Log2 fold change (log2(group2_mean / group1_mean)) If more
    than two groups, two additional columns `ref_group` and `test_group`
    will be added.

- `raw_result`: The raw result (same as tidy_result for this function)

- `meta_data`: A list containing metadata from the input experiment
