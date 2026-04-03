# Calculate fold change

Calculate fold change for a given expression matrix and group
information. It could only be used for 2-group analysis. When you run
this function, you will see message about "Ref Group" and "Test Group".
"Ref Group" is the reference group, and "Test Group" is the
test/treatment/case group.

## Usage

``` r
gly_fold_change(exp, group_col = "group", add_info = TRUE)

gly_fold_change_(expr_mat, groups)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- group_col:

  (Only for `gly_fold_change()`) The column name of the group
  information in the sample information.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the fold
  change results are returned. Only applicable to `gly_fold_change()`.

- expr_mat:

  (Only for `gly_fold_change_()`) A numeric matrix with variables as
  rows and samples as columns.

- groups:

  (Only for `gly_fold_change_()`) A factor or character vector
  specifying group membership for each sample. Character vectors will be
  automatically converted to factors. If two groups, the first level is
  the reference group. If more than two groups, pairwise comparisons
  will be performed, with levels coming first as reference groups.

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

## Details

`gly_fold_change()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_fold_change_()` is the underlying API that works with matrices and
factor vectors directly, providing more flexibility for users who don't
use the glyexp package.
