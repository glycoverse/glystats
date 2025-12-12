# Partial Least Squares Discriminant Analysis (PLS-DA)

Perform partial least squares discriminant analysis on the expression
data. The function uses
[`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html) to perform
PLS-DA and returns tidy results.

## Usage

``` r
gly_plsda(
  exp,
  group_col = "group",
  ncomp = 2,
  scale = TRUE,
  add_info = TRUE,
  ...
)

gly_plsda_(expr_mat, groups, ncomp = 2, scale = TRUE, ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- group_col:

  A character string specifying the column name in sample information
  that contains group labels. Default is "group".

- ncomp:

  An integer indicating the number of components to include. Default is
  2.

- scale:

  A logical indicating whether to scale the data. Default is TRUE.

- add_info:

  A logical value. If TRUE (default), sample and variable information
  from the experiment will be added to the result tibbles. If FALSE,
  only the PLS-DA results are returned. Only applicable to
  `gly_plsda()`.

- ...:

  Additional arguments passed to
  [`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html).

- expr_mat:

  A numeric matrix with variables as rows and samples as columns.

- groups:

  A factor or character vector specifying group membership for each
  sample. Character vectors will be automatically converted to factors.

## Value

A list containing:

- `tidy_result`: A list of tibbles with PLS-DA results:

  - `samples`: PLS-DA scores for each sample containing the following
    columns:

    - `sample`: Sample name

    - `group`: Group assignment

    - `p1`, `p2`, etc.: PLS-DA component scores

  - `variables`: PLS-DA loadings for each variable containing the
    following columns:

    - `variable`: Variable name

    - `p1`, `p2`, etc.: PLS-DA component loadings

    - `pcorr1`, `pcorr2`, etc.: Correlation between each variable and
      component

  - `variance`: PLS-DA explained variance containing the following
    columns:

    - `component`: Component name (p1, p2, etc.)

    - `prop_var_explained`: Proportion of variance explained by each
      component

    - `cumulative_prop_var`: Cumulative proportion of variance explained

  - `vip`: Variable Importance in Projection scores containing the
    following columns:

    - `variable`: Variable name

    - `vip`: VIP score

  - `perm_test`: Permutation test results containing the following
    columns:

    - `model`: Model type ("Original" for the original model,
      "Permutation" for permuted models)

    - `perm_id`: Permutation ID (0 for original model, 1+ for
      permutations)

    - Additional columns from the permutation test matrix (e.g., R2X,
      R2Y, Q2, etc.)

- `raw_result`: The raw ropls opls object from
  [`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html)

## Details

`gly_plsda()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_plsda_()` is the underlying API that works with matrices and factor
vectors directly, providing more flexibility for users who don't use the
glyexp package.

## Required packages

This function requires the following packages to be installed:

- `ropls` for PLS-DA analysis

## See also

[`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html)
