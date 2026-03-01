# Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)

Perform orthogonal partial least squares discriminant analysis on the
expression data. The function uses
[`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html) to perform
OPLS-DA and returns tidy results.

## Usage

``` r
gly_oplsda(
  exp,
  group_col = "group",
  pred_i = 1,
  ortho_i = NA,
  scale = TRUE,
  add_info = TRUE,
  ...
)

gly_oplsda_(expr_mat, groups, pred_i = 1, ortho_i = NA, scale = TRUE, ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- group_col:

  (Only for `gly_oplsda()`) A character string specifying the column
  name in sample information that contains group labels. Default is
  "group".

- pred_i:

  An integer indicating the number of predictive components to include.
  Default is 1.

- ortho_i:

  An integer indicating the number of orthogonal components to include.
  Default is NA (automatic).

- scale:

  A logical indicating whether to scale the data. Default is TRUE.

- add_info:

  A logical value. If TRUE (default), sample and variable information
  from the experiment will be added to the result tibbles. If FALSE,
  only the OPLS-DA results are returned. Only applicable to
  `gly_oplsda()`.

- ...:

  Additional arguments passed to
  [`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html).

- expr_mat:

  (Only for `gly_oplsda_()`) A numeric matrix with variables as rows and
  samples as columns.

- groups:

  (Only for `gly_oplsda_()`) A factor or character vector specifying
  group membership for each sample. Must have exactly 2 levels.
  Character vectors will be automatically converted to factors.

## Value

A list containing:

- `tidy_result`: A list of tibbles with OPLS-DA results:

  - `samples`: OPLS-DA scores for each sample containing the following
    columns:

    - `sample`: Sample name

    - `group`: Group assignment

    - `p1`, `p2`, etc.: Predictive component scores

    - `o1`, `o2`, etc.: Orthogonal component scores

  - `variables`: OPLS-DA loadings for each variable containing the
    following columns:

    - `variable`: Variable name

    - `p1`, `p2`, etc.: Predictive component loadings

    - `o1`, `o2`, etc.: Orthogonal component loadings

    - `pcorr1`, `pcorr2`, etc.: Correlation between each variable and
      the corresponding predictive component

  - `variance`: OPLS-DA explained variance containing the following
    columns:

    - `component`: Component name (p1, o1, etc.)

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

`gly_oplsda()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_oplsda_()` is the underlying API that works with matrices and
factor vectors directly, providing more flexibility for users who don't
use the glyexp package.

## Required packages

This function requires the following packages to be installed:

- `ropls` for OPLS-DA analysis

## See also

[`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html)
