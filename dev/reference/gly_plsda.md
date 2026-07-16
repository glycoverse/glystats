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
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object, or another `SummarizedExperiment` containing an expression
  matrix and sample information.

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
  only the PLS-DA results are returned.

- ...:

  Additional arguments passed to
  [`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html).

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

- `meta_data`: A list containing metadata from the input experiment

## Required packages

This function requires the following packages to be installed:

- `ropls` for PLS-DA analysis

## See also

[`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html)
