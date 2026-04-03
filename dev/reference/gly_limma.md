# Linear Models for Differential Expression Analysis using limma

Perform differential expression analysis using linear models with
empirical Bayes moderation from the limma package. Supports both
two-group and multi-group comparisons.

## Usage

``` r
gly_limma(
  exp,
  group_col = "group",
  covariate_cols = NULL,
  subject_col = NULL,
  p_adj_method = "BH",
  ref_group = NULL,
  contrasts = NULL,
  add_info = TRUE,
  ...
)

gly_limma_(
  expr_mat,
  groups,
  p_adj_method = "BH",
  ref_group = NULL,
  contrasts = NULL,
  covariates = NULL,
  subjects = NULL,
  ...
)
```

## Arguments

- exp:

  A `glyexp_experiment` object containing expression data and sample
  information.

- group_col:

  (Only for `gly_limma()`) A character string specifying the column name
  in sample information that contains group labels. Default is "group".

- covariate_cols:

  (Only for `gly_limma()`) A character vector specifying column names in
  sample information to include as covariates in the limma model.
  Default is NULL.

- subject_col:

  (Only for `gly_limma()`) A character string specifying the column name
  in sample information that contains subject identifiers for paired
  comparisons. Default is NULL.

- p_adj_method:

  A character string specifying the method for multiple testing
  correction. Must be one of the methods supported by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  is "BH" (Benjamini-Hochberg). Set to NULL to skip p-value adjustment.

- ref_group:

  A character string specifying the reference group. If NULL (default),
  the first level of the group factor is used as the reference. Only
  used for two-group comparisons.

- contrasts:

  A character vector specifying custom contrasts. If NULL (default), all
  pairwise comparisons are automatically generated, and the levels
  coming first in the factor will be used as the reference group.
  Supports two formats: "group1-group2" or "group1_vs_group2". Use the
  second format if group names contain hyphens. "group1" will be used as
  the reference group.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned. Only applicable to `gly_limma()`.

- ...:

  Additional arguments passed to
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html).

- expr_mat:

  (Only for `gly_limma_()`) A numeric matrix with variables as rows and
  samples as columns.

- groups:

  (Only for `gly_limma_()`) A factor or character vector specifying
  group membership for each sample. Must have at least 2 levels.
  Character vectors will be automatically converted to factors. If
  `contrasts` is not provided, the levels coming first in the factor
  will be used as the reference group.

- covariates:

  (Only for `gly_limma_()`) A data frame, matrix, or vector of
  sample-level covariates. Must have the same number of rows as
  `expr_mat` has columns. If row names are provided and match
  `colnames(expr_mat)`, they will be aligned automatically. Default is
  NULL.

- subjects:

  (Only for `gly_limma_()`) A vector or factor of subject identifiers
  for paired comparisons. Must have length equal to `ncol(expr_mat)`. If
  names or row names are provided and match `colnames(expr_mat)`, they
  will be aligned automatically. Default is NULL.

## Value

A list with three elements:

- `tidy_result`: A tibble with limma results containing the following
  columns:

  - `variable`: Variable name

  - `log2fc`: Log2 fold change

  - `AveExpr`: Average expression level

  - `t`: t-statistic

  - `p_val`: Raw p-value

  - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

  - `b`: B-statistic (log-odds of differential expression)

  - `ref_group`: Reference group

  - `test_group`: Test group

- `raw_result`: The raw limma fit object(s).

- `meta_data`: The meta data from the experiment object (only for
  `gly_limma()`).

## Details

The function performs log2 transformation on the expression data
(log2(x + 1)) before statistical testing. The analysis uses linear
models with empirical Bayes moderation to improve statistical power,
especially for small sample sizes.

`gly_limma()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_limma_()` is the underlying API that works with matrices and factor
vectors directly, providing more flexibility for users who don't use the
glyexp package.

For two or more groups, all pairwise comparisons are automatically
generated (one contrast for two groups) and performed using contrast
matrices unless custom contrasts are specified.

If covariates are provided, they are added as additional terms in the
design matrix and are not part of the group contrasts.

If subjects are provided, they are included as a blocking factor in the
design matrix using the formula `~ 0 + group + subject` (plus any
covariates).

When specifying custom contrasts, use either "group1-group2" or
"group1_vs_group2" format. If group names contain hyphens and you use
the first format, an error will be raised suggesting to use the second
format.

## Required packages

This function requires the following packages to be installed:

- `limma` for linear model fitting and empirical Bayes moderation

## See also

[`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html),
[`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html),
[`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html)
