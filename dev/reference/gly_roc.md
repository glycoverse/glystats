# ROC Analysis for Glycomics and Glycoproteomics Data

Perform Receiver Operating Characteristic (ROC) analysis for binary
classification of glycomics or glycoproteomics data. The function
calculates ROC curves and Area Under the Curve (AUC) values for each
variable to assess their discriminatory power between two groups.

## Usage

``` r
gly_roc(exp, group_col = "group", pos_class = NULL, add_info = TRUE)

gly_roc_(expr_mat, groups, pos_class = NULL)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- group_col:

  (Only for `gly_roc()`) A character string specifying the column name
  of the grouping variable in the sample information. Default is
  `"group"`. The grouping variable must have exactly 2 levels for binary
  classification.

- pos_class:

  A character string specifying which group level should be treated as
  the positive class. If `NULL` (default), the second level
  (alphabetically) will be used as the positive class.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibbles. If FALSE, only the ROC
  analysis results are returned. Only applicable to `gly_roc()`.

- expr_mat:

  (Only for `gly_roc_()`) A numeric matrix with variables as rows and
  samples as columns.

- groups:

  (Only for `gly_roc_()`) A factor or character vector specifying group
  membership for each sample. Must have exactly 2 levels. Character
  vectors will be automatically converted to factors.

## Value

A list with three elements:

- `tidy_result`: A list containing two tibbles:

  - `auc`: A tibble containing AUC values for each variable with the
    following columns:

    - `variable`: Variable name

    - `auc`: Area Under the Curve value

    - `auc_ci_low`: Lower bound of 95% confidence interval for AUC

    - `auc_ci_high`: Upper bound of 95% confidence interval for AUC

  - `coords`: A tibble containing ROC curve coordinates with the
    following columns:

    - `variable`: Variable name

    - `threshold`: Threshold value for classification

    - `specificity`: Specificity (True Negative Rate)

    - `sensitivity`: Sensitivity (True Positive Rate)

- `raw_result`: A list of `pROC` objects

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_roc_res` and `glystats_res`.

## Details

For each variable, a ROC curve is computed using the expression values
as predictor and the binary group labels as response.

The function requires exactly 2 groups in the specified grouping
variable. If more than 2 groups are present, an error will be thrown.

`gly_roc()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_roc_()` is the underlying API that works with matrices and factor
vectors directly, providing more flexibility for users who don't use the
glyexp package.

**Underlying Function:**

- ROC analysis is performed using
  [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html)

- Coordinates are extracted using
  [`pROC::coords()`](https://rdrr.io/pkg/pROC/man/coords.html)

## Required packages

This function requires the `pROC` package to be installed for ROC curve
computation.

## See also

[`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html),
[`pROC::coords()`](https://rdrr.io/pkg/pROC/man/coords.html)
