# Cox Proportional Hazards Model for Survival Analysis

Fit a Cox proportional hazards model for each variable in the expression
data, and extract p-values and hazard ratios from it.

## Usage

``` r
gly_cox(
  exp,
  time_col = "time",
  event_col = "event",
  p_adj_method = "BH",
  add_info = TRUE,
  ...
)

gly_cox_(expr_mat, time, event, p_adj_method = "BH", ...)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- time_col:

  A character string specifying the column name in sample information
  that contains survival time. Default is "time".

- event_col:

  A character string specifying the column name in sample information
  that contains event indicator (1 for event, 0 for censoring). Default
  is "event".

- p_adj_method:

  A character string specifying the method to adjust p-values. See
  `p.adjust.methods` for available methods. Default is "BH". If NULL, no
  adjustment is performed.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the Cox
  model results are returned. Only applicable to `gly_cox()`.

- ...:

  Additional arguments passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

- expr_mat:

  (Only for `gly_cox_()`) A numeric matrix with variables as rows and
  samples as columns.

- time:

  A numeric vector specifying the survival time for each sample.

- event:

  A numeric vector specifying the event indicator for each sample (1 for
  event, 0 for censoring).

## Value

A list with three elements:

- `tidy_result`: A tibble with Cox model results containing the
  following columns:

  - `variable`: Variable name

  - `coefficient`: Regression coefficient (log hazard ratio)

  - `std.error`: Standard error of the coefficient

  - `statistic`: Wald test statistic

  - `p_val`: Raw p-value from Wald test

  - `hr`: Hazard ratio (exp(coefficient))

  - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

- `raw_result`: A list of raw `coxph` model objects.

- `meta_data`: A list containing metadata from the input experiment

## Details

The function fits an Cox proportional hazards model for each variable
by:

    coxph(Surv(time, event) ~ expression_value)

P-values are adjusted by Benjamini-Hochberg method by default.

## See also

[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html),
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
