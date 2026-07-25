# Formula-based Linear Models

Fit a linear model to every variable in a glycomics or glycoproteomics
experiment using limma's empirical Bayes moderation.

## Usage

``` r
gly_linear_model(
  exp,
  formula,
  contrasts = NULL,
  p_adj_method = "BH",
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
  object, or another `SummarizedExperiment` containing expression data
  and sample information.

- formula:

  A one-sided formula describing the model in terms of columns in the
  sample information, for example `~ treatment * time + batch`.

- contrasts:

  An optional uniquely named character vector of limma-style contrast
  expressions. Use backticks around non-syntactic coefficient names,
  such as `` `treatmentB:timeT2` ``.

- p_adj_method:

  A character string specifying the method used to adjust p-values
  within each coefficient or contrast. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  available methods. If `NULL`, no adjustment is performed.

- add_info:

  A logical value. If `TRUE` (default), variable information from the
  experiment is added to the tidy result.

- ...:

  Additional arguments passed to
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html).

## Value

A list with classes `glystats_linear_model_res` and `glystats_res`
containing:

- `tidy_result`: A tibble with `variable`, `term`, `estimate`,
  `AveExpr`, `t`, `p_val`, optional `p_adj`, and `b` columns.

- `raw_result`: A list containing the moderated limma fit, design
  matrix, contrast matrix, and coefficient-name mapping.

- `meta_data`: Metadata copied from `exp`.

## Details

`gly_linear_model()` and
[`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
both use limma under the hood, but serve complementary use cases.
`gly_linear_model()` provides a general formula interface for specifying
a wide range of analysis designs and custom contrasts.
[`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
is a dedicated wrapper for common differential expression analysis
tasks, with convenient group, covariate, subject, and
pairwise-comparison arguments.

Expression values are transformed using `log2(x + 1e-6)` before
modeling. When `contrasts` is `NULL`, all non-intercept coefficients are
tested. When contrasts are supplied, only the named contrasts are
returned.

Character predictors are handled as factors by
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
Factor reference levels therefore determine the meaning of coefficients.
Set factor levels in the sample information before calling this function
when a specific reference level is required.

P-values are adjusted independently across variables for each
coefficient or contrast.

## See also

[`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md),
[`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html),
[`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html)

## Examples

``` r
exp <- glyexp::real_experiment
complete <- rowSums(is.finite(SummarizedExperiment::assay(exp))) == ncol(exp)
exp <- exp[which(complete)[seq_len(10)], ]
SummarizedExperiment::colData(exp)$time <- factor(
  rep(c("early", "late"), length.out = ncol(exp))
)

model_res <- gly_linear_model(exp, ~ group * time)

contrast_res <- gly_linear_model(
  exp,
  ~ group * time,
  contrasts = c(
    C_late_vs_H_late = "groupC + `groupC:timelate`"
  )
)
```
