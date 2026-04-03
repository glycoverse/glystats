# Correlation Analysis for Glycomics and Glycoproteomics Data

Perform pairwise correlation analysis on variables or samples in the
expression data. The function calculates correlation coefficients and
p-values for all pairs, with optional multiple testing correction.

## Usage

``` r
gly_cor(exp, on = "variable", method = "pearson", p_adj_method = "BH", ...)

gly_cor_(
  expr_mat,
  on = "variable",
  method = "pearson",
  p_adj_method = "BH",
  ...
)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object containing expression matrix and sample information.

- on:

  A character string specifying what to correlate. Either "variable"
  (default) to correlate variables/features, or "sample" to correlate
  samples/observations.

- method:

  A character string indicating which correlation coefficient is to be
  computed. One of "pearson" (default) or "spearman". Note: "kendall" is
  not supported by Hmisc::rcorr.

- p_adj_method:

  A character string specifying the method to adjust p-values. See
  `p.adjust.methods` for available methods. Default is "BH". If NULL, no
  adjustment is performed.

- ...:

  Additional arguments passed to
  [`Hmisc::rcorr()`](https://rdrr.io/pkg/Hmisc/man/rcorr.html).

- expr_mat:

  (Only for `gly_cor_()`) A numeric matrix with variables as rows and
  samples as columns.

## Value

A list with three elements:

- `tidy_result`: A tibble with correlation results containing the
  following columns:

  - `variable1` or `sample1`: First element of the pair (depending on
    `on` parameter)

  - `variable2` or `sample2`: Second element of the pair (depending on
    `on` parameter)

  - `cor`: Correlation coefficient

  - `p_val`: Raw p-value from correlation test

  - `p_adj`: Adjusted p-value (if p_adj_method is not NULL)

- `raw_result`: The raw rcorr object from Hmisc::rcorr()

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_cor_res` and `glystats_res`.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1)) before correlation analysis. When `on = "variable"`
(default), correlations are calculated between variables across samples.
When `on = "sample"`, correlations are calculated between samples across
variables.

`gly_cor()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects。

`gly_cor_()` is the underlying API that works with matrices directly,
providing more flexibility for users who don't use the glyexp package.

**Correlation Calculation:** Correlation coefficients and p-values are
calculated using
[`Hmisc::rcorr()`](https://rdrr.io/pkg/Hmisc/man/rcorr.html) which is
more efficient than pairwise
[`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html) calls.

**Multiple Testing Correction:** P-values are adjusted for multiple
testing using the method specified by `p_adj_method`.

## Required packages

This function requires the `Hmisc` package for efficient correlation
calculation.

## See also

[`Hmisc::rcorr()`](https://rdrr.io/pkg/Hmisc/man/rcorr.html),
[`stats::cor()`](https://rdrr.io/r/stats/cor.html)
