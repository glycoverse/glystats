# Test for differences between two independent groups using SAM

Apply Significance Analysis of Microarrays (SAM) to glycan or
glycopeptide abundances. Positive scores indicate higher mean abundance
in the test group on the analysis scale.

## Usage

``` r
gly_sam(
  exp,
  group_col = "group",
  ref_group = NULL,
  delta = NULL,
  n_perm = 100L,
  s0 = NULL,
  seed = 1L,
  add_info = TRUE,
  log_transform = NULL
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

- group_col:

  A character string specifying the column name in sample information
  that contains group labels. Default is "group".

- ref_group:

  A character string specifying the reference group. If NULL (default),
  the first level of the group factor is used as the reference. Only
  used for two-group comparisons.

- delta:

  A positive threshold or a vector of distinct positive thresholds for
  the departure of ordered scores from their permutation expectations.
  The default, `NULL`, returns the ten grid entries chosen by
  [`siggenes::samControl()`](https://rdrr.io/pkg/siggenes/man/samControl.html).
  Supply one value to call individual features. With a grid, `called` is
  `NA` for every feature.

- n_perm:

  Requested number of label permutations, at least two. Uses all
  allocations if fewer exist. For more than 30,000 possible allocations,
  siggenes can draw the same allocation more than once.

- s0:

  Nonnegative regularization constant in analysis units. By default
  siggenes estimates it across features. This requires at least 25
  distinct testable standard errors and variation in the scores used to
  estimate it. A prespecified value also permits analysis of smaller
  panels with at least two testable features.

- seed:

  Nonnegative integer seed for the permutations. The caller's
  random-number state is restored after fitting.

- add_info:

  A logical value. If TRUE (default), variable information from the
  experiment will be added to the result tibble. If FALSE, only the
  statistical results are returned.

- log_transform:

  `NULL` (default) uses `metadata(exp)$scale` to avoid logging data
  twice. `TRUE` requests the method's log transformation for raw data;
  if the input is already on `"log"`, `"log2"`, `"log10"`, or `"logx"`,
  it warns and uses the supplied values unchanged. `FALSE` always uses
  the supplied values unchanged. This argument must be named when it
  follows `...`.

## Value

A `glystats_sam_res` object inheriting from `glystats_res`.
[`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
returns one row per input feature in assay order.

- `.variable` is the SE feature ID (`variable` for legacy experiments).
  `ref_group` and `test_group` identify the comparison.

- `d` is the signed SAM score. `p_val` is the backend's two-sided
  p-value from the pooled permuted scores. `q_value` is its feature
  q-value computed from those p-values and the estimated null proportion
  (`q.version = 1`).

- `called` uses the inclusive lower and upper score cutoffs for a single
  supplied delta. `n_finite` counts finite observations on the analysis
  scale. Excluded rows have missing statistics and calls, and an
  `exclusion` of `"non_finite"`, `"backend_missing_code"`, or
  `"unestimable_variance"`.

`thresholds` has one row per delta, with the lower and upper score
cutoffs, number called (`n_called`), mean number of false calls under
permutation (`expected_false`), estimated null proportion (`p0`), and
the backend's permutation FDR (`fdr`, a fraction). siggenes reports an
FDR of zero when a threshold calls zero features.

`diagnostics` records the input and tested feature counts, fitted `s0`,
and the numbers of actual and distinct permutations. `settings` records
the requested comparison and tuning parameters.
[`get_raw_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
returns the fitted `SAM` object. Its feature labels are the original
one-based row positions so even long SE identifiers map unambiguously.
The result retains the input `meta_data`, the applied `preprocessing`
and the SE `identifiers` mapping described in
[gly_identity](https://glycoverse.github.io/glystats/dev/reference/gly_identity.md).

## Details

Supply two groups with at least two independent samples each. The null
distribution assumes sample labels are exchangeable and preserves group
sizes. Paired, blocked, covariate-adjusted and multigroup designs are
unsupported. For blocking or covariate adjustment, see
[`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md).

Install the optional Bioconductor backend with
`BiocManager::install("siggenes")`. Its dependency `multtest` supplies
the observed score components used to check which features are testable.

[`siggenes::sam()`](https://rdrr.io/pkg/siggenes/man/sam.html) uses the
test-minus-reference mean difference divided by the unpooled standard
error plus `s0` (its two-class Welch statistic). Raw assays use
`log2(x + 1e-6)` by default. For relative abundances, the sign describes
a change relative to the normalization reference. With `logx`, it
describes a difference in the supplied coordinates.

Only features with finite values in every sample and estimable
within-group variation enter the common permutation distribution. The
backend treats numerically zero standard errors as untestable, including
features constant within both groups. It also reserves the finite value
`-93074815` as a missing-data code on the analysis scale. A feature
containing that value is excluded before permutation. Filtering is
reported, and each excluded feature retains its ID and reason. The
p-values, q-values and threshold FDR estimates apply to the retained
family. No imputation or fold-change filter is applied.

## References

Tusher VG, Tibshirani R, Chu G (2001). Significance analysis of
microarrays applied to the ionizing radiation response. *Proceedings of
the National Academy of Sciences*, 98, 5116-5121.
[doi:10.1073/pnas.091062498](https://doi.org/10.1073/pnas.091062498) .

[siggenes on Bioconductor](https://bioconductor.org/packages/siggenes).

## Examples

``` r
exp <- glyexp::real_experiment2
group <- SummarizedExperiment::colData(exp)$group
two_groups <- exp[, group %in% c("H", "C")]
res <- gly_sam(two_groups, ref_group = "H", delta = 1, seed = 42)
#> Warning: SAM inference uses 29 of 66 features.
#> ℹ All input IDs are retained. See exclusion for the reason each feature was
#>   omitted.
res$thresholds
#> # A tibble: 1 × 7
#>   delta cutoff_lower cutoff_upper n_called expected_false    p0   fdr
#>   <dbl>        <dbl>        <dbl>    <int>          <dbl> <dbl> <dbl>
#> 1     1         -Inf         1.31       12           2.43     1 0.202
head(subset(get_tidy_result(res), called, select = c(.variable, d, q_value)))
#> # A tibble: 6 × 3
#>   .variable                          d q_value
#>   <chr>                          <dbl>   <dbl>
#> 1 Hex(5)HexNAc(2)                 1.39   0.425
#> 2 Hex(3)HexNAc(4)dHex(1)          6.03   0    
#> 3 Hex(5)HexNAc(3)NeuAc(1)         1.74   0.376
#> 4 Hex(5)HexNAc(5)dHex(1)          1.46   0.408
#> 5 Hex(5)HexNAc(5)dHex(1)NeuAc(1)  1.31   0.458
#> 6 Hex(3)HexNAc(5)dHex(1)          4.47   0    
```
