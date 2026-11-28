# Match variables to a group template

Assign each sample the template score for its experimental group and
correlate each variable's abundance with those scores. Positive and
negative coefficients identify matches in opposite directions.

## Usage

``` r
gly_pattern_hunter(
  exp,
  template,
  group_col = "group",
  method = "pearson",
  p_adj_method = "BH",
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
  object, or another `SummarizedExperiment` containing an expression
  matrix and sample information.

- template:

  A named numeric vector with one finite score for each observed group
  in `group_col`. Names must match the group labels exactly. Groups may
  share a score. At least two distinct scores are required.

- group_col:

  Name of the sample metadata column containing the group labels. Each
  sample must have a non-missing, non-empty label.

- method:

  Correlation coefficient, one of `"pearson"` (default), `"spearman"`,
  or `"kendall"`.

- p_adj_method:

  A character string specifying the method to adjust p-values. See
  `p.adjust.methods` for available methods. Default is "BH". If NULL, no
  adjustment is performed.

- add_info:

  Whether to include variable annotations in the tidy result.

- log_transform:

  `NULL` (default) uses `metadata(exp)$scale` to avoid logging data
  twice. `TRUE` requests the method's log transformation for raw data;
  if the input is already on `"log"`, `"log2"`, `"log10"`, or `"logx"`,
  it warns and uses the supplied values unchanged. `FALSE` always uses
  the supplied values unchanged. This argument must be named when it
  follows `...`.

## Value

A `glystats_pattern_hunter_res` object inheriting from `glystats_res`.
[`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
returns one row per input variable in assay order.

- `.variable` gives the SE variable ID (`variable` for legacy inputs).

- `cor` is the signed correlation coefficient.

- `n_samples`, `n_groups`, and `n_scores` count the finite samples,
  represented groups, and distinct template scores for the variable.

- `p_val` is the two-sided p-value and `p_adj` is its adjusted value.
  Setting `p_adj_method = NULL` omits `p_adj`.

[`get_raw_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
returns a list of
[`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html) results
named by variable ID, with `NULL` for untestable variables. The result
records `template`, `group_col`, `method`, `p_adj_method`, the input
`meta_data`, and the `preprocessing` used. SE inputs also retain the
`identifiers` mapping described in
[gly_identity](https://glycoverse.github.io/glystats/dev/reference/gly_identity.md).

## Details

Each finite sample observation contributes once. Groups with more
observed samples have more weight. For each variable, the test uses the
groups and scores represented among its finite abundances on the
analysis scale. If a group is entirely missing, the coefficient
describes the remaining groups. A test requires at least three such
samples, two distinct template scores, and variation in abundance. Other
rows retain their counts and receive `NA` for the coefficient and
p-values.

Raw abundances use `log2(x + 1e-6)` by default. See
[gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md)
for scale overrides. With relative abundances, a match describes changes
only relative to the normalization reference, such as the total measured
abundance. Even a strong positive correlation with an increasing
template can include a decrease between adjacent groups.

The two-sided tests assume independent samples and a prespecified
template. For fixed scores, the Pearson test is the slope test in a
linear regression with independent normal errors of constant variance.
The rank tests assume exchangeability under independence. Spearman uses
average ranks for ties and Kendall uses tau-b. With ties,
[`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html) uses
approximate rank-test p-values, which can be inaccurate with few
samples. Without ties, it uses its usual sample-size-dependent
calculation.

Multiplicity adjustment covers the testable variables for this template.
If several templates are screened, account for the full set of tests. A
template chosen after inspecting the data needs independent validation.

## References

Pavlidis P, Noble WS (2001). Analysis of strain and regional variation
in gene expression in mouse brain. *Genome Biology*, 2, research0042.
[doi:10.1186/gb-2001-2-10-research0042](https://doi.org/10.1186/gb-2001-2-10-research0042)
.

Pavlidis P (2003). Using ANOVA for gene selection from microarray
studies of the nervous system. *Methods*, 31, 282-289.
[doi:10.1016/S1046-2023(03)00157-9](https://doi.org/10.1016/S1046-2023%2803%2900157-9)
.

For a related pattern search in metabolomics, see Xia J, Mandal R,
Sinelnikov IV, Broadhurst D, Wishart DS (2012). MetaboAnalyst 2.0—a
comprehensive server for metabolomic data analysis. *Nucleic Acids
Research*, 40, W127-W133.
[doi:10.1093/nar/gks374](https://doi.org/10.1093/nar/gks374) .

## Examples

``` r
exp <- glyexp::real_experiment2
groups <- SummarizedExperiment::colData(exp)$group
exp <- exp[1:5, groups %in% c("H", "C")]
res <- gly_pattern_hunter(exp, template = c(H = 0, C = 1), add_info = FALSE)
get_tidy_result(res)
#> # A tibble: 5 × 7
#>   .variable                     cor n_samples n_groups n_scores p_val p_adj
#>   <chr>                       <dbl>     <int>    <int>    <int> <dbl> <dbl>
#> 1 Hex(3)HexNAc(3)         -0.0745          69        2        2 0.543 0.995
#> 2 Hex(3)HexNAc(7)         -0.0283          50        2        2 0.845 0.995
#> 3 Hex(5)HexNAc(2)          0.153           73        2        2 0.196 0.982
#> 4 Hex(6)HexNAc(4)NeuAc(2)  0.000808        69        2        2 0.995 0.995
#> 5 Hex(4)HexNAc(3)         -0.0261          72        2        2 0.828 0.995
```
