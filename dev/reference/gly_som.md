# Map abundance profiles with a self-organizing map

Train an unsupervised self-organizing map (SOM) and assign each profile
to its closest grid unit. Training favors similar prototypes at
neighboring units.

## Usage

``` r
gly_som(
  exp,
  on = "variable",
  grid = NULL,
  topology = "rectangular",
  rlen = 100L,
  scale = TRUE,
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
  object, or another `SummarizedExperiment` containing an expression
  matrix and sample information.

- on:

  Axis to map. `"variable"` maps each glycan or glycoform across the
  supplied samples. `"sample"` maps each sample across variables.

- grid:

  Grid dimensions as two positive integers, `c(xdim, ydim)`. By default,
  the map is square with side length `min(5, floor(sqrt(n_train)))`,
  where `n_train` counts usable profiles.

- topology:

  Arrangement of the grid units, either `"rectangular"` (default) or
  `"hexagonal"`.

- rlen:

  Number of training epochs. Each epoch makes as many online updates as
  there are training profiles.

- scale:

  Whether to center and standardize each input coordinate across the
  retained profiles, as in
  [`gly_kmeans()`](https://glycoverse.github.io/glystats/dev/reference/gly_kmeans.md).
  With at least two training profiles, each sample (variable map) or
  variable (sample map) is centered at its mean and divided by its
  sample standard deviation. Constant coordinates are centered at zero.
  A single profile is used directly.

- seed:

  Nonnegative integer seed for initialization, online training and tie
  breaking. The default restores the caller's random-number state after
  fitting, as does any supplied integer. Use `NULL` to draw from the
  caller's random-number stream.

- add_info:

  Whether to attach annotations for the mapped axis.

- log_transform:

  `NULL` (default) uses `metadata(exp)$scale` to avoid logging data
  twice. `TRUE` requests the method's log transformation for raw data;
  if the input is already on `"log"`, `"log2"`, `"log10"`, or `"logx"`,
  it warns and uses the supplied values unchanged. `FALSE` always uses
  the supplied values unchanged. This argument must be named when it
  follows `...`.

## Value

A `glystats_som_res` object inheriting from `glystats_res`.
[`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
returns a tibble in input order with one row per mapped-axis ID
(`.variable` or `.sample` for SE inputs).

- `unit`, `x`, and `y` identify the winning unit and its grid
  coordinates.

- `distance` is the Euclidean distance to that unit's fitted prototype
  in the training coordinates. Smaller values indicate a closer fit.

- `n_finite` counts finite coordinates. Excluded profiles retain their
  IDs and receive missing assignments and distances. `exclusion` is
  `"no_finite_values"` or `"incomplete_profile"` for these profiles and
  `NA` for the training profiles.

`diagnostics` contains `n_train`, `units` (occupancy of every grid unit,
including empty units), and `quantization_error` (mean training
distance).
[`get_raw_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
returns the
[`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html) fit
with its training data, fitted assignments, codebook vectors, and
per-epoch changes. With zero usable profiles, the automatic grid is
empty, quantization error is missing, and the fit is `NULL`. An explicit
grid requires at least one usable profile per requested unit.

`settings` records the axis, grid dimensions, topology, epochs,
standardization request and seed. `preprocessing` records the abundance
transformation, whether standardization was applied (`scaled`), and the
IDs of wholly non-finite coordinates (`excluded_coordinates`). The
result retains the input `meta_data` and SE
[identifiers](https://glycoverse.github.io/glystats/dev/reference/gly_identity.md).

## Details

The columns of the input assay represent individual samples. To map
group-average curves, first supply an experiment with one column per
group mean. Each supplied coordinate contributes to Euclidean distance.
Standardization gives each varying coordinate equal variance. With
`scale = FALSE`, coordinates with larger spreads can dominate the map.
With relative abundances, a change in one glycan's share can affect the
shares of other glycans. Distances describe similarity relative to the
normalization reference and depend on the abundance transformation.

Raw assays use `log2(x + 1e-6)` by default. See
[gly_scale](https://glycoverse.github.io/glystats/dev/reference/gly_scale.md)
for other scales and overrides. On the analysis scale, coordinates with
no finite values are removed first. Training uses profiles with finite
values at every remaining coordinate. A complete flat profile is
eligible because its level contributes to distance. Initial grid
prototypes are sampled without replacement from the training profiles.
Changing the grid or seed can change the assignments, and unit numbers
are local to each fit.

## References

Kohonen T (1982). Self-organized formation of topologically correct
feature maps. *Biological Cybernetics*, 43, 59-69.
[doi:10.1007/BF00337288](https://doi.org/10.1007/BF00337288) .

[kohonen on CRAN](https://CRAN.R-project.org/package=kohonen).

## Examples

``` r
exp <- glyexp::real_experiment2
res <- gly_som(exp, grid = c(2, 2), seed = 42, add_info = FALSE)
head(get_tidy_result(res))
#> # A tibble: 6 × 7
#>   .variable                unit     x     y distance n_finite exclusion         
#>   <chr>                   <int> <int> <int>    <dbl>    <int> <chr>             
#> 1 Hex(3)HexNAc(3)            NA    NA    NA    NA         138 incomplete_profile
#> 2 Hex(3)HexNAc(7)            NA    NA    NA    NA          94 incomplete_profile
#> 3 Hex(5)HexNAc(2)             2     2     1     2.46      144 NA                
#> 4 Hex(6)HexNAc(4)NeuAc(2)    NA    NA    NA    NA         136 incomplete_profile
#> 5 Hex(4)HexNAc(3)            NA    NA    NA    NA         143 incomplete_profile
#> 6 Hex(5)HexNAc(4)dHex(2)     NA    NA    NA    NA          64 incomplete_profile
res$diagnostics$units
#> # A tibble: 4 × 4
#>    unit     x     y n_profiles
#>   <int> <int> <int>      <int>
#> 1     1     1     1          2
#> 2     2     2     1          7
#> 3     3     1     2         11
#> 4     4     2     2          8
```
