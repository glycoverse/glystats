# Weighted Gene Co-expression Network Analysis (WGCNA)

Perform WGCNA analysis to identify co-expression modules and their
eigenvalues. The function uses WGCNA package to construct weighted gene
co-expression networks, detect modules, and calculate module membership
and eigenvalues.

## Usage

``` r
gly_wgcna(
  exp,
  powers = c(1:10, seq(12, 20, by = 2)),
  network_type = "unsigned",
  tom_type = "unsigned",
  min_module_size = 30,
  deep_split = 2,
  merge_cut_height = 0.25,
  add_info = TRUE,
  ...
)

gly_wgcna_(
  expr_mat,
  powers = c(1:10, seq(12, 20, by = 2)),
  network_type = "unsigned",
  tom_type = "unsigned",
  min_module_size = 30,
  deep_split = 2,
  merge_cut_height = 0.25,
  ...
)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- powers:

  A numeric vector of soft thresholding powers to test. Default is
  `c(1:10, seq(12, 20, by = 2))`.

- network_type:

  Character string specifying network type. Allowed values are
  "unsigned", "signed", "signed hybrid". Default is "unsigned". Passed
  to the ``` NetworkType`` argument of  ```WGCNA::blockwiseModules()\`.

- tom_type:

  Character string specifying topological overlap matrix type. Allowed
  values are "none", "unsigned", "signed". Default is "unsigned". Passed
  to the ``` TOMType`` argument of  ```WGCNA::blockwiseModules()\`.

- min_module_size:

  Minimum module size for module detection. Default is 30. Passed to the
  ``` minModuleSize`` argument of  ```WGCNA::blockwiseModules()\`.

- deep_split:

  Integer value between 0 and 4. Provides a simplified control over how
  sensitive module detection should be to module splitting. Default
  is 2. Passed to the
  ``` deepSplit`` argument of  ```WGCNA::blockwiseModules()\`.

- merge_cut_height:

  Dendrogram cut height for merging of modules. Default is 0.25. Passed
  to the
  ``` mergeCutHeight`` argument of  ```WGCNA::blockwiseModules()\`.

- add_info:

  A logical value indicating whether to add variable information to the
  modules tibble. Only applicable to `gly_wgcna()`.

- ...:

  Additional arguments passed to
  [`WGCNA::blockwiseModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseModules.html).

- expr_mat:

  A numeric matrix with variables as rows and samples as columns.

## Value

A list with two elements:

- `tidy_result`: A list containing two tibbles:

  - `modules`: Module assignments and membership values containing the
    following columns:

    - `variable`: Variable name

    - `module`: Module assignment (color name)

    - `membership`: Module membership value (correlation with module
      eigengene)

  - `eigenvalues`: Module eigenvalues containing the following columns:

    - `module`: Module name (color name)

    - `sample`: Sample name

    - `eigenvalue`: Module eigenvalue (first principal component of
      module expression)

- `raw_result`: The raw WGCNA blockwiseModules object The list has
  classes `glystats_wgcna_res` and `glystats_res`.

## Details

The function performs log2 transformation on the expression data
(log2(x + 1)) before WGCNA analysis.

`gly_wgcna()` is the top-level API that works with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and supports the `add_info` parameter for joining experiment
metadata.

`gly_wgcna_()` is the underlying API that works with matrices directly,
providing more flexibility for users who don't use the glyexp package.

**Analysis Steps:**

1.  Soft threshold selection using
    [`WGCNA::pickSoftThreshold()`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html)

2.  Network construction and module detection using
    [`WGCNA::blockwiseModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseModules.html)

3.  Module membership calculation based on correlation with module
    eigengenes

4.  Results formatting into standardized tibbles

## Required packages

This function requires the following packages to be installed:

- WGCNA\` for weighted gene co-expression network analysis

## See also

[`WGCNA::pickSoftThreshold()`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html),
[`WGCNA::blockwiseModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseModules.html)
