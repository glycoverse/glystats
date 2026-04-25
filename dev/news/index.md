# Changelog

## glystats (development version)

## glystats 0.10.0

### Breaking changes

- `gly_wgcna()`, `gly_wgcna_()`, `gly_consensus_clustering()`, and
  `gly_consensus_clustering_()` are removed. These functions were
  deprecated in 0.7.0. Use the `WGCNA` and `ConsensusClusterPlus`
  packages directly.

### Deprecations

- All enrichment analysis functions are now deprecated and will be
  removed in a future version:
  [`gly_enrich_go()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_go.md),
  [`gly_enrich_kegg()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_kegg.md),
  [`gly_enrich_reactome()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_reactome.md),
  [`gly_enrich_wikipathways()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_wikipathways.md),
  [`gly_enrich_do()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_do.md),
  [`gly_enrich_ncg()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_ncg.md),
  and
  [`gly_enrich_dgn()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_dgn.md).
  Use the corresponding functions in the `glyfun` package instead.

## glystats 0.9.0

### Breaking changes

- [`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
  is now type-stable: it always returns a tibble (or errors), instead of
  sometimes returning a list when `which` is omitted. When `tidy_result`
  is a list, the `which` argument is now required (#9).

### New features

- Add effect size calculations to
  [`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md),
  [`gly_wilcox()`](https://glycoverse.github.io/glystats/dev/reference/gly_wilcox.md),
  [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md),
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md).
  Effect sizes (Cohen’s d, rank-biserial correlation, eta-squared, and
  epsilon-squared) are now included in the tidy results (#8).

### Minor improvements and bug fixes

- Fix a bug that caused statistics to have incorrect signs for
  [`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md)
  and
  [`gly_wilcox()`](https://glycoverse.github.io/glystats/dev/reference/gly_wilcox.md).
- Fix a bug that caused the `estimate` sign to be inconsistent in
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md)
  results.

## glystats 0.8.0

### New features

- Add `meta_data` to all glystats results. This allows storing
  additional information with the result object that can be accessed
  later (#7).

### Minor improvements and bug fixes

- Remove duplicate reference group message in
  [`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md)
  output.
- Update vignette function list.

## glystats 0.7.0

### Breaking changes

- Refactored enrichment API
  ([`gly_enrich_go()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_go.md),
  [`gly_enrich_kegg()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_kegg.md),
  [`gly_enrich_reactome()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_reactome.md)).
  These functions now use explicit parameters instead of `...`. Users
  must update their code to use new parameters (e.g., `p_cutoff` instead
  of `pvalueCutoff`).
- `gly_wgcna()`, `gly_wgcna_()`, `gly_consensus_clustering()`, and
  `gly_consensus_clustering_()` are deprecated and will be removed in a
  future version. These functions are too interactive and complex for a
  pipeline-friendly package.

### New features

- Add
  [`gly_enrich_wikipathways()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_wikipathways.md)
  and
  [`gly_enrich_wikipathways_()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_wikipathways.md)
  for WikiPathways enrichment analysis.
- Add
  [`gly_enrich_do()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_do.md)
  and
  [`gly_enrich_do_()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_do.md)
  for Disease Ontology enrichment analysis.
- Add
  [`gly_enrich_ncg()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_ncg.md)
  and
  [`gly_enrich_ncg_()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_ncg.md)
  for Network of Cancer Genes enrichment analysis.
- Add
  [`gly_enrich_dgn()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_dgn.md)
  and
  [`gly_enrich_dgn_()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_dgn.md)
  for DisGeNET enrichment analysis.

## glystats 0.6.5

### Minor improvements and bug fixes

- Update dependency strategy to use the r-universe repo.

## glystats 0.6.4

### Minor improvements and bug fixes

- Fix the bug that
  [`gly_hclust()`](https://glycoverse.github.io/glystats/dev/reference/gly_hclust.md)
  leaves a plotting device unclosed.

## glystats 0.6.3

### Minor improvements and bug fixes

- Fix the bug that
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  fails for site-specific motif experiments.

## glystats 0.6.2

### Breaking changes

- The `subject_cols` argument in
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  has been renamed into `subject_col`.

## glystats 0.6.1

### Minor improvements and bug fixes

- Fix the bug that
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  fails to work with experiments with all-NA variables.

## glystats 0.6.0

### New features

- Add
  [`gly_ancova()`](https://glycoverse.github.io/glystats/dev/reference/gly_ancova.md)
  for Analysis of Covariance (ANCOVA).
- Add covariate support to
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  and
  [`gly_limma_()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  through the `covariate_cols` and the `covariates` arguments.

### Minor improvements and bug fixes

- Use the “limma-trend” method in
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  by default. This is more appropriate for glycomics and glycoproteomics
  data.
- Add CI to AUC values in the result of
  [`gly_roc()`](https://glycoverse.github.io/glystats/dev/reference/gly_roc.md).
- Update documentations about arguments.

## glystats 0.5.8

### Breaking changes

- [`filter_sig_vars()`](https://glycoverse.github.io/glystats/dev/reference/filter_sig_vars.md)
  now has a default `fc_cutoff` of `2` for glycoproteomics data, and
  `NULL` for other types of data.

## glystats 0.5.7

### Minor improvements and bug fixes

- Fix the bug that
  [`gly_pca()`](https://glycoverse.github.io/glystats/dev/reference/gly_pca.md)
  fails to work with experiments with constant variables. Now the
  constant variables are removed before PCA, with a warning issued.

## glystats 0.5.6

### Minor improvements and bug fixes

- Fix the bug that fold change estimation in
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  is not correct for glycomics data.

## glystats 0.5.5

### Minor improvements and bug fixes

- The `universe` argument of
  [`gly_enrich_go()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_go.md),
  [`gly_enrich_kegg()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_kegg.md),
  and
  [`gly_enrich_reactome()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_reactome.md)
  accepts a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object now.
- Optimize message output of
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md).
- The columns from sample information or variable information are joined
  right after the `sample` or `variable` column in the results now, not
  appending to the end as previously.
- Update documentation of
  [`filter_sig_vars()`](https://glycoverse.github.io/glystats/dev/reference/filter_sig_vars.md)
  to inform about using the same experiment used for DEA.

## glystats 0.5.4

### Minor improvements and bug fixes

- [`filter_sig_vars()`](https://glycoverse.github.io/glystats/dev/reference/filter_sig_vars.md)
  now doesn’t require `comparison` to be provided for
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  results. If `comparison` is not provided, a variable will be kept if
  it is significant in any comparison.

## glystats 0.5.3

### Minor improvements and fixes

- glystats now depends on the CRAN version of glyrepr.

## glystats 0.5.2

### Minor improvements and bug fixes

- Fix the bug in many functions that some parameters cannot be passed by
  `...` to the underlying statistical functions.

## glystats 0.5.1

### Minor improvements and bug fixes

- Fix bugs introduced by the breaking changes in `glyexp` 0.10.0.

## glystats 0.5.0

This is a big update! We make some breaking changes to the API, mainly
DEA functions. Now you can ensure these properties for all DEA results:

- The levels in the front are the reference group. For example, for
  groups with levels c(“A”, “B”), “A” is the reference group. For groups
  with levels c(“A”, “B”, “C”), “A” is the reference group in A-B and
  A-C comparisons, and “B” is the reference group in B-C comparisons.
- For the concept of “comparison”, either in result columns or as
  arguments, it is expected to be in the format of “ref_vs_test”.

### Breaking changes

- The `post_hoc` column in `tidy_result$main_test` of results from
  [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md)
  is now in the format of “ref_vs_test” instead of “ref-test”.
- The `group1` and `group2` columns in `tidy_result$post_hoc_test` of
  results from
  [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md)
  are renamed into `ref_group` and `test_group` for more clarity.
- In the `tidy_result` of
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md),
  the `comparison` column is replaced by the `ref_group` and
  `test_group` columns.
- The results of all functions now have a `p_val` column for raw
  p-values, a `p_adj` column for adjusted p-values, and a `log2fc`
  column for log2 fold change, if applicable. This consistence in column
  naming reduces the cognitive load.

### New features

- Add
  [`get_tidy_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
  and
  [`get_raw_result()`](https://glycoverse.github.io/glystats/dev/reference/get_tidy_result.md)
  to get the tidy result tibble and the raw result list from a glystats
  result object. These functions are useful to be used in pipes.
- Add
  [`filter_sig_vars()`](https://glycoverse.github.io/glystats/dev/reference/filter_sig_vars.md)
  to filter the experiment using the results from glystats DEA functions
  to keep only significant variables.
- [`gly_fold_change()`](https://glycoverse.github.io/glystats/dev/reference/gly_fold_change.md)
  now supports multiple groups.
- The post-hoc results from
  [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md)
  now have a `log2fc` column for log2 fold change.

### Minor improvements and bug fixes

- Update group information message in many functions. Instead of “Group
  1” and “Group 2”, now we use “Ref Group” and “Test Group” for more
  clarity.
- There are `ref_group` and `test_group` columns in the `tidy_result`
  tibble of
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  even for 2 groups now.
- Fix the bug that
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md)
  failed to work with experiments with 2 groups.
- Fix the bug that the direction of fold change is not consistent when
  with 2 groups and with multiple groups.
- Fix the bug that
  [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md)
  have NAs in the `post_hoc_test` tibble.
- Add `glyrepr` to dependencies to fix the result printing bug.

## glystats 0.4.2

### Minor improvements and bug fixes

- Update dependencies to depend on release versions of glycoverse
  packages.

## glystats 0.4.1

### Minor improvements and bug fixes

- Fix the `could not find function "%||%"` bug in
  [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md).
- [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md),
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md),
  [`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md),
  [`gly_wilcox()`](https://glycoverse.github.io/glystats/dev/reference/gly_wilcox.md),
  [`gly_cox()`](https://glycoverse.github.io/glystats/dev/reference/gly_cox.md),
  [`gly_roc()`](https://glycoverse.github.io/glystats/dev/reference/gly_roc.md),
  and their `gly_xxx_()` counterparts are now more robust. Previously,
  if the model failed to be fitted for any variable, the functions
  stopped with an error. Now, they only issue a warning and continue,
  assigning NAs in the results for those failed variables.

## glystats 0.4.0

### Breaking changes

- [`gly_plsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_plsda.md)
  now uses the `ropls` package as its backend instead of `mixOmics`.
  This changes the class of the object returned in the `raw_result` list
  element, which may affect downstream code that uses the raw result
  directly.

### New features

- Added permutation testing to
  [`gly_oplsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_oplsda.md)
  via the new `perm_test` parameter to assess model significance.
- The result from
  [`gly_oplsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_oplsda.md)
  now includes a `pcorr` column containing p-values for the correlation
  coefficients of features with the predictive component.

### Minor improvements and bug fixes

- Removed the check for n/p in
  [`gly_oplsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_oplsda.md)
  and
  [`gly_plsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_plsda.md),
  increasing flexibility for smaller datasets.
- Improved the error message when a specified grouping column is not
  found in the sample information.
- Fixed a bug in
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  that could cause duplicated messages to be printed.
- Fixed a bug in
  [`gly_oplsda()`](https://glycoverse.github.io/glystats/dev/reference/gly_oplsda.md)
  where the column for the first orthogonal component score (`o1`) could
  be missing from the result.
- Suppressed messages from dependency packages and extra blank lines
  when running
  [`gly_enrich_go()`](https://glycoverse.github.io/glystats/dev/reference/gly_enrich_go.md)
  for a cleaner console output.
- Removed redundant `glystats` S3 classes from the tibbles in the
  `tidy_result`, simplifying the output object structure.

## glystats 0.3.0

### Breaking changes

This version introduces new API for all functions. Briefly, the
`return_raw` parameter is removed, and all functions now return a list
with two elements: `tidy_result` and `raw_result`. The concrete types of
`tidy_result` and `raw_result` depend. `tidy_result` can be a tibble, or
a list of tibbles. `raw_result` can be a single object returned by the
underlying statistical functions, or a list of such objects. This update
makes `glystats` easier to use. And more importantly, it allows the
`glyvis` package to access the raw results directly.

### Minor improvements and bug fixes

- Update the documentation of all functions to include the detailed
  column descriptions for the tibbles in `tidy_result`.
- All `gly_xxx_()` functions now accept a character vector as the
  `groups` parameter.
- Fix an issue that `gly_consensus_clustering()` sends plots to the plot
  panel when `output_file` is NULL. This is an inconsistent behavior
  compared to other functions, and it has been fixed.
- Update the documentation of `gly_consensus_clustering()` to emphasize
  the importance of `output_file`.
- Add an introduction vignette.

## glystats 0.2.3

### Minor improvements and bug fixes

- Make the parameters of
  [`gly_umap()`](https://glycoverse.github.io/glystats/dev/reference/gly_umap.md)
  and
  [`gly_umap_()`](https://glycoverse.github.io/glystats/dev/reference/gly_umap.md)
  consistent.
- Fix a bug that some functions returned values with duplicated S3
  classes.
- Fix inconsistent behaviours between
  [`gly_tsne()`](https://glycoverse.github.io/glystats/dev/reference/gly_tsne.md)
  and
  [`gly_tsne_()`](https://glycoverse.github.io/glystats/dev/reference/gly_tsne.md).

## glystats 0.2.2

### Bug fixes

- Fix a bug in test.

## glystats 0.2.1

### Minor improvements

- Change S3 class of the results of
  [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md),
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md),
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md),
  [`gly_ttest()`](https://glycoverse.github.io/glystats/dev/reference/gly_ttest.md),
  and
  [`gly_wilcox()`](https://glycoverse.github.io/glystats/dev/reference/gly_wilcox.md).

## glystats 0.2.0

### Breaking changes

- [`gly_anova()`](https://glycoverse.github.io/glystats/dev/reference/gly_anova.md)
  and
  [`gly_kruskal()`](https://glycoverse.github.io/glystats/dev/reference/gly_kruskal.md)
  now return lists of two tibbles, one for main test and the other for
  post-hoc test.
- Remove the `method`, `dist_method` parameters from
  [`gly_hclust()`](https://glycoverse.github.io/glystats/dev/reference/gly_hclust.md).
- Remove the `max_iter`, `theta` parameters from
  [`gly_tsne()`](https://glycoverse.github.io/glystats/dev/reference/gly_tsne.md).
- Remove the `n_epochs`, `learning_rate`, `metric` parameters from
  [`gly_umap()`](https://glycoverse.github.io/glystats/dev/reference/gly_umap.md).

### New features

- Add a new set of API functions: all `gly_xxx()` functions now have a
  lower API counterpart `gly_xxx_()` that works with matrices directly,
  providing more flexibility for users who don’t use the glyexp package.
- Add
  [`gly_limma()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  and
  [`gly_limma_()`](https://glycoverse.github.io/glystats/dev/reference/gly_limma.md)
  to perform differential analysis using the limma package.
- Add
  [`gly_kmeans()`](https://glycoverse.github.io/glystats/dev/reference/gly_kmeans.md)
  and
  [`gly_kmeans_()`](https://glycoverse.github.io/glystats/dev/reference/gly_kmeans.md)
  to perform K-means clustering.
- Add `gly_wgcna()` and `gly_wgcna_()` to perform WGCNA analysis.
- Add `gly_consensus_clustering()` and `gly_consensus_clustering_()` to
  perform consensus clustering.
- Add
  [`gly_cox()`](https://glycoverse.github.io/glystats/dev/reference/gly_cox.md)
  and
  [`gly_cox_()`](https://glycoverse.github.io/glystats/dev/reference/gly_cox.md)
  to fit Cox proportional hazards model for survival analysis.
- Add `ref_group` parameter to `gly_ttext()` and
  [`gly_wilcox()`](https://glycoverse.github.io/glystats/dev/reference/gly_wilcox.md).

### Minor improvements and bug fixes

- Fix the bugs in many functions that parameter validations were not
  performed.
- Reduce size of the test data.
- Remove dependency on the `parameters` package.
