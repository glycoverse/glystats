# glystats (development version)

* `gly_anova()`, `gly_kruskal()`, `gly_ttest()`, and `gly_wilcox()` gain a `subject_col` argument for subject-matched paired analyses, including repeated-measures ANOVA, Friedman tests, paired post-hoc comparisons, and paired effect sizes.
* `gly_set_test()` now constructs and tests correlated variable sets in one call, retains isolated variables as singleton sets, and tests structurally rank-deficient sets in their nonredundant subspace when the observed contrast is estimable and classical sample-size requirements are met, without distinguishing identical profiles; custom sets remain available through `sets`. (#20)
* Analysis functions now accept experiments containing all-`NA` variables without errors, returning `NA` feature-level statistics, loadings, or clusters and `NULL` raw models when no variables can be fitted. (#17)
* `gly_linear_model()` now fits formula-based, limma-moderated models with interactions, adjustment variables, and named custom contrasts. (#16)

# glystats 0.11.2

* Documentation and vignettes now recommend `GlycomicSE` and `GlycoproteomicSE` containers with `SummarizedExperiment` accessors for Stage II of glycoverse/glyexp#15. (#15)

# glystats 0.11.1

* Statistical analysis functions now recognize legacy `glyexp::experiment()` and current `SummarizedExperiment` containers consistently across supported `glyexp` versions.

# glystats 0.11.0

## Breaking changes

* All `gly_*_()` matrix and vector interfaces are removed. Use the corresponding `gly_*()` function with a `glyexp::experiment()` or `SummarizedExperiment` object instead. (#12)
* `gly_enrich_go()`, `gly_enrich_kegg()`, `gly_enrich_reactome()`, `gly_enrich_wikipathways()`, `gly_enrich_do()`, `gly_enrich_ncg()`, and `gly_enrich_dgn()` are removed. These functions were deprecated in 0.10.0. Use the corresponding functions in the `glyfun` package instead. (#13)

## New features

* All `gly_*()` analysis functions and `filter_sig_vars()` now accept `SummarizedExperiment` objects in addition to `glyexp::experiment()` objects. (#14)

# glystats 0.10.1

## Minor improvements and bug fixes

* Use a smaller log2 pseudo-count across DEA and clustering functions to reduce bias for low-intensity values (#10).

# glystats 0.10.0

## Breaking changes

* `gly_wgcna()`, `gly_wgcna_()`, `gly_consensus_clustering()`, and `gly_consensus_clustering_()` are removed. These functions were deprecated in 0.7.0. Use the `WGCNA` and `ConsensusClusterPlus` packages directly.

## Deprecations

* All enrichment analysis functions are now deprecated and will be removed in a future version: `gly_enrich_go()`, `gly_enrich_kegg()`, `gly_enrich_reactome()`, `gly_enrich_wikipathways()`, `gly_enrich_do()`, `gly_enrich_ncg()`, and `gly_enrich_dgn()`. Use the corresponding functions in the `glyfun` package instead.

# glystats 0.9.0

## Breaking changes

* `get_tidy_result()` is now type-stable: it always returns a tibble (or errors), instead of sometimes returning a list when `which` is omitted. When `tidy_result` is a list, the `which` argument is now required (#9).

## New features

* Add effect size calculations to `gly_ttest()`, `gly_wilcox()`, `gly_anova()`, and `gly_kruskal()`. Effect sizes (Cohen's d, rank-biserial correlation, eta-squared, and epsilon-squared) are now included in the tidy results (#8).

## Minor improvements and bug fixes

* Fix a bug that caused statistics to have incorrect signs for `gly_ttest()` and `gly_wilcox()`.
* Fix a bug that caused the `estimate` sign to be inconsistent in `gly_kruskal()` results.

# glystats 0.8.0

## New features

* Add `meta_data` to all glystats results. This allows storing additional information with the result object that can be accessed later (#7).

## Minor improvements and bug fixes

* Remove duplicate reference group message in `gly_ttest()` output.
* Update vignette function list.

# glystats 0.7.0

## Breaking changes

* Refactored enrichment API (`gly_enrich_go()`, `gly_enrich_kegg()`, `gly_enrich_reactome()`). These functions now use explicit parameters instead of `...`. Users must update their code to use new parameters (e.g., `p_cutoff` instead of `pvalueCutoff`).
* `gly_wgcna()`, `gly_wgcna_()`, `gly_consensus_clustering()`, and `gly_consensus_clustering_()` are deprecated and will be removed in a future version. These functions are too interactive and complex for a pipeline-friendly package.

## New features

* Add `gly_enrich_wikipathways()` and `gly_enrich_wikipathways_()` for WikiPathways enrichment analysis.
* Add `gly_enrich_do()` and `gly_enrich_do_()` for Disease Ontology enrichment analysis.
* Add `gly_enrich_ncg()` and `gly_enrich_ncg_()` for Network of Cancer Genes enrichment analysis.
* Add `gly_enrich_dgn()` and `gly_enrich_dgn_()` for DisGeNET enrichment analysis.

# glystats 0.6.5

## Minor improvements and bug fixes

* Update dependency strategy to use the r-universe repo.

# glystats 0.6.4

## Minor improvements and bug fixes

* Fix the bug that `gly_hclust()` leaves a plotting device unclosed.

# glystats 0.6.3

## Minor improvements and bug fixes

* Fix the bug that `gly_limma()` fails for site-specific motif experiments.

# glystats 0.6.2

## Breaking changes

* The `subject_cols` argument in `gly_limma()` has been renamed into `subject_col`.

# glystats 0.6.1

## Minor improvements and bug fixes

* Fix the bug that `gly_limma()` fails to work with experiments with all-NA variables.

# glystats 0.6.0

## New features

* Add `gly_ancova()` for Analysis of Covariance (ANCOVA).
* Add covariate support to `gly_limma()` and `gly_limma_()` through the `covariate_cols` and the `covariates` arguments.

## Minor improvements and bug fixes

* Use the "limma-trend" method in `gly_limma()` by default. This is more appropriate for glycomics and glycoproteomics data.
* Add CI to AUC values in the result of `gly_roc()`.
* Update documentations about arguments.

# glystats 0.5.8

## Breaking changes

* `filter_sig_vars()` now has a default `fc_cutoff` of `2` for glycoproteomics data, and `NULL` for other types of data.

# glystats 0.5.7

## Minor improvements and bug fixes

* Fix the bug that `gly_pca()` fails to work with experiments with constant variables. Now the constant variables are removed before PCA, with a warning issued.

# glystats 0.5.6

## Minor improvements and bug fixes

* Fix the bug that fold change estimation in `gly_limma()` is not correct for glycomics data.

# glystats 0.5.5

## Minor improvements and bug fixes

* The `universe` argument of `gly_enrich_go()`, `gly_enrich_kegg()`, and `gly_enrich_reactome()` accepts a `glyexp::experiment()` object now.
* Optimize message output of `gly_limma()`.
* The columns from sample information or variable information are joined right after the `sample` or `variable` column in the results now, not appending to the end as previously.
* Update documentation of `filter_sig_vars()` to inform about using the same experiment used for DEA.

# glystats 0.5.4

## Minor improvements and bug fixes

* `filter_sig_vars()` now doesn't require `comparison` to be provided for `gly_limma()` results. If `comparison` is not provided, a variable will be kept if it is significant in any comparison.

# glystats 0.5.3

## Minor improvements and fixes

* glystats now depends on the CRAN version of glyrepr.

# glystats 0.5.2

## Minor improvements and bug fixes

* Fix the bug in many functions that some parameters cannot be passed by `...` to the underlying statistical functions.

# glystats 0.5.1

## Minor improvements and bug fixes

* Fix bugs introduced by the breaking changes in `glyexp` 0.10.0.

# glystats 0.5.0

This is a big update! We make some breaking changes to the API, mainly DEA functions.
Now you can ensure these properties for all DEA results:

* The levels in the front are the reference group. For example, for groups with levels c("A", "B"), "A" is the reference group. For groups with levels c("A", "B", "C"), "A" is the reference group in A-B and A-C comparisons, and "B" is the reference group in B-C comparisons.
* For the concept of "comparison", either in result columns or as arguments, it is expected to be in the format of "ref_vs_test".

## Breaking changes

* The `post_hoc` column in `tidy_result$main_test` of results from `gly_anova()` and `gly_kruskal()` is now in the format of "ref_vs_test" instead of "ref-test".
* The `group1` and `group2` columns in `tidy_result$post_hoc_test` of results from `gly_anova()` and `gly_kruskal()` are renamed into `ref_group` and `test_group` for more clarity.
* In the `tidy_result` of `gly_limma()`, the `comparison` column is replaced by the `ref_group` and `test_group` columns.
* The results of all functions now have a `p_val` column for raw p-values, a `p_adj` column for adjusted p-values, and a `log2fc` column for log2 fold change, if applicable. This consistence in column naming reduces the cognitive load.

## New features

* Add `get_tidy_result()` and `get_raw_result()` to get the tidy result tibble and the raw result list from a glystats result object. These functions are useful to be used in pipes.
* Add `filter_sig_vars()` to filter the experiment using the results from glystats DEA functions to keep only significant variables.
* `gly_fold_change()` now supports multiple groups.
* The post-hoc results from `gly_anova()` and `gly_kruskal()` now have a `log2fc` column for log2 fold change.

## Minor improvements and bug fixes

* Update group information message in many functions. Instead of "Group 1" and "Group 2", now we use "Ref Group" and "Test Group" for more clarity.
* There are `ref_group` and `test_group` columns in the `tidy_result` tibble of `gly_limma()` even for 2 groups now.
* Fix the bug that `gly_kruskal()` failed to work with experiments with 2 groups.
* Fix the bug that the direction of fold change is not consistent when with 2 groups and with multiple groups.
* Fix the bug that `gly_anova()` and `gly_kruskal()` have NAs in the `post_hoc_test` tibble.
* Add `glyrepr` to dependencies to fix the result printing bug.

# glystats 0.4.2

## Minor improvements and bug fixes

* Update dependencies to depend on release versions of glycoverse packages.

# glystats 0.4.1

## Minor improvements and bug fixes

* Fix the `could not find function "%||%"` bug in `gly_anova()`.
* `gly_anova()`, `gly_kruskal()`, `gly_ttest()`, `gly_wilcox()`, `gly_cox()`, `gly_roc()`, and their `gly_xxx_()` counterparts are now more robust. Previously, if the model failed to be fitted for any variable, the functions stopped with an error. Now, they only issue a warning and continue, assigning NAs in the results for those failed variables.

# glystats 0.4.0

## Breaking changes

* `gly_plsda()` now uses the `ropls` package as its backend instead of `mixOmics`. This changes the class of the object returned in the `raw_result` list element, which may affect downstream code that uses the raw result directly.

## New features

* Added permutation testing to `gly_oplsda()` via the new `perm_test` parameter to assess model significance.
* The result from `gly_oplsda()` now includes a `pcorr` column containing p-values for the correlation coefficients of features with the predictive component.

## Minor improvements and bug fixes

* Removed the check for n/p in `gly_oplsda()` and `gly_plsda()`, increasing flexibility for smaller datasets.
* Improved the error message when a specified grouping column is not found in the sample information.
* Fixed a bug in `gly_limma()` that could cause duplicated messages to be printed.
* Fixed a bug in `gly_oplsda()` where the column for the first orthogonal component score (`o1`) could be missing from the result.
* Suppressed messages from dependency packages and extra blank lines when running `gly_enrich_go()` for a cleaner console output.
* Removed redundant `glystats` S3 classes from the tibbles in the `tidy_result`, simplifying the output object structure.

# glystats 0.3.0

## Breaking changes

This version introduces new API for all functions.
Briefly, the `return_raw` parameter is removed,
and all functions now return a list with two elements:
`tidy_result` and `raw_result`.
The concrete types of `tidy_result` and `raw_result` depend.
`tidy_result` can be a tibble, or a list of tibbles.
`raw_result` can be a single object returned by the underlying statistical functions,
or a list of such objects.
This update makes `glystats` easier to use.
And more importantly, it allows the `glyvis` package to access the raw results directly.

## Minor improvements and bug fixes

* Update the documentation of all functions to include the detailed column descriptions
  for the tibbles in `tidy_result`.
* All `gly_xxx_()` functions now accept a character vector as the `groups` parameter.
* Fix an issue that `gly_consensus_clustering()` sends plots to the plot panel when `output_file` is NULL.
  This is an inconsistent behavior compared to other functions, and it has been fixed.
* Update the documentation of `gly_consensus_clustering()` to emphasize the importance of `output_file`.
* Add an introduction vignette.

# glystats 0.2.3

## Minor improvements and bug fixes

* Make the parameters of `gly_umap()` and `gly_umap_()` consistent.
* Fix a bug that some functions returned values with duplicated S3 classes.
* Fix inconsistent behaviours between `gly_tsne()` and `gly_tsne_()`.

# glystats 0.2.2

## Bug fixes

* Fix a bug in test.

# glystats 0.2.1

## Minor improvements

* Change S3 class of the results of `gly_anova()`, `gly_kruskal()`, `gly_limma()`, `gly_ttest()`, and `gly_wilcox()`.

# glystats 0.2.0

## Breaking changes

* `gly_anova()` and `gly_kruskal()` now return lists of two tibbles, one for main test and the other for post-hoc test.
* Remove the `method`, `dist_method` parameters from `gly_hclust()`.
* Remove the `max_iter`, `theta` parameters from `gly_tsne()`.
* Remove the `n_epochs`, `learning_rate`, `metric` parameters from `gly_umap()`.

## New features

* Add a new set of API functions: all `gly_xxx()` functions now have a lower API counterpart `gly_xxx_()`
  that works with matrices directly, providing more flexibility for users who don't use the glyexp package.
* Add `gly_limma()` and `gly_limma_()` to perform differential analysis using the limma package.
* Add `gly_kmeans()` and `gly_kmeans_()` to perform K-means clustering.
* Add `gly_wgcna()` and `gly_wgcna_()` to perform WGCNA analysis.
* Add `gly_consensus_clustering()` and `gly_consensus_clustering_()` to perform consensus clustering.
* Add `gly_cox()` and `gly_cox_()` to fit Cox proportional hazards model for survival analysis.
* Add `ref_group` parameter to `gly_ttext()` and `gly_wilcox()`.

## Minor improvements and bug fixes

* Fix the bugs in many functions that parameter validations were not performed.
* Reduce size of the test data.
* Remove dependency on the `parameters` package.
