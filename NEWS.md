# glystats (development version)

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
