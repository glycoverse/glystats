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
