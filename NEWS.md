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