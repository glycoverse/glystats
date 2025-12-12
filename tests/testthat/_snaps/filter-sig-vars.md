# filter_sig_vars fails for ANOVA with fc_cutoff on main test

    Code
      filter_sig_vars(exp, result, fc_cutoff = 2)
    Condition
      Error in `.check_filter_sig_vars_args_anova()`:
      ! `fc_cutoff` can't be used with `gly_anova()` or `gly_kruskal()` when `on` is "main_test".
      i Please set `fc_cutoff` to NULL or `on` to "post_hoc_test".

# filter_sig_vars fails for ANOVA with comparison on main test

    Code
      filter_sig_vars(exp, result, comparison = "C_vs_H")
    Condition
      Error in `.check_filter_sig_vars_args_anova()`:
      ! `comparison` can't be used with `gly_anova()` or `gly_kruskal()` when `on` is "main_test".
      i Please set `comparison` to NULL or `on` to "post_hoc_test".

# filter_sig_vars fails for ANOVA with both p_val_cutoff and p_adj_cutoff

    Code
      filter_sig_vars(exp, result, p_val_cutoff = 0.05, p_adj_cutoff = 0.05)
    Condition
      Error in `.check_filter_sig_vars_args()`:
      ! Only one of `p_adj_cutoff` or `p_val_cutoff` can be provided.
      x Both are provided.

# filter_sig_vars fails for ANOVA with neither p_val_cutoff, p_adj_cutoff

    Code
      filter_sig_vars(exp, result, p_val_cutoff = NULL, p_adj_cutoff = NULL)
    Condition
      Error in `.check_filter_sig_vars_args()`:
      ! At least one of `p_adj_cutoff`, `p_val_cutoff`, or `fc_cutoff` must be provided.
      x All are NULL.

# filter_sig_vars fails for ANOVA on post-hoc results without comparison

    Code
      filter_sig_vars(exp, result, on = "post_hoc_test")
    Condition
      Error in `.check_filter_sig_vars_args_anova()`:
      ! `comparison` must be provided when `on` is "post_hoc_test".
      i Please set `comparison` to a string with the format of "group1_vs_group2".

# filter_sig_vars fails for ANOVA with invalid comparison format

    Code
      filter_sig_vars(exp, result, on = "post_hoc_test", comparison = "H_vs_C_vs_M")
    Condition
      Error in `.check_filter_sig_vars_args()`:
      ! `comparison` must be in the format of "group1_vs_group2".
      x Invalid format: "H_vs_C_vs_M".

# filter_sig_vars fails for ANOVA with non-existent comparison

    Code
      filter_sig_vars(exp, result, on = "post_hoc_test", comparison = "H_vs_D")
    Condition
      Error in `.check_filter_sig_vars_args_anova()`:
      ! Can't find comparison: "H_vs_D".
      i Available comparisons: "H_vs_M", "H_vs_Y", "H_vs_C", "M_vs_Y", "M_vs_C", and "Y_vs_C".

# filter_sig_vars fails for ANOVA with comparison of wrong order

    Code
      filter_sig_vars(exp, result, on = "post_hoc_test", comparison = "C_vs_H")
    Condition
      Error in `.check_filter_sig_vars_args_anova()`:
      ! Can't find comparison: "C_vs_H".
      i Available comparisons: "H_vs_M", "H_vs_Y", "H_vs_C", "M_vs_Y", "M_vs_C", and "Y_vs_C".
      i Did you mean "H_vs_C"?

