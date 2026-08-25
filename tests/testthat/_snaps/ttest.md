# paired analyses reject ambiguous subject layouts

    Code
      .match_paired_samples(expression, groups, subjects = c("S1", "S1", "S2", "S1",
        "S2", "S3"), method = "t-test")
    Condition
      Error in `.match_paired_samples()`:
      ! Each subject must have at most one sample in each group for paired t-test.

---

    Code
      .match_paired_samples(expression, groups, subjects = c("S1", "S2", "S3", "S3",
        "S4", "S5"), method = "t-test")
    Condition
      Error in `.match_paired_samples()`:
      ! At least 2 subjects must have samples in every group for paired t-test.
