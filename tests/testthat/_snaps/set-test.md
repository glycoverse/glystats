# set construction and testing validate metadata contracts

    Code
      gly_set_test(exp, within = "unknown")
    Condition
      Error in `.set_strata()`:
      ! Columns in `within` were not found in variable information.
      x Missing columns: unknown.
      i Available columns: variable and site.

---

    Code
      gly_set_test(exp, list(signal = c("A", "unknown")))
    Condition
      Error in `.normalize_set_definitions()`:
      ! Some variables in `sets` were not found in `exp`.
      x Missing variables: "unknown".

---

    Code
      gly_set_test(paired, list(signal = c("A", "B")), subject_col = "subject")
    Message
      i Ref Group: "control"
      i Test Group: "case"
    Condition
      Error in `.set_test_design()`:
      ! Each subject must have at most one sample in each group for paired set testing.
