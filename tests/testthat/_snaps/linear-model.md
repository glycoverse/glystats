# gly_linear_model validates formulas and designs

    Code
      gly_linear_model(exp, expression ~ treatment)
    Condition
      Error in `.validate_linear_model_formula()`:
      ! `formula` must be one-sided because expression values come from `exp`.

---

    Code
      gly_linear_model(exp, ~missing_predictor)
    Condition
      Error in `.build_linear_model_design()`:
      ! Failed to evaluate `formula` against sample information.
      x object 'missing_predictor' not found
      i Available columns: sample, treatment, time, age, and batch.
      Caused by error:
      ! object 'missing_predictor' not found

---

    Code
      gly_linear_model(exp, ~1)
    Condition
      Error in `.analyze_linear_model()`:
      ! The model must contain at least one non-intercept coefficient when `contrasts` is NULL.

---

    Code
      gly_linear_model(missing_exp, ~age)
    Condition
      Error in `.build_linear_model_design()`:
      ! Failed to evaluate `formula` against sample information.
      x missing values in object
      i Available columns: sample, treatment, time, age, and batch.
      Caused by error in `na.fail.default()`:
      ! missing values in object

---

    Code
      gly_linear_model(exp, ~ treatment + duplicate)
    Condition
      Error in `.build_linear_model_design()`:
      ! The model design is rank deficient.
      x Non-estimable coefficients: "duplicateB".
      i Remove redundant terms or unused factor levels from `formula`.

# gly_linear_model validates contrasts

    Code
      gly_linear_model(exp, ~treatment, contrasts = "treatmentB")
    Condition
      Error in `.validate_linear_model_contrasts()`:
      ! `contrasts` must be a named character vector.

---

    Code
      gly_linear_model(exp, ~treatment, contrasts = c(test = "treatmentB", test = "-treatmentB"))
    Condition
      Error in `.validate_linear_model_contrasts()`:
      ! `contrasts` names must be unique.

---

    Code
      gly_linear_model(exp, ~ treatment * time, contrasts = c(test = "unknown_coefficient"))
    Condition
      Error in `.analyze_linear_model()`:
      ! Failed to construct `contrasts`.
      x object 'unknown_coefficient' not found
      i Available coefficients: "(Intercept)", "treatmentB", "timeT2", and "treatmentB:timeT2".
      Caused by error:
      ! object 'unknown_coefficient' not found

---

    Code
      gly_linear_model(exp, ~ treatment * time, contrasts = c(test = "treatmentB + treatmentB:timeT2"))
    Condition
      Error in `.translate_linear_model_contrasts()`:
      ! Non-syntactic coefficient "treatmentB:timeT2" must be wrapped in backticks in `contrasts`.

# filter_sig_vars validates linear model filtering

    Code
      filter_sig_vars(exp, result, term = "unknown")
    Condition
      Error in `filter_sig_vars()`:
      ! Can't find term: "unknown".
      i Available terms: "treatmentB".

---

    Code
      filter_sig_vars(exp, result, fc_cutoff = 2)
    Condition
      Error in `filter_sig_vars()`:
      ! `fc_cutoff` cannot be used with `gly_linear_model()` results because model estimates are not necessarily fold changes.

---

    Code
      filter_sig_vars(exp, unadjusted)
    Condition
      Error in `filter_sig_vars()`:
      ! `p_adj_cutoff` cannot be used because this result has no adjusted p-values.
      i Set `p_adj_cutoff` to NULL and provide `p_val_cutoff`.
