#' Abundance scales and automatic transformations
#'
#' @name gly_scale
#' @param log_transform `NULL` (default) uses `metadata(exp)$scale` to avoid
#'   logging data twice. `TRUE` forces the method's log transformation on the
#'   supplied values, even if they are already logged; `FALSE` uses the supplied
#'   values unchanged. This argument must be named when it follows `...`.
#'
#' @details
#' The supported scales are `"raw"`, `"log"` (natural logarithm), `"log2"`,
#' and `"log10"`. Missing scale metadata is treated as `"raw"` for compatibility
#' with older experiments. An explicitly supplied invalid scale is an error.
#' Older development objects labeled `"log"` must be relabeled if the data
#' were actually transformed with base 2 or 10. The base is never inferred.
#'
#' With `log_transform = NULL`, known logged data are used unchanged. Raw data
#' retain the method's historical default: PCA, t-SNE, UMAP, PLS-DA, and OPLS-DA
#' use `log(x + 1)`; differential tests, linear models, set tests, correlation,
#' and clustering use `log2(x + 1e-6)`. ROC, Cox, and fold change use raw values.
#' For these last three methods, `TRUE` requests `log2(x + 1e-6)`.
#'
#' Results retain the input metadata in `meta_data` and describe the actual
#' analysis in `preprocessing`: `input_scale`, `scale_missing`,
#' `requested_log_transform`, `log_transform` (whether a new log was applied),
#' `analysis_scale`, and `transformation`. The input object is never modified.
#' Model estimates, standard errors, intervals, and raw fits use analysis units.
#'
#' For raw inputs, fold changes retain `log2(mean(test) / mean(reference))`.
#' For logged inputs, `log2fc` is the difference of group means converted to
#' base 2: divide natural-log differences by `log(2)` and multiply base-10
#' differences by `log2(10)`. This estimates a ratio of geometric means.
#' Differential tests calculate fold changes from the input scale independently
#' of the matrix used for testing. In `gly_fold_change(log_transform = TRUE)`,
#' fold changes instead use the newly logged matrix, including its pseudo-count.
#' Limma contrasts on a logged analysis scale are converted to `log2fc`;
#' native model estimates and fits stay in analysis units. With an unlogged raw
#' limma analysis, contrasts are returned as `estimate` and `log2fc` is `NA`,
#' because an adjusted raw-scale difference is not a log fold change.
#'
#' @examples
#' values <- matrix(c(1, 2, 4, 8), nrow = 1,
#'                  dimnames = list("feature", paste0("s", 1:4)))
#' exp <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(abundance = log10(values)),
#'   colData = S4Vectors::DataFrame(group = c("A", "A", "B", "B")),
#'   metadata = list(scale = "log10")
#' )
#' res <- gly_fold_change(exp)
#' res$tidy_result
#' res$preprocessing
NULL

.resolve_transformation <- function(
  exp,
  log_transform,
  default = TRUE,
  method = "log2",
  meta_data = .get_meta_data(exp)
) {
  checkmate::assert_logical(
    log_transform,
    len = 1,
    any.missing = FALSE,
    null.ok = TRUE
  )
  meta <- meta_data
  scale_missing <- !"scale" %in% names(meta)
  input_scale <- if (scale_missing) "raw" else meta[["scale"]]
  if (!checkmate::test_choice(input_scale, c("raw", "log", "log2", "log10"))) {
    cli::cli_abort('Metadata `scale` must be "raw", "log", "log2", or "log10".')
  }
  apply_log <- if (is.null(log_transform)) {
    identical(input_scale, "raw") && default
  } else {
    log_transform
  }
  list(
    input_scale = input_scale,
    scale_missing = scale_missing,
    requested_log_transform = log_transform,
    log_transform = apply_log,
    analysis_scale = if (apply_log) method else input_scale,
    transformation = if (!apply_log) {
      "none"
    } else if (method == "log") {
      "log(x + 1)"
    } else {
      "log2(x + 1e-6)"
    }
  )
}


.log2_factor <- function(scale) {
  switch(scale, log = 1 / log(2), log2 = 1, log10 = log2(10))
}

.mean_log2fc <- function(reference, test, input_scale = "raw") {
  if (input_scale == "raw") {
    log2(test / reference)
  } else {
    (test - reference) * .log2_factor(input_scale)
  }
}

.limma_effect_scale <- function(result, transformation) {
  if (is.null(transformation) || transformation$analysis_scale == "log2") {
    return(result)
  }
  result$tidy_result$estimate <- result$tidy_result$log2fc
  result$tidy_result$log2fc <- if (transformation$analysis_scale == "raw") {
    rep(NA_real_, nrow(result$tidy_result))
  } else {
    result$tidy_result$estimate * .log2_factor(transformation$analysis_scale)
  }
  result
}
