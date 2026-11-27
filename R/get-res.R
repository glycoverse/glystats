#' Get the result of a glystats analysis
#'
#' Syntax sugar for `$tidy_result` and `$raw_result` elements of a glystats result object.
#' It's useful to be used in pipes.
#'
#' @param res A glystats result object.
#' @param which Used to specify which element to get, when the result is a list.
#'   For glystats results with only one tidy result (e.g. [gly_ttest()]), this argument is not needed.
#'   For others (e.g. [gly_anova()]), this argument is required to specify which tidy result to get.
#'
#' @returns `get_tidy_result()` returns a tibble. `get_raw_result()` returns
#'   the stored raw result or its selected component, which may be a fitted
#'   model, a list of models, or `NULL` when no model was fitted.
#'
#' @examples
#' library(glyexp)
#' library(glyclean)
#' library(dplyr)
#'
#' exp <- auto_clean(real_experiment) |>
#'   glyexp::slice_head_row(n = 10)
#'
#' # Using a pipe
#' sig_res <- exp |>
#'   gly_anova() |>
#'   get_tidy_result("main_test") |>
#'   filter(p_adj < 0.05)
#'
#' # Equivalent to
#' anova_res <- gly_anova(exp)
#' sig_res <- anova_res$tidy_result$main_test |>
#'   filter(p_adj < 0.05)
#'
#' @export
get_tidy_result <- function(res, which = NULL) {
  checkmate::assert_class(res, "glystats_res")
  checkmate::assert_string(which, null.ok = TRUE)

  result <- res$tidy_result
  if (tibble::is_tibble(result)) {
    if (!is.null(which)) {
      cli::cli_abort(c(
        "This analysis has a single tidy result.",
        "i" = "Omit {.arg which} to retrieve the tibble."
      ))
    }
    return(result)
  }
  if (is.null(which)) {
    cli::cli_abort(c(
      "{.arg which} must be provided for {.cls {class(res)[[1]]}} result.",
      "i" = "Available tibbles: {.val {names(result)}}"
    ))
  }
  if (!which %in% names(result)) {
    cli::cli_abort(c(
      "{.arg which} must be one of {.val {names(result)}}"
    ))
  }
  result <- result[[which]]
  if (!tibble::is_tibble(result)) {
    cli::cli_abort("The selected tidy result must be a tibble.")
  }
  result
}

#' @rdname get_tidy_result
#' @export
get_raw_result <- function(res, which = NULL) {
  checkmate::assert_class(res, "glystats_res")
  checkmate::assert_string(which, null.ok = TRUE)
  if (is.null(which)) {
    return(res$raw_result)
  }
  if (!which %in% names(res$raw_result)) {
    cli::cli_abort(c(
      "{.arg which} must be one of {.val {names(res$raw_result)}}"
    ))
  }
  return(res$raw_result[[which]])
}
