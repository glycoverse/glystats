#' Get the result of a glystats analysis
#'
#' Syntax sugar for `$tidy_result` and `$raw_result` elements of a glystats result object.
#' It's useful to be used in pipes.
#'
#' @param res A glystats result object.
#' @param which Used to specify which element to get, when the result is a list.
#'   If NULL, the whole result (tibble or list) will be returned.
#'
#' @returns A tibble.
#'
#' @examples
#' library(glyexp)
#' library(glyclean)
#' library(dplyr)
#'
#' exp <- auto_clean(real_experiment) |>
#'   glyexp::slice_head_var(n = 10)
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
  if (is.null(which)) {
    return(res$tidy_result)
  }
  if (!which %in% names(res$tidy_result)) {
    cli::cli_abort(c(
      "{.arg which} must be one of {.val {names(res$tidy_result)}}"
    ))
  }
  return(res$tidy_result[[which]])
}

#' @rdname get_tidy_result
#' @export
get_raw_result <- function(res, which = NULL) {
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