#' Formula-based Linear Models
#'
#' Fit a linear model to every variable in a glycomics or glycoproteomics
#' experiment using limma's empirical Bayes moderation.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object,
#'   or another `SummarizedExperiment` containing expression data and sample
#'   information.
#' @param formula A one-sided formula describing the model in terms of columns
#'   in the sample information, for example `~ treatment * time + batch`.
#' @param contrasts An optional uniquely named character vector of limma-style
#'   contrast expressions. Use backticks around non-syntactic coefficient names,
#'   such as `` `treatmentB:timeT2` ``.
#' @param p_adj_method A character string specifying the method used to adjust
#'   p-values within each coefficient or contrast. See [stats::p.adjust()] for
#'   available methods. If `NULL`, no adjustment is performed.
#' @param add_info A logical value. If `TRUE` (default), variable information
#'   from the experiment is added to the tidy result.
#' @param ... Additional arguments passed to [limma::lmFit()].
#'
#' @details
#' `gly_linear_model()` and [gly_limma()] both use limma under the hood, but
#' serve complementary use cases. `gly_linear_model()` provides a general
#' formula interface for specifying a wide range of analysis designs and custom
#' contrasts. [gly_limma()] is a dedicated wrapper for common differential
#' expression analysis tasks, with convenient group, covariate, subject, and
#' pairwise-comparison arguments.
#'
#' Expression values are transformed using `log2(x + 1e-6)` before modeling.
#' When `contrasts` is `NULL`, all non-intercept coefficients are tested. When
#' contrasts are supplied, only the named contrasts are returned.
#'
#' Character predictors are handled as factors by [stats::model.matrix()].
#' Factor reference levels therefore determine the meaning of coefficients.
#' Set factor levels in the sample information before calling this function when
#' a specific reference level is required.
#'
#' P-values are adjusted independently across variables for each coefficient or
#' contrast.
#'
#' @returns A list with classes `glystats_linear_model_res` and `glystats_res`
#'   containing:
#'   - `tidy_result`: A tibble with `variable`, `term`, `estimate`, `AveExpr`,
#'     `t`, `p_val`, optional `p_adj`, and `b` columns.
#'   - `raw_result`: A list containing the moderated limma fit, design matrix,
#'     contrast matrix, and coefficient-name mapping.
#'   - `meta_data`: Metadata copied from `exp`.
#'
#' @examplesIf requireNamespace("limma", quietly = TRUE)
#' exp <- glyexp::real_experiment
#' complete <- rowSums(is.finite(SummarizedExperiment::assay(exp))) == ncol(exp)
#' exp <- exp[which(complete)[seq_len(10)], ]
#' SummarizedExperiment::colData(exp)$time <- factor(
#'   rep(c("early", "late"), length.out = ncol(exp))
#' )
#'
#' model_res <- gly_linear_model(exp, ~ group * time)
#'
#' contrast_res <- gly_linear_model(
#'   exp,
#'   ~ group * time,
#'   contrasts = c(
#'     C_late_vs_H_late = "groupC + `groupC:timelate`"
#'   )
#' )
#'
#' @seealso [gly_limma()], [limma::lmFit()], [limma::eBayes()]
#' @export
gly_linear_model <- function(
  exp,
  formula,
  contrasts = NULL,
  p_adj_method = "BH",
  add_info = TRUE,
  ...
) {
  .assert_data_container(exp)
  .validate_linear_model_formula(formula)
  .validate_linear_model_contrasts(contrasts)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  checkmate::assert_logical(add_info, len = 1)
  rlang::check_installed("limma")

  result <- .analyze_linear_model(
    expr_mat = .get_expr_mat(exp),
    sample_info = .get_sample_info(exp),
    formula = formula,
    contrasts = contrasts,
    p_adj_method = p_adj_method,
    ...
  )
  result$tidy_result$.glystats_order <- seq_len(nrow(result$tidy_result))
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  ) |>
    dplyr::arrange(.data$.glystats_order) |>
    dplyr::select(-".glystats_order")
  result$meta_data <- .get_meta_data(exp)
  result
}

.analyze_linear_model <- function(
  expr_mat,
  sample_info,
  formula,
  contrasts = NULL,
  p_adj_method = "BH",
  ...
) {
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_data_frame(sample_info, nrows = ncol(expr_mat))
  .validate_linear_model_formula(formula)
  .validate_linear_model_contrasts(contrasts)
  checkmate::assert_choice(
    p_adj_method,
    stats::p.adjust.methods,
    null.ok = TRUE
  )
  rlang::check_installed("limma")

  design_info <- .build_linear_model_design(formula, sample_info)
  design <- design_info$design
  coefficient_mapping <- design_info$coefficient_mapping

  if (is.null(contrasts)) {
    report_original <- setdiff(
      coefficient_mapping$original,
      "(Intercept)"
    )
    if (length(report_original) == 0) {
      cli::cli_abort(
        "The model must contain at least one non-intercept coefficient when {.arg contrasts} is NULL."
      )
    }
    report_internal <- coefficient_mapping$internal[
      match(report_original, coefficient_mapping$original)
    ]
    contrast_matrix <- NULL
  } else {
    translated_contrasts <- .translate_linear_model_contrasts(
      contrasts,
      coefficient_mapping
    )
    contrast_matrix <- tryCatch(
      limma::makeContrasts(
        contrasts = unname(translated_contrasts),
        levels = design
      ),
      error = function(err) {
        cli::cli_abort(
          c(
            "Failed to construct {.arg contrasts}.",
            "x" = conditionMessage(err),
            "i" = "Available coefficients: {.val {coefficient_mapping$original}}."
          ),
          parent = err
        )
      }
    )
    colnames(contrast_matrix) <- names(contrasts)
    report_original <- names(contrasts)
    report_internal <- names(contrasts)
  }

  dots <- rlang::list2(...)
  disallowed_args <- intersect(names(dots), c("object", "design"))
  if (length(disallowed_args) > 0) {
    cli::cli_abort(
      "Arguments {cli::format_inline(disallowed_args)} should not be supplied through `...`; they are controlled internally."
    )
  }

  log_expr_mat <- .log_transform_expr_mat(expr_mat)
  row_has_data <- rowSums(is.finite(log_expr_mat)) > 0
  if (!any(row_has_data)) {
    cli::cli_abort(
      "No finite expression values remain after log2 transformation."
    )
  }
  model_expr_mat <- log_expr_mat[row_has_data, , drop = FALSE]
  dots <- .subset_linear_model_weights(
    dots,
    row_has_data = row_has_data,
    expression_dim = dim(log_expr_mat)
  )

  fit <- rlang::exec(
    limma::lmFit,
    object = model_expr_mat,
    design = design,
    !!!dots
  )
  if (!is.null(contrast_matrix)) {
    fit <- limma::contrasts.fit(fit, contrast_matrix)
  }
  fit <- .moderate_linear_model_fit(fit)

  tidy_result <- .tidy_linear_model_fit(
    fit = fit,
    variable_names = rownames(expr_mat),
    term_names = report_original,
    coefficient_names = report_internal,
    p_adj_method = p_adj_method
  )

  structure(
    list(
      tidy_result = tidy_result,
      raw_result = list(
        fit = fit,
        design = design,
        contrast_matrix = contrast_matrix,
        coefficient_mapping = coefficient_mapping
      )
    ),
    class = c("glystats_linear_model_res", "glystats_res")
  )
}

.validate_linear_model_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort("{.arg formula} must be a formula.")
  }
  if (length(formula) != 2L) {
    cli::cli_abort(
      "{.arg formula} must be one-sided because expression values come from {.arg exp}."
    )
  }
  invisible(NULL)
}

.validate_linear_model_contrasts <- function(contrasts) {
  checkmate::assert_character(contrasts, min.len = 1, null.ok = TRUE)
  if (is.null(contrasts)) {
    return(invisible(NULL))
  }

  contrast_names <- names(contrasts)
  if (
    is.null(contrast_names) ||
      any(is.na(contrast_names)) ||
      any(contrast_names == "")
  ) {
    cli::cli_abort("{.arg contrasts} must be a named character vector.")
  }
  if (anyDuplicated(contrast_names) > 0) {
    cli::cli_abort("{.arg contrasts} names must be unique.")
  }
  invisible(NULL)
}

.build_linear_model_design <- function(formula, sample_info) {
  model_frame <- tryCatch(
    stats::model.frame(
      formula,
      data = as.data.frame(sample_info),
      na.action = stats::na.fail,
      drop.unused.levels = TRUE
    ),
    error = function(err) {
      cli::cli_abort(
        c(
          "Failed to evaluate {.arg formula} against sample information.",
          "x" = conditionMessage(err),
          "i" = "Available columns: {.field {colnames(sample_info)}}."
        ),
        parent = err
      )
    }
  )
  design <- stats::model.matrix(formula, data = model_frame)

  if (qr(design)$rank < ncol(design)) {
    non_estimable <- limma::nonEstimable(design)
    cli::cli_abort(c(
      "The model design is rank deficient.",
      "x" = "Non-estimable coefficients: {.val {non_estimable}}.",
      "i" = "Remove redundant terms or unused factor levels from {.arg formula}."
    ))
  }

  original_names <- colnames(design)
  internal_names <- make.names(original_names, unique = TRUE)
  colnames(design) <- internal_names
  coefficient_mapping <- tibble::tibble(
    original = original_names,
    internal = internal_names
  )

  list(
    design = design,
    coefficient_mapping = coefficient_mapping
  )
}

.translate_linear_model_contrasts <- function(contrasts, mapping) {
  translated_contrasts <- contrasts
  for (contrast_index in seq_along(contrasts)) {
    original_contrast <- contrasts[[contrast_index]]
    translated <- original_contrast
    for (i in seq_len(nrow(mapping))) {
      original <- mapping$original[[i]]
      internal <- mapping$internal[[i]]
      if (!identical(original, internal)) {
        backticked <- paste0("`", original, "`")
        if (
          .has_unquoted_linear_model_coefficient(
            original_contrast,
            original
          )
        ) {
          cli::cli_abort(
            "Non-syntactic coefficient {.val {original}} must be wrapped in backticks in {.arg contrasts}."
          )
        }
        translated <- gsub(
          backticked,
          internal,
          translated,
          fixed = TRUE
        )
      }
    }
    translated_contrasts[[contrast_index]] <- translated
  }
  translated_contrasts
}

.has_unquoted_linear_model_coefficient <- function(expression, coefficient) {
  unquoted_expression <- gsub(
    "`[^`]*`",
    "",
    expression,
    perl = TRUE
  )
  positions <- gregexpr(
    coefficient,
    unquoted_expression,
    fixed = TRUE
  )[[1]]
  if (identical(positions[[1]], -1L)) {
    return(FALSE)
  }

  coefficient_length <- nchar(coefficient)
  expression_length <- nchar(unquoted_expression)
  any(vapply(
    positions,
    function(position) {
      before <- if (position > 1L) {
        substr(unquoted_expression, position - 1L, position - 1L)
      } else {
        ""
      }
      after_position <- position + coefficient_length
      after <- if (after_position <= expression_length) {
        substr(unquoted_expression, after_position, after_position)
      } else {
        ""
      }
      !grepl("[[:alnum:]_.:]", before) &&
        !grepl("[[:alnum:]_.:]", after)
    },
    logical(1)
  ))
}

.subset_linear_model_weights <- function(
  dots,
  row_has_data,
  expression_dim
) {
  weights <- dots$weights
  if (is.null(weights)) {
    return(dots)
  }

  if (is.matrix(weights)) {
    if (
      nrow(weights) == expression_dim[[1]] &&
        ncol(weights) %in% c(1L, expression_dim[[2]])
    ) {
      dots$weights <- weights[row_has_data, , drop = FALSE]
    }
  } else if (is.numeric(weights) && length(weights) == expression_dim[[1]]) {
    dots$weights <- weights[row_has_data]
  }

  dots
}

.moderate_linear_model_fit <- function(fit) {
  tryCatch(
    limma::eBayes(fit, trend = TRUE),
    error = function(err) {
      msg <- conditionMessage(err)
      if (grepl("covariate", msg, ignore.case = TRUE)) {
        cli::cli_warn(
          "Limma trend fitting failed due to covariate issues; re-running eBayes with trend = FALSE."
        )
        limma::eBayes(fit, trend = FALSE)
      } else {
        stop(err)
      }
    }
  )
}

.tidy_linear_model_fit <- function(
  fit,
  variable_names,
  term_names,
  coefficient_names,
  p_adj_method
) {
  purrr::map2_dfr(
    term_names,
    coefficient_names,
    function(term, coefficient) {
      top_table <- limma::topTable(
        fit,
        coef = coefficient,
        number = Inf,
        adjust.method = if (is.null(p_adj_method)) "none" else p_adj_method,
        sort.by = "none"
      )
      result <- tibble::as_tibble(
        top_table,
        rownames = "variable"
      ) |>
        dplyr::rename(
          estimate = "logFC",
          p_val = "P.Value",
          p_adj = "adj.P.Val",
          b = "B"
        )

      if (is.null(p_adj_method)) {
        result <- dplyr::select(result, -"p_adj")
      }

      result <- tibble::tibble(variable = variable_names) |>
        dplyr::left_join(result, by = "variable") |>
        dplyr::mutate(term = .env$term, .after = "variable")

      result_columns <- c(
        "variable",
        "term",
        "estimate",
        "AveExpr",
        "t",
        "p_val"
      )
      if (!is.null(p_adj_method)) {
        result_columns <- c(result_columns, "p_adj")
      }
      result_columns <- c(result_columns, "b")
      dplyr::select(result, dplyr::all_of(result_columns))
    }
  )
}
