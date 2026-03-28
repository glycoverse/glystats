#' Weighted Gene Co-expression Network Analysis (WGCNA)
#'
#' `r lifecycle::badge("deprecated")`
#'
#' Perform WGCNA analysis to identify co-expression modules and their eigenvalues.
#' The function uses WGCNA package to construct weighted gene co-expression networks,
#' detect modules, and calculate module membership and eigenvalues.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param expr_mat (Only for [gly_wgcna_()]) A numeric matrix with variables as rows and samples as columns.
#' @param powers A numeric vector of soft thresholding powers to test.
#'   Default is `c(1:10, seq(12, 20, by = 2))`.
#' @param network_type Character string specifying network type.
#'   Allowed values are "unsigned", "signed", "signed hybrid".
#'   Default is "unsigned". Passed to the `NetworkType`` argument of `WGCNA::blockwiseModules()`.
#' @param tom_type Character string specifying topological overlap matrix type.
#'   Allowed values are "none", "unsigned", "signed".
#'   Default is "unsigned". Passed to the `TOMType`` argument of `WGCNA::blockwiseModules()`.
#' @param min_module_size Minimum module size for module detection.
#'   Default is 30. Passed to the `minModuleSize`` argument of `WGCNA::blockwiseModules()`.
#' @param deep_split Integer value between 0 and 4. Provides a simplified control
#'   over how sensitive module detection should be to module splitting.
#'   Default is 2. Passed to the `deepSplit`` argument of `WGCNA::blockwiseModules()`.
#' @param merge_cut_height Dendrogram cut height for merging of modules.
#'   Default is 0.25. Passed to the `mergeCutHeight`` argument of `WGCNA::blockwiseModules()`.
#' @param add_info A logical value indicating whether to add variable information
#'   to the modules tibble. Only applicable to `gly_wgcna()`.
#' @param ... Additional arguments passed to `WGCNA::blockwiseModules()`.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - WGCNA` for weighted gene co-expression network analysis
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' WGCNA analysis.
#'
#' `gly_wgcna()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_wgcna_()` is the underlying API that works with matrices directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' **Analysis Steps:**
#' 1. Soft threshold selection using `WGCNA::pickSoftThreshold()`
#' 2. Network construction and module detection using `WGCNA::blockwiseModules()`
#' 3. Module membership calculation based on correlation with module eigengenes
#' 4. Results formatting into standardized tibbles
#'
#' @return A list with two elements:
#' - `tidy_result`: A list containing two tibbles:
#'   - `modules`: Module assignments and membership values containing the following columns:
#'     - `variable`: Variable name
#'     - `module`: Module assignment (color name)
#'     - `membership`: Module membership value (correlation with module eigengene)
#'   - `eigenvalues`: Module eigenvalues containing the following columns:
#'     - `module`: Module name (color name)
#'     - `sample`: Sample name
#'     - `eigenvalue`: Module eigenvalue (first principal component of module expression)
#' - `raw_result`: The raw WGCNA blockwiseModules object
#' The list has classes `glystats_wgcna_res` and `glystats_res`.
#'
#' @seealso [WGCNA::pickSoftThreshold()], [WGCNA::blockwiseModules()]
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
gly_wgcna <- function(
  exp,
  powers = c(1:10, seq(12, 20, by = 2)),
  network_type = "unsigned",
  tom_type = "unsigned",
  min_module_size = 30,
  deep_split = 2,
  merge_cut_height = 0.25,
  add_info = TRUE,
  ...
) {
  lifecycle::deprecate_warn(
    when = "0.7.0",
    what = "gly_wgcna()",
    details = "This function will be removed in a future version."
  )
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  # Other argments are validated in gly_wgcna_().

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)

  # Call the underlying API
  result <- gly_wgcna_(
    expr_mat,
    powers,
    network_type,
    tom_type,
    min_module_size,
    deep_split,
    merge_cut_height,
    ...
  )

  # Process results with add_info logic
  result$tidy_result <- .process_results_add_info(
    result$tidy_result,
    exp,
    add_info
  )
  result
}

#' @rdname gly_wgcna
#' @export
gly_wgcna_ <- function(
  expr_mat,
  powers = c(1:10, seq(12, 20, by = 2)),
  network_type = "unsigned",
  tom_type = "unsigned",
  min_module_size = 30,
  deep_split = 2,
  merge_cut_height = 0.25,
  ...
) {
  lifecycle::deprecate_warn(
    when = "0.7.0",
    what = "gly_wgcna_()",
    details = "This function will be removed in a future version."
  )
  # Check if WGCNA package is available
  rlang::check_installed("WGCNA")

  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_numeric(powers, min.len = 1)
  checkmate::assert_choice(
    network_type,
    c("unsigned", "signed", "signed hybrid")
  )
  checkmate::assert_choice(tom_type, c("none", "unsigned", "signed"))
  checkmate::assert_int(min_module_size, lower = 1)
  checkmate::assert_int(deep_split, lower = 0, upper = 4)
  checkmate::assert_number(merge_cut_height, lower = 0, upper = 1)

  dots <- rlang::list2(...)
  if ("datExpr" %in% names(dots)) {
    cli::cli_abort(
      "{.field datExpr} should not be supplied through `...`; data comes from the function inputs."
    )
  }

  # Prepare data for WGCNA (samples as rows, variables as columns)
  # Apply log2 transformation
  expr_mat <- log2(t(expr_mat) + 1)

  if (!"power" %in% names(dots)) {
    # Step 1: Choose soft threshold power
    # Suppress all output from pickSoftThreshold (including direct console output)
    threshold_res <- suppressMessages(suppressWarnings({
      utils::capture.output(
        {
          threshold_res_temp <- WGCNA::pickSoftThreshold(
            expr_mat,
            powerVector = powers,
            verbose = 0
          )
        },
        type = "output"
      )
      threshold_res_temp
    }))

    # Extract the recommended power
    power <- threshold_res$powerEstimate
    if (is.na(power)) {
      cli::cli_warn(c(
        "Can't determine optimal soft thresholding power.",
        "x" = "This implies the data is not suitable for WGCNA.",
        "i" = "Falling back to the first power that gives R^2 > 0.8 or the power with highest R^2."
      ))
      # If no power is recommended, use the first power that gives R^2 > 0.8
      sft_data <- threshold_res$fitIndices
      power_candidates <- sft_data$Power[sft_data$SFT.R.sq > 0.8]
      if (length(power_candidates) > 0) {
        power <- min(power_candidates)
      } else {
        # If still no suitable power, use the power with highest R^2
        power <- sft_data$Power[which.max(sft_data$SFT.R.sq)]
      }
    }

    # Inform user about the selected power
    cli::cli_alert_info("Selected soft threshold power: {.val {power}}")
  } else {
    power <- dots$power
    cli::cli_alert_info(
      "Using user-supplied soft threshold power: {.val {power}}"
    )
  }

  # Step 2: Network construction and module detection
  # **Note about `cor`**
  # In source code of `WGCNA`, they defined their own version of `cor()`.
  # Latter in `WGCNA::blockwiseModules()`,
  # they used `corFun <- match.fun("cor")` to fetch the correlation function.
  # `match.fun()` tries to find `cor()` in the `parent.frame()`.
  # In an interactive R console,
  # `library(WGCNA)` will attatch "package:WGCNA" in the searching path,
  # so `match.fun("cor")` returns `WGCNA::cor`.
  # However, in our case, we call `WGCNA::blockwiseModules()` directly.
  # `match.fun("cor")` can only find `stats::cor()` in this case.
  # We blame `WGCNA` for not using `::` to specify the package name.
  # To resolve this, we first extract `WGCNA::cor()` from the package namespace.
  # Then, we use `local()` to temporarily assign `cor` to `WGCNA::cor` in the function scope.
  wgcna_cor <- utils::getFromNamespace("cor", "WGCNA")
  net <- local({
    cor <- wgcna_cor
    suppressMessages(suppressWarnings(
      {
        call_args <- c(
          list(datExpr = expr_mat),
          dots
        )
        if (!"power" %in% names(call_args)) {
          call_args$power <- power
        }
        if (!"networkType" %in% names(call_args)) {
          call_args$networkType <- network_type
        }
        if (!"TOMType" %in% names(call_args)) {
          call_args$TOMType <- tom_type
        }
        if (!"minModuleSize" %in% names(call_args)) {
          call_args$minModuleSize <- min_module_size
        }
        if (!"deepSplit" %in% names(call_args)) {
          call_args$deepSplit <- deep_split
        }
        if (!"mergeCutHeight" %in% names(call_args)) {
          call_args$mergeCutHeight <- merge_cut_height
        }
        if (!"verbose" %in% names(call_args)) {
          call_args$verbose <- 0
        }
        do.call(WGCNA::blockwiseModules, call_args)
      }
    ))
  })

  # Step 3: Process results
  tidy_result <- .process_wgcna_results(net, expr_mat, rownames(expr_mat))

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = net
    ),
    class = c("glystats_wgcna_res", "glystats_res")
  )
}

# Internal helper function to process WGCNA results
.process_wgcna_results <- function(net, expr_mat, variable_names) {
  # Extract modules and eigenvalues
  modules_df <- tibble::tibble(
    variable = names(net$colors),
    module = net$colors
  )
  eigenvalue_df <- net$MEs %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      -all_of("sample"),
      names_to = "module",
      values_to = "eigenvalue",
      names_prefix = "ME"
    )

  # Calculate memberships
  membership_df <- stats::cor(
    expr_mat,
    net$MEs,
    use = "pairwise.complete.obs"
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("variable") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      -all_of("variable"),
      names_to = "module",
      values_to = "membership",
      names_prefix = "ME"
    )
  modules_df <- modules_df %>%
    dplyr::left_join(membership_df, by = c("variable", "module"))

  # Create result list
  list(
    "modules" = modules_df,
    "eigenvalues" = eigenvalue_df
  )
}
