#' Consensus Clustering for Glycomics and Glycoproteomics Data
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated because we realized that the interactive nature of consensus clustering
#' is not ideal for a programmatic API of glystats.
#' Users are encouraged to use the `ConsensusClusterPlus` package directly.
#'
#' Perform consensus clustering on the expression data using ConsensusClusterPlus.
#' The function returns cluster assignments for all k values from 2 to max_k in long format.
#' Setting `output_file` to visualize consensus CDF curves and consensus matrices,
#' which helps you decide the optimal number of clusters.
#' See [this tutorial](https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf)
#' for more information.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat (Only for [gly_consensus_clustering_()]) A numeric matrix with variables as rows and samples as columns.
#' @param on A character string specifying what to cluster. Either "sample" (default) to cluster
#'   samples/observations, or "variable" to cluster variables/features.
#' @param max_k Maximum number of clusters to test. Default is 9.
#'   Passed to the `maxK` argument of `ConsensusClusterPlus::ConsensusClusterPlus()`.
#' @param reps Number of resampling iterations. Default is 1000.
#'   Passed to the `reps` argument of `ConsensusClusterPlus::ConsensusClusterPlus()`.
#' @param p_item Proportion of items to sample. Default is 0.8.
#'   Passed to the `pItem` argument of `ConsensusClusterPlus::ConsensusClusterPlus()`.
#' @param cluster_alg Clustering algorithm to use. Default is "hc" (hierarchical clustering).
#'   Other options include "km" (k-means), "pam" (partitioning around medoids).
#'   Passed to the `clusterAlg` argument of `ConsensusClusterPlus::ConsensusClusterPlus()`.
#' @param scale A logical indicating whether to scale the data. Default is TRUE.
#' @param output_file A character string specifying the output file path for plots.
#'   If NULL (default), no plots are saved and no files are generated. If provided,
#'   the consensus clustering plot will be saved directly to this file path.
#'   The function will create any necessary parent directories.
#' @param add_info A logical value. If TRUE (default), sample and variable information from the experiment
#'  will be added to the result tibbles. If FALSE, only the consensus clustering results are returned.
#'  Only applicable to `gly_consensus_clustering()`.
#' @param ... Additional arguments passed to `ConsensusClusterPlus::ConsensusClusterPlus()`.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `ConsensusClusterPlus` for consensus clustering
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' consensus clustering. When `on = "sample"` (default), samples are clustered based
#' on their expression profiles across variables. When `on = "variable"`, variables are clustered based on their
#' expression patterns across samples.
#'
#' `gly_consensus_clustering()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_consensus_clustering_()` is the underlying API that works with matrices directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' **Consensus Clustering Process:**
#' 1. Perform consensus clustering for k = 2 to max_k using ConsensusClusterPlus
#' 2. Return cluster assignments for all k values in long format
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with consensus clustering results containing the following columns:
#'    - `variable` or `sample`: Variable or sample name (depending on `on` parameter)
#'    - `k`: Number of clusters
#'    - `cluster`: Cluster assignment for the corresponding k
#'  - `raw_result`: The raw ConsensusClusterPlus object
#'  - `meta_data` (only for [gly_consensus_clustering()]): A list containing metadata from the input experiment
#' The list has classes `glystats_cc_res` and `glystats_res`.
#'
#' @seealso [ConsensusClusterPlus::ConsensusClusterPlus()]
#' @keywords internal
#' @export
gly_consensus_clustering <- function(
  exp,
  on = "sample",
  max_k = 9,
  reps = 1000,
  p_item = 0.8,
  cluster_alg = "hc",
  scale = TRUE,
  output_file = NULL,
  add_info = TRUE,
  ...
) {
  lifecycle::deprecate_warn(
    when = "0.7.0",
    what = "gly_consensus_clustering()",
    details = "Please use the `ConsensusClusterPlus` package directly for consensus clustering analysis."
  )
  checkmate::assert_class(exp, "glyexp_experiment")

  expr_mat <- glyexp::get_expr_mat(exp)
  result <- gly_consensus_clustering_(
    expr_mat,
    on,
    max_k,
    reps,
    p_item,
    cluster_alg,
    scale,
    output_file,
    ...
  )
  result <- .process_results_add_info(result, exp, add_info)

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' @rdname gly_consensus_clustering
#' @export
gly_consensus_clustering_ <- function(
  expr_mat,
  on = "sample",
  max_k = 9,
  reps = 1000,
  p_item = 0.8,
  cluster_alg = "hc",
  scale = TRUE,
  output_file = NULL,
  ...
) {
  lifecycle::deprecate_warn(
    when = "0.7.0",
    what = "gly_consensus_clustering_()",
    details = "This function will be removed in a future version."
  )
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_choice(on, c("sample", "variable"))
  checkmate::assert_integerish(max_k, lower = 2, len = 1)
  checkmate::assert_integerish(reps, lower = 1, len = 1)
  checkmate::assert_number(p_item, lower = 0, upper = 1)
  checkmate::assert_choice(cluster_alg, c("hc", "km", "pam"))
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_character(output_file, len = 1, null.ok = TRUE)

  # Check if required packages are available
  rlang::check_installed("ConsensusClusterPlus")

  # Prepare data for clustering based on 'on' parameter
  if (on == "sample") {
    # Cluster samples: samples as rows, variables as columns
    mat <- t(expr_mat)
  } else {
    # Cluster variables: variables as rows, samples as columns
    mat <- expr_mat
  }

  # Apply log transformation
  mat <- log2(mat + 1)
  # Scale data if requested
  if (scale) {
    mat <- scale(mat)
  }

  # Perform consensus clustering
  # Suppress output from ConsensusClusterPlus
  withr::with_options(
    list(warn = -1),
    {
      if (is.null(output_file)) {
        # When no output file is specified, disable plotting completely
        # Use a temporary directory to avoid any file generation
        cc_res <- .cc_no_output(
          mat,
          maxK = max_k,
          reps = reps,
          pItem = p_item,
          clusterAlg = cluster_alg,
          ...
        )
      } else {
        # When output file is specified, use a temporary directory and then move the result
        cc_res <- .cc_with_output(
          mat,
          output_file,
          maxK = max_k,
          reps = reps,
          pItem = p_item,
          clusterAlg = cluster_alg,
          ...
        )
      }
    }
  )

  # Extract cluster assignments for all k values
  cluster_results <- purrr::map_dfr(2:max_k, function(k) {
    clusters <- cc_res[[k]]$consensusClass
    tibble::tibble(
      item = names(clusters),
      k = k,
      cluster = as.integer(clusters)
    )
  })

  # Set the appropriate column name based on clustering type
  colnames(cluster_results)[1] <- on

  # Add S3 class to tidy result
  tidy_result <- structure(
    cluster_results,
    class = c("glystats_cc_res", "glystats_res", class(cluster_results))
  )

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = cc_res
    ),
    class = c("glystats_cc_res", "glystats_res")
  )
}

.cc_no_output <- function(mat, ...) {
  temp_dir <- withr::local_tempdir("consensus_temp")

  dots <- rlang::list2(...)
  if ("d" %in% names(dots)) {
    cli::cli_abort(
      "{.field d} should not be supplied through `...`; data comes from the function inputs."
    )
  }

  # Capture and suppress all graphics output
  suppressMessages({
    withr::with_dir(temp_dir, {
      # Temporarily redirect graphics to a null device
      temp_pdf <- tempfile(fileext = ".pdf")
      withr::with_pdf(temp_pdf, {
        call_args <- c(
          list(d = t(mat)),
          dots
        )
        if (!"title" %in% names(call_args)) {
          call_args$title <- "temp_consensus"
        }
        if (!"plot" %in% names(call_args)) {
          call_args$plot <- "pdf"
        }
        if (!"verbose" %in% names(call_args)) {
          call_args$verbose <- FALSE
        }
        do.call(ConsensusClusterPlus::ConsensusClusterPlus, call_args)
      })
    })
  })
}

.cc_with_output <- function(mat, output_file, ...) {
  temp_dir <- withr::local_tempdir("consensus_temp")
  dots <- rlang::list2(...)
  if ("d" %in% names(dots)) {
    cli::cli_abort(
      "{.field d} should not be supplied through `...`; data comes from the function inputs."
    )
  }
  cc_res <- suppressMessages({
    call_args <- c(
      list(d = t(mat)),
      dots
    )
    if (!"title" %in% names(call_args)) {
      call_args$title <- file.path(temp_dir, "consensus_output")
    }
    if (!"plot" %in% names(call_args)) {
      call_args$plot <- "pdf"
    }
    if (!"verbose" %in% names(call_args)) {
      call_args$verbose <- FALSE
    }
    do.call(ConsensusClusterPlus::ConsensusClusterPlus, call_args)
  })

  # Move the generated consensus.pdf to the specified output file
  consensus_pdf <- file.path(temp_dir, "consensus_output", "consensus.pdf")
  if (file.exists(consensus_pdf)) {
    # Ensure the output directory exists
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir) && output_dir != ".") {
      dir.create(output_dir, recursive = TRUE)
    }
    file.copy(consensus_pdf, output_file, overwrite = TRUE)
  }

  cc_res
}
