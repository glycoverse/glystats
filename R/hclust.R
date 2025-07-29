#' Hierarchical Clustering for Glycomics and Glycoproteomics Data
#'
#' Perform hierarchical clustering on the expression data.
#' The function uses `stats::hclust()` to perform clustering and provides
#' tidy results including cluster assignments, dendrogram data for plotting,
#' and merge heights.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param on A character string specifying what to cluster. Either "variable" (default) to cluster
#'   variables/features, or "sample" to cluster samples/observations.
#' @param k_values A numeric vector specifying the number of clusters to cut the tree into.
#'   Default is c(2, 3, 4, 5). If NULL, no cluster assignments are returned.
#' @param scale A logical indicating whether to scale the data before clustering. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'   will be added to the result tibbles. If FALSE, only the clustering results are returned.
#'   Only applicable to `gly_hclust()`.
#' @param ... Additional arguments passed to `stats::dist()` and `stats::hclust()`.
#'   Note: if both functions need a `method` parameter, use `dist.method` for distance
#'   and `hclust.method` for clustering method.
#'
#' @section Required packages:
#' This function only uses base R packages and does not require additional dependencies.
#' For enhanced dendrogram plotting capabilities, the `ggdendro` package is recommended
#' but not required.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' clustering. When `on = "variable"` (default), variables are clustered based on their
#' expression patterns across samples. When `on = "sample"`, samples are clustered based
#' on their expression profiles across variables.
#'
#' `gly_hclust()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_hclust_()` is the underlying API that works with matrices directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' **Distance Calculation:**
#' Distance is calculated using `stats::dist()` with the specified method.
#'
#' **Clustering Method:**
#' Hierarchical clustering is performed using `stats::hclust()` with the specified method.
#'
#' **Cluster Assignment:**
#' The dendrogram is cut at different heights to produce cluster assignments for
#' the specified k values using `stats::cutree()`.
#'
#' @return A list containing:
#'  - `tidy_result`: A list of tibbles with clustering results:
#'    - `clusters`: Cluster assignments containing the following columns:
#'      - `variable` or `sample`: Variable or sample name (depending on `on` parameter)
#'      - `cluster_k2`, `cluster_k3`, etc.: Cluster assignments for different k values
#'    - `dendrogram`: Dendrogram segment data for plotting (if ggdendro is available) containing:
#'      - `x`, `y`, `xend`, `yend`: Segment coordinates for plotting
#'    - `heights`: Merge heights and steps containing the following columns:
#'      - `merge_step`: Step number in the clustering process
#'      - `height`: Height at which clusters are merged
#'      - `n_clusters`: Number of clusters remaining after this merge
#'    - `labels`: Labels and their positions (if ggdendro is available) containing:
#'      - `x`, `y`: Position coordinates
#'      - `label`: Label text
#'  - `raw_result`: The raw hclust object from `stats::hclust()`
#'
#' @seealso [stats::hclust()], [stats::dist()], [stats::cutree()]
#' @export
gly_hclust <- function(
  exp,
  on = "variable",
  k_values = c(2, 3, 4, 5),
  scale = TRUE,
  add_info = TRUE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_choice(on, c("variable", "sample"))
  if (!is.null(k_values)) {
    checkmate::assert_integerish(k_values, lower = 2, min.len = 1)
  }
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)

  # Call the underlying API
  result <- gly_hclust_(expr_mat, on, k_values, scale, ...)
  result$tidy_result <- .process_results_add_info(result$tidy_result, exp, add_info)
  result
}

#' @rdname gly_hclust
#' @export
gly_hclust_ <- function(
  expr_mat,
  on = "variable",
  k_values = c(2, 3, 4, 5),
  scale = TRUE,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_choice(on, c("variable", "sample"))
  if (!is.null(k_values)) {
    checkmate::assert_integerish(k_values, lower = 2, min.len = 1)
  }
  checkmate::assert_logical(scale, len = 1)

  # Prepare data for clustering based on 'on' parameter
  if (on == "sample") {
    # Cluster samples: samples as rows, variables as columns
    mat <- t(expr_mat)
    cluster_type <- "sample"
  } else {
    # Cluster variables: variables as rows, samples as columns
    mat <- expr_mat
    cluster_type <- "variable"
  }

  # Apply log transformation
  mat <- log2(mat + 1)
  # Scale data if requested
  if (scale) {
    mat <- scale(mat)
  }

  # Extract arguments for dist() and hclust()
  dots <- list(...)

  # Handle method parameter conflicts
  dist_args <- dots
  hclust_args <- dots

  # If dist.method is specified, use it for dist() and remove from hclust_args
  if ("dist.method" %in% names(dots)) {
    dist_args$method <- dots$dist.method
    dist_args$dist.method <- NULL
    hclust_args$dist.method <- NULL
  }

  # If hclust.method is specified, use it for hclust() and remove from dist_args
  if ("hclust.method" %in% names(dots)) {
    hclust_args$method <- dots$hclust.method
    hclust_args$hclust.method <- NULL
    dist_args$hclust.method <- NULL
  }

  # Calculate distance matrix
  dist_mat <- do.call(stats::dist, c(list(x = mat), dist_args))

  # Perform hierarchical clustering
  hclust_res <- do.call(stats::hclust, c(list(d = dist_mat), hclust_args))

  # Initialize result list
  tidy_result <- list()

  # 1. Cluster assignments for different k values
  if (!is.null(k_values)) {
    cluster_assignments <- purrr::map_dfc(k_values, function(k) {
      clusters <- stats::cutree(hclust_res, k = k)
      col_name <- paste0("cluster_k", k)
      result_tbl <- tibble::tibble(x = clusters)
      colnames(result_tbl) <- col_name
      result_tbl
    })

    # Create clusters tibble with appropriate column name
    clusters_tbl <- tibble::tibble(x = rownames(mat)) %>%
      dplyr::bind_cols(cluster_assignments)

    # Set the appropriate column name based on clustering type
    colnames(clusters_tbl)[1] <- cluster_type
    tidy_result$clusters <- clusters_tbl
  }

  # 2. Heights data
  tidy_result$heights <- tibble::tibble(
    merge_step = seq_along(hclust_res$height),
    height = hclust_res$height,
    n_clusters = length(hclust_res$labels) - seq_along(hclust_res$height) + 1
  )

  # 3. Try to create dendrogram data using ggdendro if available
  if (requireNamespace("ggdendro", quietly = TRUE)) {
    tryCatch({
      # Prevent opening graphics device window on macOS using withr
      withr::with_options(
        list(device = function() grDevices::pdf(NULL)),
        dendro_data <- ggdendro::dendro_data(hclust_res)
      )
      tidy_result$dendrogram <- tibble::as_tibble(dendro_data$segments)
      tidy_result$labels <- tibble::as_tibble(dendro_data$labels)
    }, error = function(e) {
      cli::cli_warn("Failed to extract dendrogram data using ggdendro: {e$message}")
    })
  } else {
    # Create basic dendrogram data without ggdendro
    # This is a simplified version that provides basic information
    tidy_result$dendrogram <- tibble::tibble(
      note = "Install 'ggdendro' package for enhanced dendrogram plotting data"
    )

    tidy_result$labels <- tibble::tibble(
      label = hclust_res$labels,
      x = seq_along(hclust_res$labels),
      y = 0
    )
  }

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = hclust_res
    ),
    class = c("glystats_hclust_res", "glystats_res")
  )
}
