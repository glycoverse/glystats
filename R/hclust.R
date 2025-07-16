#' Hierarchical Clustering for Glycomics and Glycoproteomics Data
#'
#' Perform hierarchical clustering on the expression data.
#' The function uses `stats::hclust()` to perform clustering and provides
#' tidy results including cluster assignments, dendrogram data for plotting,
#' and merge heights.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param on A character string specifying what to cluster. Either "variable" (default) to cluster
#'   variables/features, or "sample" to cluster samples/observations.
#' @param method The agglomeration method to be used. This should be one of
#'   "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid".
#'   Default is "complete".
#' @param dist_method The distance measure to be used. This must be one of
#'   "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski".
#'   Default is "euclidean".
#' @param k_values A numeric vector specifying the number of clusters to cut the tree into.
#'   Default is c(2, 3, 4, 5). If NULL, no cluster assignments are returned.
#' @param scale A logical indicating whether to scale the data before clustering. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'   will be added to the result tibbles. If FALSE, only the clustering results are returned.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw hclust object.
#' @param ... Additional arguments passed to `stats::dist()` or `stats::hclust()`.
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
#' @return A list containing multiple tibbles (when return_raw = FALSE):
#'  - `clusters`: Cluster assignments for different k values (variables or samples depending on `on` parameter)
#'  - `dendrogram`: Dendrogram segment data for plotting (if ggdendro is available)
#'  - `heights`: Merge heights and steps for the clustering process
#'  - `labels`: Labels and their positions (if ggdendro is available)
#' When return_raw = TRUE, returns the raw hclust object.
#'
#' @seealso [stats::hclust()], [stats::dist()], [stats::cutree()]
#' @export
gly_hclust <- function(
  exp,
  on = "variable",
  method = "complete",
  dist_method = "euclidean",
  k_values = c(2, 3, 4, 5),
  scale = TRUE,
  add_info = TRUE,
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_choice(on, c("variable", "sample"))
  checkmate::assert_string(method)
  checkmate::assert_string(dist_method)
  if (!is.null(k_values)) {
    checkmate::assert_integerish(k_values, lower = 2, min.len = 1)
  }
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(add_info, len = 1)
  checkmate::assert_logical(return_raw, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)

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

  # Calculate distance matrix
  dist_mat <- stats::dist(mat, method = dist_method, ...)

  # Perform hierarchical clustering
  hclust_res <- stats::hclust(dist_mat, method = method, ...)

  # Return raw results if requested
  if (return_raw) {
    return(hclust_res)
  }

  # Initialize result list
  result <- list()

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
    result$clusters <- clusters_tbl
  }

  # 2. Heights data
  result$heights <- tibble::tibble(
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
      result$dendrogram <- tibble::as_tibble(dendro_data$segments)
      result$labels <- tibble::as_tibble(dendro_data$labels)
    }, error = function(e) {
      cli::cli_warn("Failed to extract dendrogram data using ggdendro: {e$message}")
    })
  } else {
    # Create basic dendrogram data without ggdendro
    # This is a simplified version that provides basic information
    result$dendrogram <- tibble::tibble(
      note = "Install 'ggdendro' package for enhanced dendrogram plotting data"
    )

    result$labels <- tibble::tibble(
      label = hclust_res$labels,
      x = seq_along(hclust_res$labels),
      y = 0
    )
  }

  result <- .process_results_add_info(result, exp, add_info)
  structure(result, class = c("glystats_hclust_res", "glystats_res"))
}
