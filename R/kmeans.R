#' K-means Clustering for Glycomics and Glycoproteomics Data
#'
#' Perform k-means clustering on the expression data.
#' The function uses `stats::kmeans()` to perform clustering and provides
#' tidy results with cluster assignments.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param expr_mat A numeric matrix with variables as rows and samples as columns.
#' @param on A character string specifying what to cluster. Either "variable" (default) to cluster
#'   variables/features, or "sample" to cluster samples/observations.
#' @param centers Either the number of clusters (integer) or a set of initial cluster centers.
#'   Default is 3.
#' @param scale A logical indicating whether to scale the data before clustering. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'   will be added to the result tibbles. If FALSE, only the clustering results are returned.
#'   Only applicable to `gly_kmeans()`.
#' @param ... Additional arguments passed to `stats::kmeans()`.
#'
#' @section Required packages:
#' This function only uses base R packages and does not require additional dependencies.
#'
#' @details
#' The function performs log2 transformation on the expression data (log2(x + 1)) before
#' clustering. When `on = "variable"` (default), variables are clustered based on their
#' expression patterns across samples. When `on = "sample"`, samples are clustered based
#' on their expression profiles across variables.
#'
#' `gly_kmeans()` is the top-level API that works with `glyexp::experiment()` objects and supports
#' the `add_info` parameter for joining experiment metadata.
#'
#' `gly_kmeans_()` is the underlying API that works with matrices directly,
#' providing more flexibility for users who don't use the glyexp package.
#'
#' **Data Preparation:**
#' Data is log2-transformed and optionally scaled before clustering.
#'
#' **Clustering Method:**
#' K-means clustering is performed using `stats::kmeans()` with the specified parameters.
#'
#' @returns A list with two elements:
#'  - `tidy_result`: A tibble with cluster assignments containing the following columns:
#'    - `variable` or `sample`: Variable or sample name (depending on `on` parameter)
#'    - `cluster`: Cluster assignment
#'  - `raw_result`: The raw kmeans object from `stats::kmeans()`.
#'
#' @seealso [stats::kmeans()]
#' @export
gly_kmeans <- function(
  exp,
  on = "variable",
  centers = 3,
  scale = TRUE,
  add_info = TRUE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_choice(on, c("variable", "sample"))
  checkmate::assert(
    checkmate::check_integerish(centers, lower = 1, len = 1),
    checkmate::check_matrix(centers)
  )
  checkmate::assert_logical(scale, len = 1)
  checkmate::assert_logical(add_info, len = 1)

  # Extract data from experiment object
  expr_mat <- glyexp::get_expr_mat(exp)

  # Call the underlying API
  result <- gly_kmeans_(expr_mat, on, centers, scale, ...)
  result$tidy_result <- .process_results_add_info(result$tidy_result, exp, add_info)
  result
}

#' @rdname gly_kmeans
#' @export
gly_kmeans_ <- function(
  expr_mat,
  on = "variable",
  centers = 3,
  scale = TRUE,
  ...
) {
  # Validate inputs
  checkmate::assert_matrix(expr_mat, mode = "numeric")
  checkmate::assert_choice(on, c("variable", "sample"))
  checkmate::assert(
    checkmate::check_integerish(centers, lower = 1, len = 1),
    checkmate::check_matrix(centers)
  )
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

  # Perform k-means clustering
  kmeans_res <- stats::kmeans(
    mat,
    centers = centers,
    ...
  )

  # Create result tibble
  result_tbl <- tibble::tibble(
    x = rownames(mat),
    cluster = kmeans_res$cluster
  )

  # Set the appropriate column name based on clustering type
  colnames(result_tbl)[1] <- cluster_type

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = result_tbl,
      raw_result = kmeans_res
    ),
    class = c("glystats_kmeans_res", "glystats_res")
  )
}
