#' K-means Clustering for Glycomics and Glycoproteomics Data
#'
#' Perform k-means clustering on the expression data.
#' The function uses `stats::kmeans()` to perform clustering and provides
#' tidy results with cluster assignments.
#'
#' @param exp A `glyexp::experiment()` object containing expression matrix and sample information.
#' @param on A character string specifying what to cluster. Either "variable" (default) to cluster
#'   variables/features, or "sample" to cluster samples/observations.
#' @param centers Either the number of clusters (integer) or a set of initial cluster centers.
#'   Default is 3.
#' @param iter.max The maximum number of iterations allowed. Default is 10.
#' @param nstart How many random sets should be chosen if centers is a number. Default is 1.
#' @param algorithm Character string specifying the algorithm to use. May be abbreviated.
#'   Default is "Hartigan-Wong".
#' @param scale A logical indicating whether to scale the data before clustering. Default is TRUE.
#' @param add_info A logical value. If TRUE (default), sample information from the experiment
#'   will be added to the result tibbles. If FALSE, only the clustering results are returned.
#' @param return_raw A logical value. If FALSE (default), returns processed tibble results.
#'   If TRUE, returns raw kmeans object.
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
#' **Data Preparation:**
#' Data is log2-transformed and optionally scaled before clustering.
#'
#' **Clustering Method:**
#' K-means clustering is performed using `stats::kmeans()` with the specified parameters.
#'
#' @return A tibble with two columns (when return_raw = FALSE):
#'  - First column: variable or sample names (depending on `on` parameter)
#'  - Second column: cluster assignments
#' When return_raw = TRUE, returns the raw kmeans object.
#'
#' @seealso [stats::kmeans()]
#' @export
gly_kmeans <- function(
  exp,
  on = "variable",
  centers = 3,
  iter.max = 10,
  nstart = 1,
  algorithm = "Hartigan-Wong",
  scale = TRUE,
  add_info = TRUE,
  return_raw = FALSE,
  ...
) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_choice(on, c("variable", "sample"))
  checkmate::assert(
    checkmate::check_integerish(centers, lower = 1, len = 1),
    checkmate::check_matrix(centers)
  )
  checkmate::assert_integerish(iter.max, lower = 1, len = 1)
  checkmate::assert_integerish(nstart, lower = 1, len = 1)
  checkmate::assert_string(algorithm)
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

  # Perform k-means clustering
  kmeans_res <- stats::kmeans(
    mat, 
    centers = centers,
    iter.max = iter.max,
    nstart = nstart,
    algorithm = algorithm,
    ...
  )

  # Return raw results if requested
  if (return_raw) {
    return(kmeans_res)
  }

  # Create result tibble
  result_tbl <- tibble::tibble(
    x = rownames(mat),
    cluster = kmeans_res$cluster
  )

  # Set the appropriate column name based on clustering type
  colnames(result_tbl)[1] <- cluster_type

  # Add info if requested
  if (add_info) {
    if (cluster_type == "sample") {
      sample_info <- glyexp::get_sample_info(exp)
      if (!is.null(sample_info) && nrow(sample_info) > 0) {
        result_tbl <- result_tbl %>%
          dplyr::left_join(sample_info, by = "sample")
      }
    } else {
      var_info <- glyexp::get_var_info(exp)
      if (!is.null(var_info) && nrow(var_info) > 0) {
        result_tbl <- result_tbl %>%
          dplyr::left_join(var_info, by = "variable")
      }
    }
  }

  return(result_tbl)
}
