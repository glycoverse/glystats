#' Map abundance profiles with a self-organizing map
#'
#' Train an unsupervised self-organizing map (SOM) and assign each profile
#' to its closest grid unit. Training favors similar prototypes at neighboring
#' units.
#'
#' @inheritParams gly_kmeans
#' @inheritParams gly_scale
#' @param on Axis to map. `"variable"` maps each glycan or glycoform across
#'   the supplied samples. `"sample"` maps each sample across variables.
#' @param grid Grid dimensions as two positive integers, `c(xdim, ydim)`.
#'   By default, the map is square with side length
#'   `min(5, floor(sqrt(n_train)))`, where `n_train` counts usable profiles.
#' @param topology Arrangement of the grid units, either `"rectangular"`
#'   (default) or `"hexagonal"`.
#' @param rlen Number of training epochs. Each epoch makes as many online
#'   updates as there are training profiles.
#' @param scale Whether to center and standardize each input coordinate
#'   across the retained profiles, as in [gly_kmeans()]. With at least two
#'   training profiles, each sample (variable map) or variable (sample map)
#'   is centered at its mean and divided by its sample standard deviation.
#'   Constant coordinates are centered at zero. A single profile is used
#'   directly.
#' @param seed Nonnegative integer seed for initialization, online training
#'   and tie breaking. The default restores the caller's random-number state
#'   after fitting, as does any supplied integer. Use `NULL` to draw from the
#'   caller's random-number stream.
#' @param add_info Whether to attach annotations for the mapped axis.
#'
#' @details
#' The columns of the input assay represent individual samples. To map
#' group-average curves, first supply an experiment with one column per
#' group mean. Each supplied coordinate contributes to Euclidean distance.
#' Standardization gives each varying coordinate equal variance. With
#' `scale = FALSE`, coordinates with larger spreads can dominate the map.
#' With relative abundances, a change in one glycan's share can affect the
#' shares of other glycans. Distances describe similarity relative to the
#' normalization reference and depend on the abundance transformation.
#'
#' Raw assays use `log2(x + 1e-6)` by default. See [gly_scale] for other
#' scales and overrides. On the analysis scale, coordinates with no finite
#' values are removed first. Training uses profiles with finite values at
#' every remaining coordinate. A complete flat profile is eligible because
#' its level contributes to distance. Initial grid prototypes are sampled
#' without replacement from the training profiles. Changing the grid or seed
#' can change the assignments, and unit numbers are local to each fit.
#'
#' @returns A `glystats_som_res` object inheriting from `glystats_res`.
#'   [get_tidy_result()] returns a tibble in input order with one row per
#'   mapped-axis ID (`.variable` or `.sample` for SE inputs).
#'
#'   * `unit`, `x`, and `y` identify the winning unit and its grid coordinates.
#'   * `distance` is the Euclidean distance to that unit's fitted prototype
#'     in the training coordinates. Smaller values indicate a closer fit.
#'   * `n_finite` counts finite coordinates. Excluded profiles retain their
#'     IDs and receive missing assignments and distances. `exclusion` is
#'     `"no_finite_values"` or `"incomplete_profile"` for these profiles
#'     and `NA` for the training profiles.
#'
#'   `diagnostics` contains `n_train`, `units` (occupancy of every grid unit,
#'   including empty units), and `quantization_error` (mean training distance).
#'   [get_raw_result()] returns the [kohonen::som()] fit with its training data,
#'   fitted assignments, codebook vectors, and per-epoch changes. With zero
#'   usable profiles, the automatic grid is empty, quantization error is
#'   missing, and the fit is `NULL`.
#'   An explicit grid requires at least one usable profile per requested unit.
#'
#'   `settings` records the axis, grid dimensions, topology, epochs,
#'   standardization request and seed. `preprocessing` records the abundance
#'   transformation, whether standardization was applied (`scaled`), and the
#'   IDs of wholly non-finite coordinates (`excluded_coordinates`). The
#'   result retains the input `meta_data` and SE [identifiers][gly_identity].
#'
#' @references
#' Kohonen T (1982). Self-organized formation of topologically correct
#' feature maps. *Biological Cybernetics*, 43, 59-69.
#' \doi{10.1007/BF00337288}.
#'
#' [kohonen on CRAN](https://CRAN.R-project.org/package=kohonen).
#'
#' @examplesIf requireNamespace("kohonen", quietly = TRUE)
#' exp <- glyexp::real_experiment2
#' res <- gly_som(exp, grid = c(2, 2), seed = 42, add_info = FALSE)
#' head(get_tidy_result(res))
#' res$diagnostics$units
#' @export
gly_som <- function(
  exp,
  on = "variable",
  grid = NULL,
  topology = "rectangular",
  rlen = 100L,
  scale = TRUE,
  seed = 1L,
  add_info = TRUE,
  log_transform = NULL
) {
  .assert_data_container(exp)
  checkmate::assert_choice(on, c("variable", "sample"))
  checkmate::assert_integerish(
    grid,
    len = 2L,
    lower = 1L,
    upper = .Machine$integer.max,
    any.missing = FALSE,
    null.ok = TRUE
  )
  checkmate::assert_choice(topology, c("rectangular", "hexagonal"))
  checkmate::assert_int(rlen, lower = 1L)
  checkmate::assert_flag(scale)
  checkmate::assert_int(seed, lower = 0L, null.ok = TRUE)
  checkmate::assert_flag(add_info)
  rlang::check_installed("kohonen", version = "3.0.0")

  meta_data <- .get_meta_data(exp)
  preprocessing <- .resolve_transformation(
    exp,
    log_transform,
    meta_data = meta_data
  )
  mat <- .get_expr_mat(exp)
  checkmate::assert_matrix(mat, mode = "numeric")
  mat <- .log_transform_expr_mat(mat, preprocessing$log_transform)
  if (on == "sample") {
    mat <- t(mat)
  }

  finite <- is.finite(mat)
  coordinates <- colSums(finite) > 0L
  n_finite <- as.integer(rowSums(finite))
  usable <- n_finite > 0L & n_finite == sum(coordinates)
  training <- mat[usable, coordinates, drop = FALSE]
  n_train <- nrow(training)
  preprocessing$excluded_coordinates <- as.character(colnames(mat)[
    !coordinates
  ])

  if (is.null(grid)) {
    grid <- rep(min(5L, floor(sqrt(n_train))), 2L)
  }
  grid <- as.integer(grid)
  n_units <- prod(grid)
  if (n_units > n_train) {
    cli::cli_abort(
      "The grid requests {n_units} unit{?s}. The data contain {n_train} usable training profile{?s}."
    )
  }

  preprocessing$scaled <- scale && n_train > 1L
  if (preprocessing$scaled) {
    training <- .scale_cluster_data(training)
  }
  fit <- NULL
  positions <- matrix(numeric(), ncol = 2L)
  unit_id <- rep(NA_integer_, nrow(mat))
  distance <- rep(NA_real_, nrow(mat))
  if (n_train > 0L) {
    som_grid <- kohonen::somgrid(grid[1L], grid[2L], topo = topology)
    if (!is.null(seed)) {
      withr::local_seed(seed)
    }
    fit <- kohonen::som(
      training,
      grid = som_grid,
      rlen = rlen,
      dist.fcts = "euclidean",
      keep.data = TRUE
    )
    positions <- som_grid$pts
    unit_id[usable] <- as.integer(fit$unit.classif)
    distance[usable] <- fit$distances
  }

  exclusion <- rep(NA_character_, nrow(mat))
  exclusion[!usable] <- "incomplete_profile"
  exclusion[n_finite == 0L] <- "no_finite_values"
  tidy <- tibble::tibble(
    id = as.character(rownames(mat)),
    unit = unit_id,
    x = positions[unit_id, 1L],
    y = positions[unit_id, 2L],
    distance = distance,
    n_finite = n_finite,
    exclusion = exclusion
  )
  names(tidy)[1L] <- on
  units <- tibble::tibble(
    unit = seq_len(n_units),
    x = positions[, 1L],
    y = positions[, 2L],
    n_profiles = tabulate(unit_id, nbins = n_units)
  )
  result <- list(
    tidy_result = .process_results_add_info(tidy, exp, add_info),
    raw_result = fit,
    diagnostics = list(
      n_train = n_train,
      units = units,
      quantization_error = if (n_train > 0L) {
        mean(distance[usable])
      } else {
        NA_real_
      }
    ),
    settings = list(
      on = on,
      grid = grid,
      topology = topology,
      rlen = rlen,
      scale = scale,
      seed = seed
    ),
    meta_data = meta_data,
    preprocessing = preprocessing
  )
  class(result) <- c("glystats_som_res", "glystats_res")
  .record_se_identity(result, exp)
}
