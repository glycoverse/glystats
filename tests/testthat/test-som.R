som_test_se <- function() {
  values <- rbind(
    g4 = c(1, 4, 7, 3),
    g1 = c(8, 5, 2, 11),
    g6 = c(9, 3, 4, 8),
    g2 = c(2, 6, 1, 5),
    g5 = c(7, 8, 5, 2),
    g3 = c(4, 3, 9, 6)
  )
  colnames(values) <- c("s3", "s1", "s4", "s2")
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = values),
    rowData = S4Vectors::DataFrame(variable = letters[1:6], unit = 11:16),
    colData = S4Vectors::DataFrame(
      sample = LETTERS[1:4],
      unit = 21:24,
      group = c("A", "A", "B", "B")
    )
  )
}

test_that("SOM maps the six variables by default and the four samples on request", {
  skip_if_not_installed("kohonen", "3.0.0")
  se <- som_test_se()
  values <- log2(SummarizedExperiment::assay(se) + 1e-6)
  for (on in c("variable", "sample")) {
    if (on == "variable") {
      res <- gly_som(
        se,
        grid = c(2, 1),
        topology = "hexagonal",
        rlen = 12,
        seed = 9
      )
      mat <- scale(values)
      annotation <- SummarizedExperiment::rowData(se)
    } else {
      res <- gly_som(
        se,
        on = "sample",
        grid = c(2, 1),
        topology = "hexagonal",
        rlen = 12,
        seed = 9
      )
      mat <- scale(t(values))
      annotation <- SummarizedExperiment::colData(se)
    }
    expected <- withr::with_seed(
      9,
      kohonen::som(
        mat,
        grid = kohonen::somgrid(2, 1, "hexagonal"),
        rlen = 12,
        dist.fcts = "euclidean"
      )
    )
    tidy <- get_tidy_result(res)
    fit <- get_raw_result(res)
    expect_s3_class(res, "glystats_som_res")
    expect_identical(tidy[[paste0(".", on)]], rownames(mat))
    expect_equal(tidy[[on]], annotation[[on]])
    expect_identical(tidy$unit_metadata, annotation$unit)
    expect_identical(tidy$unit, as.integer(expected$unit.classif))
    expect_equal(fit$data[[1L]], mat)
    expect_equal(kohonen::getCodes(fit), kohonen::getCodes(expected))
    coords <- expected$grid$pts[tidy$unit, , drop = FALSE]
    expect_equal(as.matrix(tidy[c("x", "y")]), coords, ignore_attr = TRUE)
    prototype <- kohonen::getCodes(expected)[tidy$unit, , drop = FALSE]
    expect_equal(tidy$distance, unname(sqrt(rowSums((mat - prototype)^2))))
    expect_equal(res$diagnostics$quantization_error, mean(tidy$distance))
    expect_identical(res$diagnostics$units$n_profiles, tabulate(tidy$unit, 2L))
    expect_identical(res$preprocessing$analysis_scale, "log2")
    expect_false(any(c("p_val", "p_adj") %in% names(tidy)))
  }
})

test_that("SOM filters before initialization and traces exclusions on either axis", {
  skip_if_not_installed("kohonen", "3.0.0")
  for (on in c("variable", "sample")) {
    se <- som_test_se()
    rownames(se) <- rep("duplicate", nrow(se))
    colnames(se) <- NULL
    mat <- SummarizedExperiment::assay(se)
    if (on == "sample") {
      mat <- t(mat)
    }
    mat[1, ] <- 9 # Complete flat profile
    mat[, 1] <- 9 # Constant coordinate among training profiles
    mat[, 2] <- NaN # Wholly unavailable coordinate
    mat[2, ] <- NA_real_ # Empty profile
    mat[4, 3] <- Inf # Incomplete profile
    SummarizedExperiment::assay(se) <- if (on == "variable") mat else t(mat)
    kept <- setdiff(seq_len(nrow(mat)), c(2L, 4L))
    res <- gly_som(se, on = on, log_transform = FALSE, add_info = FALSE)
    tidy <- get_tidy_result(res)
    fit <- get_raw_result(res)
    ids <- res$identifiers[[on]]$data
    coord_axis <- if (on == "variable") "sample" else "variable"
    coord_ids <- res$identifiers[[coord_axis]]$data$id
    expect_identical(ids$position, seq_len(nrow(mat)))
    expect_identical(ids$id, paste0(".glystats_", on, "_", ids$position))
    expect_identical(tidy[[paste0(".", on)]], ids$id)
    expect_identical(rownames(fit$data[[1L]]), ids$id[kept])
    expect_named(
      tidy,
      c(paste0(".", on), "unit", "x", "y", "distance", "n_finite", "exclusion")
    )
    expect_identical(res$preprocessing$excluded_coordinates, coord_ids[2L])
    expect_identical(colnames(fit$data[[1L]]), coord_ids[-2L])
    expect_equal(unname(fit$data[[1L]][, 1L]), rep(0, length(kept)))
    expect_identical(tidy$n_finite, as.integer(rowSums(is.finite(mat))))
    expect_identical(
      tidy$exclusion[c(2, 4)],
      c("no_finite_values", "incomplete_profile")
    )
    expect_true(all(is.na(tidy$exclusion[kept])))
    expect_true(all(is.na(tidy[c(2, 4), c("unit", "x", "y", "distance")])))
    expect_identical(tidy$unit[kept], as.integer(fit$unit.classif))
    expect_identical(res$diagnostics$n_train, length(kept))
    expect_equal(res$settings$grid, rep(floor(sqrt(length(kept))), 2L))
    expect_error(
      gly_som(
        se,
        on = on,
        grid = c(length(kept) + 1L, 1L),
        log_transform = FALSE
      ),
      paste(length(kept), "usable training profiles")
    )
  }
})

test_that("automatic SOM grids handle one and zero usable profiles", {
  skip_if_not_installed("kohonen", "3.0.0")
  se <- som_test_se()
  SummarizedExperiment::assay(se)[-3, ] <- NA_real_
  res <- gly_som(se, add_info = FALSE)
  tidy <- get_tidy_result(res)
  expect_identical(res$settings$grid, c(1L, 1L))
  expect_false(res$preprocessing$scaled)
  expect_identical(which(!is.na(tidy$unit)), 3L)
  expect_identical(tidy$unit[3], 1L)
  expect_equal(tidy$distance[3], 0)
  expect_error(gly_som(se, grid = c(2, 1)), "1 usable training profile")

  SummarizedExperiment::assay(se)[3, ] <- NA_real_
  empty <- gly_som(se)
  expect_null(get_raw_result(empty))
  expect_identical(get_tidy_result(empty)$.variable, rownames(se))
  expect_true(all(get_tidy_result(empty)$exclusion == "no_finite_values"))
  expect_identical(empty$diagnostics$n_train, 0L)
  expect_equal(nrow(empty$diagnostics$units), 0)
  expect_identical(empty$diagnostics$quantization_error, NA_real_)
  expect_identical(empty$preprocessing$excluded_coordinates, colnames(se))
  expect_error(gly_som(se, grid = c(1, 1)), "0 usable training profiles")
})

test_that("SOM honors transformed coordinates, standardization and the seed", {
  skip_if_not_installed("kohonen", "3.0.0")
  se <- som_test_se()
  S4Vectors::metadata(se)$scale <- "logx"
  SummarizedExperiment::assay(se) <- SummarizedExperiment::assay(se) - 5
  withr::local_seed(381)
  rng <- .Random.seed
  res <- gly_som(se, scale = FALSE)
  expect_identical(.Random.seed, rng)
  expect_identical(res$preprocessing$input_scale, "logx")
  expect_false(res$preprocessing$log_transform)
  expect_false(res$preprocessing$scaled)
  expect_equal(get_raw_result(res)$data[[1L]], SummarizedExperiment::assay(se))
  stats::runif(1)
  expect_warning(
    again <- gly_som(se, scale = FALSE, log_transform = TRUE),
    class = "glystats_already_logged"
  )
  expect_identical(get_tidy_result(res), get_tidy_result(again))
  expect_identical(get_raw_result(res)$codes, get_raw_result(again)$codes)

  withr::local_seed(72)
  rng <- .Random.seed
  from_caller <- gly_som(se, seed = NULL)
  expect_false(identical(.Random.seed, rng))
  replay <- withr::with_seed(72, gly_som(se, seed = NULL))
  expect_identical(get_tidy_result(from_caller), get_tidy_result(replay))
})

test_that("SOM maps real glycomics and glycoproteomics experiments", {
  skip_if_not_installed("kohonen", "3.0.0")
  se <- glyexp::real_experiment2
  res <- gly_som(se, rlen = 10)
  tidy <- get_tidy_result(res)
  values <- SummarizedExperiment::assay(se)
  complete <- rowSums(!is.finite(values)) == 0L
  expect_true(any(complete) && any(!complete))
  expect_identical(tidy$.variable, rownames(se))
  expect_identical(!is.na(tidy$unit), unname(complete))
  expect_equal(
    tidy$glycan_composition,
    SummarizedExperiment::rowData(se)$glycan_composition
  )
  expect_identical(res$meta_data, S4Vectors::metadata(se))

  samples <- gly_som(test_gp_se, on = "sample", rlen = 10)
  expect_identical(get_tidy_result(samples)$.sample, colnames(test_gp_se))
  expect_identical(
    get_tidy_result(samples)$group,
    SummarizedExperiment::colData(test_gp_se)$group
  )
  expect_true(all(is.finite(get_tidy_result(samples)$distance)))
})

test_that("SOM validates grid dimensions and protects fitted assignments", {
  se <- som_test_se()
  for (grid in list(2L, c(1, 0), c(1, NA), c(1, 2.5))) {
    expect_error(gly_som(se, grid = grid), "grid")
  }
  expect_error(gly_som(se, on = "group"), "on")
  expect_error(gly_som(se, topology = "circle"), "topology")
  expect_error(gly_som(se, rlen = 0), "rlen")
  expect_error(gly_som(se, seed = NA_integer_), "seed")
  expect_error(gly_som(se, keep.data = FALSE), "unused argument.*keep.data")
  expect_error(
    gly_som(se, maxNA.fraction = 1),
    "unused argument.*maxNA.fraction"
  )
})

test_that("SOM identifies the missing optional backend", {
  skip_if(
    rlang::is_installed("kohonen", version = "3.0.0"),
    "requires a library without kohonen"
  )
  expect_error(
    gly_som(som_test_se()),
    "kohonen",
    class = "rlib_error_package_not_found"
  )
})
