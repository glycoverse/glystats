#skip("Too long to run.")
skip_if_not_installed("clusterProfiler")
skip_if_not_installed("org.Hs.eg.db")

# Integration test for GO enrichment
test_that("gly_enrich_go works with protein column (integration)", {
  # Create test data with protein column
  exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)

  result <- suppressMessages(gly_enrich_go(exp_with_protein))
  expect_s3_class(result, "glystats_go_ora_res")
  expect_s3_class(result, "glystats_res")
  expect_true(is.list(result))
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))
})

test_that("gly_enrich_go works with protein column and returns properly formatted results", {
  # Mock the enrichGO function - create a mock data frame
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings({
    exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)
    result <- gly_enrich_go(exp_with_protein)
    expect_s3_class(result, "glystats_go_ora_res")
    expect_s3_class(result, "glystats_res")
    expect_true(is.list(result))
    expect_named(result, c("tidy_result", "raw_result"))
    expect_true("id" %in% colnames(result$tidy_result))  # Check column name cleaning works
  },
  enrichGO = function(...) mock_result,
  .package = "clusterProfiler"
  )
})

test_that("gly_enrich_go handles missing protein information", {
  # Use smaller subset
  exp_small <- test_gp_exp |>
    glyexp::slice_head_var(n = 5)

  # Remove protein column
  exp_no_proteins <- exp_small
  if ("protein" %in% colnames(exp_no_proteins$var_info)) {
    exp_no_proteins$var_info <- exp_no_proteins$var_info |> dplyr::select(-protein)
  }

  # Should error when no protein information is available
  expect_error(
    suppressMessages(gly_enrich_go(exp_no_proteins)),
    "protein.*column not found"
  )
})

test_that("gly_enrich_go filters out NA proteins", {
  # Mock the enrichGO function to capture gene list
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  captured_genes <- NULL

  with_mocked_bindings({
    # Create protein column with some NAs
    exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)
    exp_with_nas <- exp_with_protein |>
      glyexp::mutate_var(protein = ifelse(seq_len(glyexp::n_variables(exp_with_protein)) %% 2 == 0,
                                          NA_character_, protein))

    result <- gly_enrich_go(exp_with_nas)

    expect_s3_class(result, "glystats_go_ora_res")
    expect_s3_class(result, "glystats_res")
    expect_true(is.list(result))
    expect_named(result, c("tidy_result", "raw_result"))
    # Check that NA values were filtered out
    expect_true(all(!is.na(captured_genes)))
  },
  enrichGO = function(gene, ...) {
    captured_genes <<- gene
    mock_result
  },
  .package = "clusterProfiler"
  )
})



# Integration test for KEGG enrichment
test_that("gly_enrich_kegg works with protein column (integration)", {
  # Skip if required packages are not available
  skip_if_not_installed("clusterProfiler")

  # Skip if network connection is not available
  skip_if_offline()

  # Create test data with protein column
  exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)

  # Try KEGG enrichment - might fail due to network issues
  result <- try(suppressMessages(gly_enrich_kegg(exp_with_protein)), silent = TRUE)

  # Skip test if network connection fails
  skip_if(inherits(result, "try-error"), "Network connection to KEGG failed")

  expect_s3_class(result, "glystats_kegg_ora_res")
  expect_s3_class(result, "glystats_res")
  expect_true(is.list(result))
  expect_named(result, c("tidy_result", "raw_result"))
})

test_that("gly_enrich_kegg works with protein column and returns properly formatted results", {
  # Mock the enrichKEGG function - create a mock data frame
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings({
    exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)
    result <- gly_enrich_kegg(exp_with_protein)
    expect_s3_class(result, "glystats_kegg_ora_res")
    expect_s3_class(result, "glystats_res")
    expect_true(is.list(result))
    expect_named(result, c("tidy_result", "raw_result"))
    expect_true("id" %in% colnames(result$tidy_result))  # Check column name cleaning works
  },
  enrichKEGG = function(...) mock_result,
  .package = "clusterProfiler"
  )
})

test_that("gly_enrich_kegg handles missing protein information", {
  # Use smaller subset
  exp_small <- test_gp_exp |>
    glyexp::slice_head_var(n = 5)

  # Remove protein column
  exp_no_proteins <- exp_small
  if ("protein" %in% colnames(exp_no_proteins$var_info)) {
    exp_no_proteins$var_info <- exp_no_proteins$var_info |> dplyr::select(-protein)
  }

  # Should error when no protein information is available
  expect_error(
    suppressMessages(gly_enrich_kegg(exp_no_proteins)),
    "protein.*column not found"
  )
})

test_that("gly_enrich_kegg filters out NA proteins", {
  # Mock the enrichKEGG function to capture gene list
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  captured_genes <- NULL

  with_mocked_bindings({
    # Create protein column with some NAs
    exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)
    exp_with_nas <- exp_with_protein |>
      glyexp::mutate_var(protein = ifelse(seq_len(glyexp::n_variables(exp_with_protein)) %% 2 == 0,
                                          NA_character_, protein))

    result <- gly_enrich_kegg(exp_with_nas)

    expect_s3_class(result, "glystats_kegg_ora_res")
    expect_s3_class(result, "glystats_res")
    expect_true(is.list(result))
    expect_named(result, c("tidy_result", "raw_result"))
    # Check that NA values were filtered out
    expect_true(all(!is.na(captured_genes)))
  },
  enrichKEGG = function(gene, ...) {
    captured_genes <<- gene
    mock_result
  },
  .package = "clusterProfiler"
  )
})

# Integration test for Reactome enrichment
test_that("gly_enrich_reactome works with protein column (integration)", {
  # Skip if required packages are not available
  skip_if_not_installed("ReactomePA")

  # Skip if network connection is not available
  skip_if_offline()

  # Create test data with protein column
  exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)

  # Try Reactome enrichment - might fail due to network issues
  result <- try(suppressMessages(gly_enrich_reactome(exp_with_protein)), silent = TRUE)

  # Skip test if network connection fails
  skip_if(inherits(result, "try-error"), "Network connection to Reactome failed")

  expect_s3_class(result, "glystats_reactome_ora_res")
  expect_s3_class(result, "glystats_res")
  expect_true(is.list(result))
  expect_named(result, c("tidy_result", "raw_result"))
})

test_that("gly_enrich_reactome works with protein column and returns properly formatted results", {
  # Mock the enrichPathway function - create a mock data frame
  mock_result <- new("enrichResult", result = data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  ))

  with_mocked_bindings({
    exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)
    result <- suppressMessages(gly_enrich_reactome(exp_with_protein))
    expect_s3_class(result, "glystats_reactome_ora_res")
    expect_s3_class(result, "glystats_res")
    expect_true(is.list(result))
    expect_named(result, c("tidy_result", "raw_result"))
    expect_true("id" %in% colnames(result$tidy_result))  # Check column name cleaning works
  },
  enrichPathway = function(...) mock_result,
  .package = "ReactomePA"
  )
})

test_that("gly_enrich_reactome handles missing protein information", {
  # Use smaller subset
  exp_small <- test_gp_exp |>
    glyexp::slice_head_var(n = 5)

  # Remove protein column
  exp_no_proteins <- exp_small
  if ("protein" %in% colnames(exp_no_proteins$var_info)) {
    exp_no_proteins$var_info <- exp_no_proteins$var_info |> dplyr::select(-protein)
  }

  # Should error when no protein information is available
  expect_error(
    suppressMessages(gly_enrich_reactome(exp_no_proteins)),
    "protein.*column not found"
  )
})

test_that("gly_enrich_reactome filters out NA proteins", {
  # Mock bitr to capture gene list
  mock_enrich_result <- new("enrichResult", result = data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  ))
  captured_uniprot_ids <- NULL

  with_mocked_bindings({
    with_mocked_bindings({
      # Create protein column with some NAs
      exp_with_protein <- test_gp_exp |> glyexp::slice_head_var(n = 5)
      exp_with_nas <- exp_with_protein |>
        glyexp::mutate_var(protein = ifelse(seq_len(glyexp::n_variables(exp_with_protein)) %% 2 == 0,
                                            NA_character_, protein))

      suppressMessages(result <- gly_enrich_reactome(exp_with_nas))

      expect_s3_class(result, "glystats_reactome_ora_res")
      expect_s3_class(result, "glystats_res")
      expect_true(is.list(result))
      expect_named(result, c("tidy_result", "raw_result"))
      # Check that NA values were filtered out
      expect_true(all(!is.na(captured_uniprot_ids)))
    },
    enrichPathway = function(gene, ...) {
      # check if gene is empty as bitr returns empty df
      if (length(gene) == 0) {
        return(mock_enrich_result)
      }
      # fail if not empty
      stop("enrichPathway should have received an empty gene vector")
    },
    .package = "ReactomePA"
    )
  },
  bitr = function(gene_id, ...) {
    captured_uniprot_ids <<- gene_id
    data.frame(UNIPROT = character(0), ENTREZID = character(0), stringsAsFactors = FALSE)
  },
  .package = "clusterProfiler"
  )
})

test_that("gly_enrich_go_ works correctly", {
  # Create test proteins
  proteins <- c("P12345", "Q67890", "R11111")

  # Mock the enrichGO function
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings({
    result <- suppressMessages(gly_enrich_go_(proteins))
    expect_s3_class(result, "glystats_go_ora_res")
    expect_s3_class(result, "glystats_res")
    expect_true(is.list(result))
    expect_named(result, c("tidy_result", "raw_result"))
    expect_true("id" %in% colnames(result$tidy_result))
  },
  enrichGO = function(...) mock_result,
  .package = "clusterProfiler"
  )
})

test_that("gly_enrich_kegg_ works correctly", {
  # Create test proteins
  proteins <- c("P12345", "Q67890", "R11111")

  # Mock the enrichKEGG function
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings({
    result <- suppressMessages(gly_enrich_kegg_(proteins))
    expect_s3_class(result, "glystats_kegg_ora_res")
    expect_s3_class(result, "glystats_res")
    expect_true(is.list(result))
    expect_named(result, c("tidy_result", "raw_result"))
    expect_true("id" %in% colnames(result$tidy_result))
  },
  enrichKEGG = function(...) mock_result,
  .package = "clusterProfiler"
  )
})

test_that("gly_enrich_reactome_ works correctly", {
  # Create test proteins
  proteins <- c("P12345", "Q67890", "R11111")

  # Mock the functions
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings({
    with_mocked_bindings({
      result <- suppressMessages(gly_enrich_reactome_(proteins))
      expect_s3_class(result, "glystats_reactome_ora_res")
      expect_s3_class(result, "glystats_res")
      expect_true(is.list(result))
      expect_named(result, c("tidy_result", "raw_result"))
      expect_true("id" %in% colnames(result$tidy_result))
    },
    enrichPathway = function(...) mock_result,
    .package = "ReactomePA"
    )
  },
  bitr = function(...) data.frame(UNIPROT = character(0), ENTREZID = character(0), stringsAsFactors = FALSE),
  .package = "clusterProfiler"
  )
})

test_that("gly_enrich_go_ supports custom OrgDb", {
  skip_if_not_installed("org.Mm.eg.db")
  proteins <- c("Q9WV54", "Q02788", "Q3V3R4", "Q61207", "P97821", "P02468", "P09470", "Q9CYA0")
  result <- suppressMessages(gly_enrich_go_(proteins, OrgDb = "org.Mm.eg.db"))
  expect_s3_class(result, "glystats_go_ora_res")
  expect_s3_class(result, "glystats_res")
  expect_true(is.list(result))
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true("id" %in% colnames(result$tidy_result))
})