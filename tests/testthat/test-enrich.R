# Integration test for GO enrichment
test_that("enrich_go works with protein column (integration)", {
  # Skip if required packages are not available
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Use glyclean::infer_protein to convert genes column to protein column
  exp_with_protein <- glyclean::infer_protein(test_gp_exp |> glyexp::slice_head_var(n = 5))
  
  result <- suppressMessages(enrich_go(exp_with_protein))
  expect_s3_class(result, "glystats_go_ora_res")
})

test_that("enrich_go works with genes column and infer_protein", {
  # Mock the enrichGO function - create a mock data frame
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  
  with_mocked_bindings({
    result <- enrich_go(test_gp_exp |> glyexp::slice_head_var(n = 5))
    expect_s3_class(result, "glystats_go_ora_res")
    expect_true("id" %in% colnames(result))  # Check column name cleaning works
  },
  enrichGO = function(...) mock_result,
  .package = "clusterProfiler"
  )
})

test_that("enrich_go handles missing protein information", {
  # Use smaller subset
  exp_small <- test_gp_exp |> 
    glyexp::slice_head_var(n = 5)
  
  # Remove both protein and proteins columns
  exp_no_proteins <- exp_small
  if ("protein" %in% colnames(exp_no_proteins$var_info)) {
    exp_no_proteins$var_info <- exp_no_proteins$var_info |> dplyr::select(-protein)
  }
  if ("proteins" %in% colnames(exp_no_proteins$var_info)) {
    exp_no_proteins$var_info <- exp_no_proteins$var_info |> dplyr::select(-proteins)
  }
  
  # Should error when no protein information is available
  expect_error(
    enrich_go(exp_no_proteins),
    "protein.*or.*proteins.*column not found"
  )
})

test_that("enrich_go passes additional arguments to enrichGO", {
  # Mock the enrichGO function to capture arguments
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  enrichGO_args <- NULL
  
  with_mocked_bindings({
    exp_with_protein <- glyclean::infer_protein(test_gp_exp |> glyexp::slice_head_var(n = 3))
    result <- enrich_go(exp_with_protein, pvalueCutoff = 0.01, ont = "BP")
    
    expect_s3_class(result, "glystats_go_ora_res")
    expect_equal(enrichGO_args$pvalueCutoff, 0.01)
    expect_equal(enrichGO_args$ont, "BP")
    expect_equal(enrichGO_args$keyType, "UNIPROT")
  },
  enrichGO = function(...) {
    enrichGO_args <<- list(...)
    mock_result
  },
  .package = "clusterProfiler"
  )
})

test_that("enrich_go filters out NA proteins", {
  # Mock the enrichGO function to capture gene list
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  captured_genes <- NULL
  
  with_mocked_bindings({
    # Convert to protein column first, then add some NAs
    exp_with_protein <- glyclean::infer_protein(test_gp_exp |> glyexp::slice_head_var(n = 5))
    exp_with_nas <- exp_with_protein |>
      glyexp::mutate_var(protein = ifelse(seq_len(glyexp::n_variables(exp_with_protein)) %% 2 == 0, 
                                         NA_character_, protein))
    
    result <- enrich_go(exp_with_nas)
    
    expect_s3_class(result, "glystats_go_ora_res")
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
test_that("enrich_kegg works with protein column (integration)", {
  # Skip if required packages are not available
  skip_if_not_installed("clusterProfiler")
  
  # Skip if network connection is not available
  skip_if_offline()
  
  # Use glyclean::infer_protein to convert genes column to protein column
  exp_with_protein <- glyclean::infer_protein(test_gp_exp |> glyexp::slice_head_var(n = 5))
  
  # Try KEGG enrichment - might fail due to network issues
  result <- try(suppressMessages(enrich_kegg(exp_with_protein)), silent = TRUE)
  
  # Skip test if network connection fails
  skip_if(inherits(result, "try-error"), "Network connection to KEGG failed")
  
  expect_s3_class(result, "glystats_kegg_ora_res")
})

test_that("enrich_kegg works with genes column and infer_protein", {
  # Mock the enrichKEGG function - create a mock data frame
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  
  with_mocked_bindings({
    result <- enrich_kegg(test_gp_exp |> glyexp::slice_head_var(n = 5))
    expect_s3_class(result, "glystats_kegg_ora_res")
    expect_true("id" %in% colnames(result))  # Check column name cleaning works
  },
  enrichKEGG = function(...) mock_result,
  .package = "clusterProfiler"
  )
})

test_that("enrich_kegg handles missing protein information", {
  # Use smaller subset
  exp_small <- test_gp_exp |> 
    glyexp::slice_head_var(n = 5)
  
  # Remove both protein and proteins columns
  exp_no_proteins <- exp_small
  if ("protein" %in% colnames(exp_no_proteins$var_info)) {
    exp_no_proteins$var_info <- exp_no_proteins$var_info |> dplyr::select(-protein)
  }
  if ("proteins" %in% colnames(exp_no_proteins$var_info)) {
    exp_no_proteins$var_info <- exp_no_proteins$var_info |> dplyr::select(-proteins)
  }
  
  # Should error when no protein information is available
  expect_error(
    enrich_kegg(exp_no_proteins),
    "protein.*or.*proteins.*column not found"
  )
})

test_that("enrich_kegg passes additional arguments to enrichKEGG", {
  # Mock the enrichKEGG function to capture arguments
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  enrichKEGG_args <- NULL
  
  with_mocked_bindings({
    exp_with_protein <- glyclean::infer_protein(test_gp_exp |> glyexp::slice_head_var(n = 3))
    result <- enrich_kegg(exp_with_protein, pvalueCutoff = 0.001, qvalueCutoff = 0.1)
    
    expect_s3_class(result, "glystats_kegg_ora_res")
    expect_equal(enrichKEGG_args$pvalueCutoff, 0.001)
    expect_equal(enrichKEGG_args$qvalueCutoff, 0.1)
  },
  enrichKEGG = function(...) {
    enrichKEGG_args <<- list(...)
    mock_result
  },
  .package = "clusterProfiler"
  )
})

test_that("enrich_kegg filters out NA proteins", {
  # Mock the enrichKEGG function to capture gene list
  mock_result <- data.frame(
    ID = character(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  captured_genes <- NULL
  
  with_mocked_bindings({
    # Convert to protein column first, then add some NAs
    exp_with_protein <- glyclean::infer_protein(test_gp_exp |> glyexp::slice_head_var(n = 5))
    exp_with_nas <- exp_with_protein |>
      glyexp::mutate_var(protein = ifelse(seq_len(glyexp::n_variables(exp_with_protein)) %% 2 == 0, 
                                         NA_character_, protein))
    
    result <- enrich_kegg(exp_with_nas)
    
    expect_s3_class(result, "glystats_kegg_ora_res")
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
