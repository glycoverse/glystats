#' GO over-representation analysis (ORA)
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function was deprecated because we decided to move all enrichment analysis functions
#' to the separate `glyfun` package, which has more features and better API design.
#'
#' Perform GO ORA for protein UniProt accessions using [clusterProfiler::enrichGO()].
#' - `gly_enrich_go()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::enrichGO()].
#' @param ont Passed to `ont` of [clusterProfiler::enrichGO()]. "BP", "MF", "CC", or "ALL". Defaults to "MF".
#' @param universe A `glyexp::experiment()` defining the background proteins, or `NULL` to use the default background.
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichGO()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichGO()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichGO()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#' - `org.Hs.eg.db` for human gene annotation or other OrgDb packages
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw clusterProfiler enrichResult object
#'  - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_go_ora_res` and `glystats_res`.
#' @seealso [clusterProfiler::enrichGO()]
#' @keywords internal
#' @export
gly_enrich_go <- function(
  exp,
  add_info = TRUE,
  orgdb = "org.Hs.eg.db",
  ont = "MF",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  lifecycle::deprecate_warn(
    "0.10.0",
    "gly_enrich_go()",
    "glyfun::enrich_ora_go()"
  )
  rlang::check_installed("clusterProfiler")
  checkmate::assert_logical(add_info, len = 1)
  proteins <- .extract_uniprot_from_exp(exp)
  result <- .analyze_enrich_go(
    proteins,
    orgdb = orgdb,
    ont = ont,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_enrich_go()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_enrich_go <- function(
  proteins,
  orgdb = "org.Hs.eg.db",
  keytype = "UNIPROT",
  ont = "MF",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  rlang::check_installed("clusterProfiler")

  # Validate and prepare arguments
  checkmate::assert_character(proteins, min.len = 1)
  universe <- .prepare_universe_uniprot(universe)

  # Perform analysis
  suppressMessages(
    raw_result <- clusterProfiler::enrichGO(
      proteins,
      OrgDb = orgdb,
      keyType = keytype,
      ont = ont,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      universe = universe,
      pAdjustMethod = p_adj_method,
      readable = TRUE
    )
  )

  .package_enrich_result(raw_result, "glystats_go_ora_res")
}

#' KEGG over-representation analysis (ORA)
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function was deprecated because we decided to move all enrichment analysis functions
#' to the separate `glyfun` package, which has more features and better API design.
#'
#' Perform KEGG ORA for protein UniProt accessions using [clusterProfiler::enrichKEGG()].
#' - `gly_enrich_kegg()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#' @param organism Passed to `organism` of [clusterProfiler::enrichKEGG()].
#'   KEGG organism code (e.g., "hsa" for human, "mmu" for mouse). Defaults to "hsa".
#' @param universe A `glyexp::experiment()` defining the background proteins, or `NULL` to use the default background.
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichKEGG()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichKEGG()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichKEGG()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID (e.g., hsa:XXXXX)
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw clusterProfiler enrichResult object
#'  - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_kegg_ora_res` and `glystats_res`.
#' @seealso [clusterProfiler::enrichKEGG()]
#' @keywords internal
#' @export
gly_enrich_kegg <- function(
  exp,
  add_info = TRUE,
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  lifecycle::deprecate_warn(
    "0.10.0",
    "gly_enrich_kegg()",
    "glyfun::enrich_ora_kegg()"
  )
  rlang::check_installed("clusterProfiler")
  checkmate::assert_logical(add_info, len = 1)
  proteins <- .extract_uniprot_from_exp(exp)
  result <- .analyze_enrich_kegg(
    proteins,
    organism = organism,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_enrich_kegg()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_enrich_kegg <- function(
  proteins,
  keytype = "uniprot",
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  rlang::check_installed("clusterProfiler")

  # Validate and prepare arguments
  checkmate::assert_character(proteins, min.len = 1)
  universe <- .prepare_universe_uniprot(universe)

  suppressMessages(
    raw_result <- clusterProfiler::enrichKEGG(
      proteins,
      organism = organism,
      keyType = keytype,
      universe = universe,
      pAdjustMethod = p_adj_method,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff
    )
  )

  .package_enrich_result(raw_result, "glystats_kegg_ora_res")
}

#' Reactome pathway over-representation analysis (ORA)
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function was deprecated because we decided to move all enrichment analysis functions
#' to the separate `glyfun` package, which has more features and better API design.
#'
#' Perform Reactome pathway ORA for protein UniProt accessions using [ReactomePA::enrichPathway()].
#' - `gly_enrich_reactome()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#'
#' As [ReactomePA::enrichPathway()] only accepts Entrez IDs,
#' the UniProt IDs will be first transformed into Entrez IDs with [clusterProfiler::bitr()].
#'
#' @param exp A `glyexp::experiment()` object.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::bitr()].
#'   Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to "org.Hs.eg.db".
#' @param organism Passed to `organism` of [ReactomePA::enrichPathway()].
#'   Species name (e.g., "human", "mouse", "rat"). Defaults to "human".
#' @param universe A `glyexp::experiment()` defining the background proteins, or `NULL` to use the default background.
#' @param p_adj_method Passed to `pAdjustMethod` of [ReactomePA::enrichPathway()].
#' @param p_cutoff Passed to `pvalueCutoff` of [ReactomePA::enrichPathway()].
#' @param q_cutoff Passed to `qvalueCutoff` of [ReactomePA::enrichPathway()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for ID conversion
#' - `ReactomePA` for enrichment analysis
#' - `org.Hs.eg.db` for human gene annotation or other OrgDb packages
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID (e.g., R-HSA-XXXXXX)
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw ReactomePA enrichPathway result object
#'  - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_reactome_ora_res` and `glystats_res`.
#' @seealso [ReactomePA::enrichPathway()]
#' @keywords internal
#' @export
gly_enrich_reactome <- function(
  exp,
  add_info = TRUE,
  orgdb = "org.Hs.eg.db",
  organism = "human",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  lifecycle::deprecate_warn(
    "0.10.0",
    "gly_enrich_reactome()",
    "glyfun::enrich_ora_reactome()"
  )
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("ReactomePA")
  checkmate::assert_logical(add_info, len = 1)

  proteins <- .extract_uniprot_from_exp(exp)
  result <- .analyze_enrich_reactome(
    proteins,
    orgdb = orgdb,
    organism = organism,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_enrich_reactome()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_enrich_reactome <- function(
  proteins,
  orgdb = "org.Hs.eg.db",
  organism = "human",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("ReactomePA")

  # Validate arguments
  checkmate::assert_character(proteins, min.len = 1)

  # Convert foreground proteins to Entrez IDs
  cli::cli_alert_info("Converting foreground proteins to Entrez IDs")
  fg_genes <- .uniprot_to_entrez(proteins, orgdb)

  # Handle universe if provided
  bg_genes <- .prepare_universe_entrez(universe, orgdb)

  # Perform Reactome pathway analysis
  suppressMessages(
    raw_result <- ReactomePA::enrichPathway(
      gene = fg_genes,
      organism = organism,
      universe = bg_genes,
      pAdjustMethod = p_adj_method,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      readable = TRUE
    )
  )

  .package_enrich_result(raw_result, "glystats_reactome_ora_res")
}

#' WikiPathways over-representation analysis (ORA)
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function was deprecated because we decided to move all enrichment analysis functions
#' to the separate `glyfun` package, which has more features and better API design.
#'
#' Perform WikiPathways ORA for protein UniProt accessions using [clusterProfiler::enrichWP()].
#' - `gly_enrich_wikipathways()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#'
#' As [clusterProfiler::enrichWP()] only accepts Entrez IDs,
#' the UniProt IDs will be first transformed into Entrez IDs with [clusterProfiler::bitr()].
#'
#' @param exp A `glyexp::experiment()` object.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#' @param organism Passed to `organism` of [clusterProfiler::enrichWP()].
#'   Species name (e.g., "Homo sapiens", "Mus musculus", "Rattus norvegicus"). Defaults to "Homo sapiens".
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::bitr()].
#'   Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to "org.Hs.eg.db".
#' @param universe A `glyexp::experiment()` defining the background proteins, or `NULL` to use the default background.
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichWP()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichWP()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichWP()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis and ID conversion
#' - `org.Hs.eg.db` for human gene annotation or other OrgDb packages
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID (e.g., WPXXXXX)
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw clusterProfiler enrichResult object
#'  - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_wikipathways_ora_res` and `glystats_res`.
#' @seealso [clusterProfiler::enrichWP()]
#' @keywords internal
#' @export
gly_enrich_wikipathways <- function(
  exp,
  add_info = TRUE,
  organism = "Homo sapiens",
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  lifecycle::deprecate_warn(
    "0.10.0",
    "gly_enrich_wikipathways()",
    "glyfun::enrich_ora_wp()"
  )
  rlang::check_installed("clusterProfiler")
  checkmate::assert_logical(add_info, len = 1)

  proteins <- .extract_uniprot_from_exp(exp)
  result <- .analyze_enrich_wikipathways(
    proteins,
    organism = organism,
    orgdb = orgdb,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_enrich_wikipathways()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_enrich_wikipathways <- function(
  proteins,
  organism = "Homo sapiens",
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  rlang::check_installed("clusterProfiler")

  # Validate arguments
  checkmate::assert_character(proteins, min.len = 1)

  # Convert foreground proteins to Entrez IDs
  cli::cli_alert_info("Converting foreground proteins to Entrez IDs")
  fg_genes <- .uniprot_to_entrez(proteins, orgdb)

  # Handle universe if provided
  bg_genes <- .prepare_universe_entrez(universe, orgdb)

  # Perform WikiPathways analysis
  suppressMessages(
    raw_result <- clusterProfiler::enrichWP(
      gene = fg_genes,
      organism = organism,
      universe = bg_genes,
      pAdjustMethod = p_adj_method,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff
    )
  )
  .package_enrich_result(raw_result, "glystats_wikipathways_ora_res")
}

#' Disease Ontology over-representation analysis (ORA)
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function was deprecated because we decided to move all enrichment analysis functions
#' to the separate `glyfun` package, which has more features and better API design.
#'
#' Perform Disease Ontology ORA for protein UniProt accessions using [DOSE::enrichDO()].
#' - `gly_enrich_do()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#'
#' As [DOSE::enrichDO()] only accepts Entrez IDs,
#' the UniProt IDs will be first transformed into Entrez IDs with [clusterProfiler::bitr()].
#'
#' @param exp A `glyexp::experiment()` object.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#' @param ont Passed to `ont` of [DOSE::enrichDO()].
#'   Disease Ontology type. Options: `"HDO"` (Human Disease Ontology), `"HPO"` (Human Phenotype Ontology),
#'   `"MPO"` (Mammalian Phenotype Ontology). Defaults to `"HDO"`.
#' @param organism Passed to `organism` of [DOSE::enrichDO()].
#'   Organism code. `"hsa"` for human (Homo sapiens), `"mmu"` for mouse (Mus musculus). Defaults to `"hsa"`.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::bitr()].
#'   Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to "org.Hs.eg.db".
#' @param universe A `glyexp::experiment()` defining the background proteins, or `NULL` to use the default background.
#' @param p_adj_method Passed to `pAdjustMethod` of [DOSE::enrichDO()].
#' @param p_cutoff Passed to `pvalueCutoff` of [DOSE::enrichDO()].
#' @param q_cutoff Passed to `qvalueCutoff` of [DOSE::enrichDO()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for ID conversion
#' - `DOSE` for enrichment analysis
#' - `org.Hs.eg.db` for human gene annotation or other OrgDb packages
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID (e.g., DOID:XXXX)
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw DOSE enrichResult object
#'  - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_do_ora_res` and `glystats_res`.
#' @seealso [DOSE::enrichDO()]
#' @keywords internal
#' @export
gly_enrich_do <- function(
  exp,
  add_info = TRUE,
  ont = "HDO",
  organism = "hsa",
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  lifecycle::deprecate_warn(
    "0.10.0",
    "gly_enrich_do()",
    "glyfun::enrich_ora_do()"
  )
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("DOSE")
  checkmate::assert_logical(add_info, len = 1)

  proteins <- .extract_uniprot_from_exp(exp)
  result <- .analyze_enrich_do(
    proteins,
    ont = ont,
    organism = organism,
    orgdb = orgdb,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_enrich_do()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_enrich_do <- function(
  proteins,
  ont = "HDO",
  organism = "hsa",
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("DOSE")

  # Validate arguments
  checkmate::assert_character(proteins, min.len = 1)

  # Convert foreground proteins to Entrez IDs
  cli::cli_alert_info("Converting foreground proteins to Entrez IDs")
  fg_genes <- .uniprot_to_entrez(proteins, orgdb)

  # Handle universe if provided
  bg_genes <- .prepare_universe_entrez(universe, orgdb)

  # Perform Disease Ontology analysis
  suppressMessages(
    raw_result <- DOSE::enrichDO(
      gene = fg_genes,
      ont = ont,
      organism = organism,
      universe = bg_genes,
      pAdjustMethod = p_adj_method,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      readable = TRUE
    )
  )
  .package_enrich_result(raw_result, "glystats_do_ora_res")
}

#' Network of Cancer Genes (NCG) over-representation analysis (ORA)
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function was deprecated because we decided to move all enrichment analysis functions
#' to the separate `glyfun` package, which has more features and better API design.
#'
#' Perform NCG ORA for protein UniProt accessions using [DOSE::enrichNCG()].
#' - `gly_enrich_ncg()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#'
#' As [DOSE::enrichNCG()] only accepts Entrez IDs,
#' the UniProt IDs will be first transformed into Entrez IDs with [clusterProfiler::bitr()].
#'
#' @param exp A `glyexp::experiment()` object.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::bitr()].
#'   Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to "org.Hs.eg.db".
#' @param universe A `glyexp::experiment()` defining the background proteins, or `NULL` to use the default background.
#' @param p_adj_method Passed to `pAdjustMethod` of [DOSE::enrichNCG()].
#' @param p_cutoff Passed to `pvalueCutoff` of [DOSE::enrichNCG()].
#' @param q_cutoff Passed to `qvalueCutoff` of [DOSE::enrichNCG()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for ID conversion
#' - `DOSE` for enrichment analysis
#' - `org.Hs.eg.db` for human gene annotation or other OrgDb packages
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw DOSE enrichResult object
#'  - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_ncg_ora_res` and `glystats_res`.
#' @seealso [DOSE::enrichNCG()]
#' @keywords internal
#' @export
gly_enrich_ncg <- function(
  exp,
  add_info = TRUE,
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  lifecycle::deprecate_warn(
    "0.10.0",
    "gly_enrich_ncg()",
    "glyfun::enrich_ora_ncg()"
  )
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("DOSE")
  checkmate::assert_logical(add_info, len = 1)

  proteins <- .extract_uniprot_from_exp(exp)
  result <- .analyze_enrich_ncg(
    proteins,
    orgdb = orgdb,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_enrich_ncg()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_enrich_ncg <- function(
  proteins,
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("DOSE")

  # Validate arguments
  checkmate::assert_character(proteins, min.len = 1)

  # Convert foreground proteins to Entrez IDs
  cli::cli_alert_info("Converting foreground proteins to Entrez IDs")
  fg_genes <- .uniprot_to_entrez(proteins, orgdb)

  # Handle universe if provided
  bg_genes <- .prepare_universe_entrez(universe, orgdb)

  # Perform NCG analysis
  suppressMessages(
    raw_result <- DOSE::enrichNCG(
      gene = fg_genes,
      universe = bg_genes,
      pAdjustMethod = p_adj_method,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      readable = TRUE
    )
  )
  .package_enrich_result(raw_result, "glystats_ncg_ora_res")
}

#' DisGeNET over-representation analysis (ORA)
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function was deprecated because we decided to move all enrichment analysis functions
#' to the separate `glyfun` package, which has more features and better API design.
#'
#' Perform DisGeNET ORA for protein UniProt accessions using [DOSE::enrichDGN()].
#' - `gly_enrich_dgn()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#'
#' As [DOSE::enrichDGN()] only accepts Entrez IDs,
#' the UniProt IDs will be first transformed into Entrez IDs with [clusterProfiler::bitr()].
#'
#' @param exp A `glyexp::experiment()` object.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::bitr()].
#'   Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to "org.Hs.eg.db".
#' @param universe A `glyexp::experiment()` defining the background proteins, or `NULL` to use the default background.
#' @param p_adj_method Passed to `pAdjustMethod` of [DOSE::enrichDGN()].
#' @param p_cutoff Passed to `pvalueCutoff` of [DOSE::enrichDGN()].
#' @param q_cutoff Passed to `qvalueCutoff` of [DOSE::enrichDGN()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for ID conversion
#' - `DOSE` for enrichment analysis
#' - `org.Hs.eg.db` for human gene annotation or other OrgDb packages
#'
#' @return A list with three elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw DOSE enrichResult object
#'  - `meta_data`: A list containing metadata from the input experiment
#' The list has classes `glystats_dgn_ora_res` and `glystats_res`.
#' @seealso [DOSE::enrichDGN()]
#' @keywords internal
#' @export
gly_enrich_dgn <- function(
  exp,
  add_info = TRUE,
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  lifecycle::deprecate_warn("0.10.0", "gly_enrich_dgn()")
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("DOSE")
  checkmate::assert_logical(add_info, len = 1)

  proteins <- .extract_uniprot_from_exp(exp)
  result <- .analyze_enrich_dgn(
    proteins,
    orgdb = orgdb,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )

  # Add meta_data from experiment
  result$meta_data <- glyexp::get_meta_data(exp)

  result
}

#' Run the internal computation for `gly_enrich_dgn()`
#'
#' This helper receives components extracted from a validated experiment.
#'
#' @noRd
.analyze_enrich_dgn <- function(
  proteins,
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("DOSE")

  # Validate arguments
  checkmate::assert_character(proteins, min.len = 1)

  # Convert foreground proteins to Entrez IDs
  cli::cli_alert_info("Converting foreground proteins to Entrez IDs")
  fg_genes <- .uniprot_to_entrez(proteins, orgdb)

  # Handle universe if provided
  bg_genes <- .prepare_universe_entrez(universe, orgdb)

  # Perform DisGeNET analysis
  suppressMessages(
    raw_result <- DOSE::enrichDGN(
      gene = fg_genes,
      universe = bg_genes,
      pAdjustMethod = p_adj_method,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      readable = TRUE
    )
  )
  .package_enrich_result(raw_result, "glystats_dgn_ora_res")
}

#' @param exp The experiment.
#' @returns A character vector of UniProt IDs.
#' @noRd
.extract_uniprot_from_exp <- function(exp) {
  checkmate::assert_class(exp, "glyexp_experiment")
  var_info <- glyexp::get_var_info(exp)
  if ("protein" %in% colnames(var_info)) {
    uniprot <- var_info$protein
  } else {
    cli::cli_abort(
      "{.field protein} column not found in the variable information tibble."
    )
  }
  unique(uniprot[!is.na(uniprot)])
}

#' Convert UniProt IDs to Entrez IDs
#' @param uniprot The UniProt IDs.
#' @param orgdb An OrgDb object.
#' @returns The transformed Entrez IDs.
#' @noRd
.uniprot_to_entrez <- function(uniprot, orgdb) {
  suppressWarnings(suppressMessages(
    entrez_ids <- clusterProfiler::bitr(
      uniprot,
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = orgdb
    )$ENTREZID
  ))
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  n_failed <- length(uniprot) - length(entrez_ids)
  if (n_failed > 0) {
    pct_failed <- round(n_failed / length(uniprot) * 100, 1)
    cli::cli_alert_warning(
      "{.val {n_failed}} of {.val {length(uniprot)}} ({.val {pct_failed}}%) proteins failed to map to Entrez IDs."
    )
  }
  entrez_ids
}

#' Tidy and package enrichment result into standard format
#'
#' @param raw_result The raw enrichment result from clusterProfiler or ReactomePA.
#' @param result_class The specific result class to assign (e.g., "glystats_go_ora_res").
#'
#' @returns A structured list with tidy_result and raw_result.
#' @noRd
.package_enrich_result <- function(raw_result, result_class) {
  # Tidy result
  tidy_result <- tibble::as_tibble(raw_result)
  tidy_result <- janitor::clean_names(tidy_result)

  # Rename p_value to p_val and p_adjust to p_adj for consistency
  if ("p_value" %in% colnames(tidy_result)) {
    tidy_result <- dplyr::rename(tidy_result, p_val = "p_value")
  }
  if ("p_adjust" %in% colnames(tidy_result)) {
    tidy_result <- dplyr::rename(tidy_result, p_adj = "p_adjust")
  }

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result
    ),
    class = c(result_class, "glystats_res")
  )
}

#' Prepare universe for enrichment analysis (UniProt path)
#'
#' Extracts UniProt IDs from an experiment.
#'
#' @param universe A glyexp experiment or `NULL`.
#' @returns Character vector of UniProt IDs or NULL.
#' @noRd
.prepare_universe_uniprot <- function(universe) {
  if (is.null(universe)) {
    return(NULL)
  }
  .extract_uniprot_from_exp(universe)
}

#' Prepare universe for enrichment analysis (Entrez path)
#'
#' Extracts UniProt IDs from an experiment, then converts them to Entrez IDs.
#'
#' @param universe A glyexp experiment or `NULL`.
#' @param orgdb An OrgDb object for ID conversion.
#' @returns Character vector of Entrez IDs or NULL.
#' @noRd
.prepare_universe_entrez <- function(universe, orgdb) {
  if (is.null(universe)) {
    return(NULL)
  }
  universe <- .extract_uniprot_from_exp(universe)
  cli::cli_alert_info("Converting background proteins to Entrez IDs")
  .uniprot_to_entrez(universe, orgdb)
}
