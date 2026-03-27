#' GO over-representation analysis (ORA)
#'
#' @description
#' Perform GO ORA for protein UniProt accessions using [clusterProfiler::enrichGO()].
#' - `gly_enrich_go()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#' - `gly_enrich_go_()` accepts a character vector of Uniprot IDs.
#'
#' @param exp (Only for [gly_enrich_go()]) A `glyexp::experiment()` object.
#' @param proteins (Only for [gly_enrich_go_()]) A character vector of UniProt accession IDs.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#'  Only applicable to top-level APIs.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::enrichGO()].
#' @param keytype Passed to `keyType` of [clusterProfiler::enrichGO()].
#' @param ont Passed to `ont` of [clusterProfiler::enrichGO()]. "BP", "MF", "CC", or "ALL". Defaults to "MF".
#' @param universe Background genes. If a character vector, directly passed to `universe` of [clusterProfiler::enrichGO()].
#'   You can also provide a [glyexp::experiment()] object with "glycoproteomics" type.
#'   In this case all detected proteins in this experiment will be extracted and passed to
#'   `universe` of [clusterProfiler::enrichGO()].
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichGO()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichGO()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichGO()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#' - `org.Hs.eg.db` for human gene annotation or other OrgDb packages
#'
#' @return A list with two elements:
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
#' The list has classes `glystats_go_ora_res` and `glystats_res`.
#' @seealso [clusterProfiler::enrichGO()]
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
  rlang::check_installed("clusterProfiler")
  checkmate::assert_logical(add_info, len = 1)
  proteins <- .extract_uniprot_from_exp(exp)
  gly_enrich_go_(
    proteins,
    orgdb = orgdb,
    ont = ont,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_go_ <- function(
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
#' Perform KEGG ORA for protein UniProt accessions using [clusterProfiler::enrichKEGG()].
#' - `gly_enrich_kegg()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#' - `gly_enrich_kegg_()` accepts a character vector of Uniprot IDs.
#'
#' @param exp (Only for [gly_enrich_kegg()]) A `glyexp::experiment()` object.
#' @param proteins (Only for [gly_enrich_kegg_()]) A character vector of UniProt accession IDs.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#'  Only applicable to top-level APIs.
#' @param organism Passed to `organism` of [clusterProfiler::enrichKEGG()].
#'   KEGG organism code (e.g., "hsa" for human, "mmu" for mouse). Defaults to "hsa".
#' @param keytype Passed to `keyType` of [clusterProfiler::enrichKEGG()].
#'   Defaults to "uniprot".
#' @param universe Background genes. If a character vector, directly passed to `universe` of [clusterProfiler::enrichKEGG()].
#'   You can also provide a [glyexp::experiment()] object with "glycoproteomics" type.
#'   In this case all detected proteins in this experiment will be extracted and passed to
#'   `universe` of [clusterProfiler::enrichKEGG()].
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichKEGG()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichKEGG()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichKEGG()].
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#'
#' @return A list with two elements:
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
#' The list has classes `glystats_kegg_ora_res` and `glystats_res`.
#' @seealso [clusterProfiler::enrichKEGG()]
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
  rlang::check_installed("clusterProfiler")
  checkmate::assert_logical(add_info, len = 1)
  proteins <- .extract_uniprot_from_exp(exp)
  gly_enrich_kegg_(
    proteins,
    organism = organism,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )
}

#' @rdname gly_enrich_kegg
#' @export
gly_enrich_kegg_ <- function(
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
#' Perform Reactome pathway ORA for protein UniProt accessions using [ReactomePA::enrichPathway()].
#' - `gly_enrich_reactome()` accepts a [glyexp::experiment()] and extracts protein information
#' from the "protein" column in the variable information tibble.
#' - `gly_enrich_reactome_()` accepts a character vector of Uniprot IDs.
#'
#' As [ReactomePA::enrichPathway()] only accepts Entrez IDs,
#' the Uniprot IDs will be firsted transformed into Entrez IDs with [clusterProfiler::bitr()].
#'
#' @param exp (Only for [gly_enrich_reactome()]) A `glyexp::experiment()` object.
#' @param proteins (Only for [gly_enrich_reactome_()]) A character vector of UniProt accession IDs.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#'  Only applicable to top-level APIs.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::bitr()].
#'   Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to "org.Hs.eg.db".
#' @param organism Passed to `organism` of [ReactomePA::enrichPathway()].
#'   Species name (e.g., "human", "mouse", "rat"). Defaults to "human".
#' @param universe Background genes. If a character vector, directly passed to `universe` of [ReactomePA::enrichPathway()].
#'   You can also provide a [glyexp::experiment()] object with "glycoproteomics" type.
#'   In this case all detected proteins in this experiment will be extracted and passed to
#'   `universe` of [ReactomePA::enrichPathway()].
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
#' @return A list with two elements:
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
#' The list has classes `glystats_reactome_ora_res` and `glystats_res`.
#' @seealso [ReactomePA::enrichPathway()]
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
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("ReactomePA")
  checkmate::assert_logical(add_info, len = 1)

  proteins <- .extract_uniprot_from_exp(exp)
  gly_enrich_reactome_(
    proteins,
    orgdb = orgdb,
    organism = organism,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )
}

#' @rdname gly_enrich_reactome
#' @export
gly_enrich_reactome_ <- function(
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

#' Extract UniProt accessions from experiment object helper function
#' @param exp The experiment.
#' @returns A character vector of Uniprot IDs.
#' @noRd
.extract_uniprot_from_exp <- function(exp) {
  if ("protein" %in% colnames(exp$var_info)) {
    uniprot <- exp$var_info$protein
  } else {
    cli::cli_abort(
      "{.field protein} column not found in the variable information tibble."
    )
  }
  unique(uniprot[!is.na(uniprot)])
}

#' Convert Uniprot IDs to Entrez IDs
#' @param uniprot The Uniprot IDs.
#' @param orgdb An OrgDb object.
#' @returns The transformed Ectrez IDs.
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
#' Extracts UniProt IDs from experiment if needed.
#'
#' @param universe Background genes. Can be NULL, a character vector, or a glyexp experiment.
#' @returns Character vector of UniProt IDs or NULL.
#' @noRd
.prepare_universe_uniprot <- function(universe) {
  if (is.null(universe)) {
    return(NULL)
  }
  if (glyexp::is_experiment(universe)) {
    universe <- .extract_uniprot_from_exp(universe)
  }
  universe
}

#' Prepare universe for enrichment analysis (Entrez path)
#'
#' Extracts UniProt IDs from experiment if needed, then converts to Entrez IDs.
#'
#' @param universe Background genes. Can be NULL, a character vector of UniProt IDs, or a glyexp experiment.
#' @param orgdb An OrgDb object for ID conversion.
#' @returns Character vector of Entrez IDs or NULL.
#' @noRd
.prepare_universe_entrez <- function(universe, orgdb) {
  if (is.null(universe)) {
    return(NULL)
  }
  if (glyexp::is_experiment(universe)) {
    universe <- .extract_uniprot_from_exp(universe)
  }
  cli::cli_alert_info("Converting background proteins to Entrez IDs")
  .uniprot_to_entrez(universe, orgdb)
}
