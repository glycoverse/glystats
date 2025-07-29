#' GO, KEGG, and Reactome over-representation analysis (ORA)
#'
#' @description
#' Perform GO, KEGG, and Reactome ORA for proteins/genes.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param proteins A character vector of UniProt accession IDs.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#'  Only applicable to top-level APIs.
#' @param ... Additional arguments passed to `clusterProfiler::enrichGO()`, `clusterProfiler::enrichKEGG()`, or `ReactomePA::enrichPathway()`.
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#' - `ReactomePA` for Reactome pathway analysis
#' - `org.Hs.eg.db` for human gene annotation (GO analysis only)
#'
#' @details
#' These functions perform over-representation analysis using the specified database.
#'
#' `gly_enrich_go()`, `gly_enrich_kegg()`, and `gly_enrich_reactome()` are the top-level APIs
#' that work with `glyexp::experiment()` objects and extract protein information automatically
#' from the "protein" column in the variable information tibble.
#'
#' `gly_enrich_go_()`, `gly_enrich_kegg_()`, and `gly_enrich_reactome_()` are the underlying APIs
#' that work with protein vectors directly, providing more flexibility for users who don't use the glyexp package.
#'
#' **Gene Extraction (top-level APIs only):**
#' Proteins are extracted from the experiment's variable information. The function
#' looks for columns containing protein identifiers and uses them for enrichment analysis.
#' Protein identifiers should be UniProt accessions.
#'
#' **GO Analysis:**
#' Uses `clusterProfiler::enrichGO()` with UniProt IDs as input.
#'
#' **KEGG Analysis:**
#' Uses `clusterProfiler::enrichKEGG()` with UniProt IDs as input.
#'
#' **Reactome Analysis:**
#' Converts UniProt IDs to Entrez IDs and uses `ReactomePA::enrichPathway()`.
#'
#' @return A list with two elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: Term ID (GO:XXXXXXX, hsa:XXXXX, or R-HSA-XXXXX)
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_value`: Raw p-value from hypergeometric test
#'    - `p_adjust`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw clusterProfiler enrichResult object
#' The list has classes `glystats_go_ora_res`/`glystats_kegg_ora_res`/`glystats_reactome_ora_res` and `glystats_res`.
#' @seealso [clusterProfiler::enrichGO()], [clusterProfiler::enrichKEGG()], [ReactomePA::enrichPathway()]
#' @export
gly_enrich_go <- function(exp, add_info = TRUE, ...) {
  .check_pkg_available("clusterProfiler")
  .check_pkg_available("org.Hs.eg.db")

  checkmate::assert_logical(add_info, len = 1)

  genes <- .extract_genes_from_exp(exp)
  gly_enrich_go_(genes, ...)
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_go_ <- function(proteins, ...) {
  .check_pkg_available("clusterProfiler")
  .check_pkg_available("org.Hs.eg.db")

  checkmate::assert_character(proteins, min.len = 1)

  raw_result <- clusterProfiler::enrichGO(
    gene = proteins,
    OrgDb = "org.Hs.eg.db",
    keyType = "UNIPROT",
    readable = TRUE,
    ...
  )

  tidy_result <- tibble::as_tibble(raw_result)
  tidy_result <- janitor::clean_names(tidy_result)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result
    ),
    class = c("glystats_go_ora_res", "glystats_res")
  )
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_kegg <- function(exp, add_info = TRUE, ...) {
  .check_pkg_available("clusterProfiler")

  checkmate::assert_logical(add_info, len = 1)

  genes <- .extract_genes_from_exp(exp)
  gly_enrich_kegg_(genes, ...)
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_kegg_ <- function(proteins, ...) {
  .check_pkg_available("clusterProfiler")

  checkmate::assert_character(proteins, min.len = 1)

  raw_result <- clusterProfiler::enrichKEGG(
    gene = proteins,
    keyType = "uniprot",
    ...
  )

  tidy_result <- tibble::as_tibble(raw_result)
  tidy_result <- janitor::clean_names(tidy_result)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result
    ),
    class = c("glystats_kegg_ora_res", "glystats_res")
  )
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_reactome <- function(exp, add_info = TRUE, ...) {
  .check_pkg_available("clusterProfiler")
  .check_pkg_available("ReactomePA")
  .check_pkg_available("org.Hs.eg.db")

  checkmate::assert_logical(add_info, len = 1)

  uniprot_ids <- .extract_genes_from_exp(exp)
  gly_enrich_reactome_(uniprot_ids, ...)
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_reactome_ <- function(proteins, ...) {
  .check_pkg_available("clusterProfiler")
  .check_pkg_available("ReactomePA")
  .check_pkg_available("org.Hs.eg.db")

  checkmate::assert_character(proteins, min.len = 1)

  # Convert UniProt to Entrez IDs
  suppressWarnings(suppressMessages(
    entrez_ids <- clusterProfiler::bitr(
      proteins,
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )$ENTREZID
  ))
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  n_failed <- length(proteins) - length(entrez_ids)
  if (n_failed > 0) {
    pct_failed <- round(n_failed / length(proteins) * 100, 1)
    cli::cli_alert_warning("{.val {n_failed}} of {.val {length(proteins)}} ({.val {pct_failed}}%) proteins failed to map to Entrez IDs.")
  }

  # Perform Reactome pathway analysis
  raw_result <- ReactomePA::enrichPathway(
    gene = entrez_ids,
    organism = "human",
    readable = TRUE,
    ...
  )

  tidy_result <- tibble::as_tibble(raw_result)
  tidy_result <- janitor::clean_names(tidy_result)

  # Return list with both tidy and raw results
  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result
    ),
    class = c("glystats_reactome_ora_res", "glystats_res")
  )
}

# Extract genes from experiment object helper function
.extract_genes_from_exp <- function(exp) {
  if ("protein" %in% colnames(exp$var_info)) {
    genes <- exp$var_info$protein  # Uniprot
  } else {
    cli::cli_abort("{.field protein} column not found in the variable information tibble.")
  }
  unique(genes[!is.na(genes)])
}
