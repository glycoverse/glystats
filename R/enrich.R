#' GO, KEGG, and Reactome over-representation analysis (ORA)
#'
#' @description
#' Perform GO, KEGG, and Reactome ORA for all proteins (genes) in a `glyexp::experiment()`.
#'
#' This function uses the "protein" column in the variable information tibble.
#' Protein identifiers should be UniProt accessions.
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#' - `ReactomePA` for Reactome pathway analysis
#' - `org.Hs.eg.db` for human gene annotation (GO analysis only)
#' 
#' @param exp A `glyexp::experiment()` object.
#' @param return_raw A logical value indicating whether to return raw clusterProfiler enrichResult objects.
#' @param ... Additional arguments passed to `clusterProfiler::enrichGO()`, `clusterProfiler::enrichKEGG()`, or `ReactomePA::enrichPathway()`.
#'
#' @return A tibble of GO, KEGG, or Reactome enrichment results (when return_raw = FALSE),
#'   or raw clusterProfiler enrichResult objects (when return_raw = TRUE).
#' @seealso [clusterProfiler::enrichGO()], [clusterProfiler::enrichKEGG()], [ReactomePA::enrichPathway()]
#' @export
enrich_go <- function(exp, return_raw = FALSE, ...) {
  .check_pkg_available("clusterProfiler")
  .check_pkg_available("org.Hs.eg.db")
  
  checkmate::check_logical(return_raw, len = 1)
  
  genes <- .extract_genes_from_exp(exp)
  res <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "UNIPROT",
    readable = TRUE,
    ...
  )
  
  # Return raw results if requested
  if (return_raw) {
    return(res)
  }
  
  res <- tibble::as_tibble(res)
  res <- janitor::clean_names(res)
  structure(res, class = c("glystats_go_ora_res", class(res)))
}

#' @rdname enrich_go
#' @export
enrich_kegg <- function(exp, return_raw = FALSE, ...) {
  .check_pkg_available("clusterProfiler")
  
  checkmate::check_logical(return_raw, len = 1)
  
  genes <- .extract_genes_from_exp(exp)
  res <- clusterProfiler::enrichKEGG(
    gene = genes,
    keyType = "uniprot",
    ...
  )
  
  # Return raw results if requested
  if (return_raw) {
    return(res)
  }
  
  res <- tibble::as_tibble(res)
  res <- janitor::clean_names(res)
  structure(res, class = c("glystats_kegg_ora_res", class(res)))
}

#' @rdname enrich_go
#' @export
enrich_reactome <- function(exp, return_raw = FALSE, ...) {
  .check_pkg_available("clusterProfiler")
  .check_pkg_available("ReactomePA")
  .check_pkg_available("org.Hs.eg.db")
  
  checkmate::check_logical(return_raw, len = 1)
  
  uniprot_ids <- .extract_genes_from_exp(exp)
  
  # Convert UniProt to Entrez IDs
  suppressWarnings(suppressMessages(
    entrez_ids <- clusterProfiler::bitr(
      uniprot_ids,
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )$ENTREZID
  ))
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  n_failed <- length(uniprot_ids) - length(entrez_ids)
  if (n_failed > 0) {
    pct_failed <- round(n_failed / length(uniprot_ids) * 100, 1)
    cli::cli_alert_warning("{.val {n_failed}} of {.val {length(uniprot_ids)}} ({.val {pct_failed}}%) proteins failed to map to Entrez IDs.")
  }
  
  # Perform Reactome pathway analysis
  res <- ReactomePA::enrichPathway(
    gene = entrez_ids,
    organism = "human",
    readable = TRUE,
    ...
  )
  
  # Return raw results if requested
  if (return_raw) {
    return(res)
  }
  
  res <- tibble::as_tibble(res)
  res <- janitor::clean_names(res)
  structure(res, class = c("glystats_reactome_ora_res", class(res)))
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
