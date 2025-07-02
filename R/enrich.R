#' GO and KEGG over-representation analysis (ORA)
#'
#' @description
#' Perform GO and KEGG ORA for all proteins (genes) in a `glyexp::experiment()`.
#'
#' This function uses the "protein" column in the variable information tibble, or
#' falls back to the "proteins" column if "protein" is not present.
#' If the "proteins" column is used, [glyclean::infer_protein()] is called to resolve
#' cases with multiple gene entries.
#' In both cases, protein identifiers should be UniProt accessions.
#' For the "proteins" column, multiple proteins should be separated by ";".
#'
#' @section Required packages:
#' This function requires the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#' - `org.Hs.eg.db` for human gene annotation (GO analysis only)
#' 
#' @param exp A `glyexp::experiment()` object.
#' @param return_raw A logical value indicating whether to return raw clusterProfiler enrichResult objects.
#' @param ... Additional arguments passed to `clusterProfiler::enrichGO()` or `clusterProfiler::enrichKEGG()`.
#'
#' @return A tibble of GO or KEGG enrichment results (when return_raw = FALSE),
#'   or raw clusterProfiler enrichResult objects (when return_raw = TRUE).
#' @seealso [clusterProfiler::enrichGO()], [clusterProfiler::enrichKEGG()]
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

# Extract genes from experiment object helper function
.extract_genes_from_exp <- function(exp) {
  if ("protein" %in% colnames(exp$var_info)) {
    genes <- exp$var_info$protein  # Uniprot
  } else if ("proteins" %in% colnames(exp$var_info)) {
    exp2 <- glyclean::infer_protein(exp)  # resolve multiple gene problem
    genes <- exp2$var_info$protein
  } else {
    cli::cli_abort("{.field protein} or {.field proteins} column not found in the variable information tibble.")
  }
  unique(genes[!is.na(genes)])
}
