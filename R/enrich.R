#' GO and KEGG enrichment analysis
#'
#' Perform GO and KEGG enrichment analysis based on all the proteins (genes)
#' in a `glyexp::experiment()`.
#'
#' @param exp A `glyexp::experiment()` object.
#' @param ... Arguments passed to [clusterProfiler::enrichGO()] or [clusterProfiler::enrichKEGG()].
#'
#' @return A tibble of GO or KEGG enrichment results.
#' @seealso [clusterProfiler::enrichGO()], [clusterProfiler::enrichKEGG()]
#' @export
enrich_go <- function(exp, ...) {
  genes <- .extract_genes_from_exp(exp)
  res <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "UNIPROT",
    readable = TRUE,
    ...
  )
  res <- tibble::as_tibble(res)
  res <- janitor::clean_names(res)
  structure(res, class = c("glystats_go_ora_res", class(res)))
}

#' @rdname enrich_go
#' @export
enrich_kegg <- function(exp, ...) {
  genes <- .extract_genes_from_exp(exp)
  res <- clusterProfiler::enrichKEGG(
    gene = genes,
    keyType = "uniprot",
    ...
  )
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
