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
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw clusterProfiler enrichResult object
#' The list has classes `glystats_go_ora_res`/`glystats_kegg_ora_res`/`glystats_reactome_ora_res` and `glystats_res`.
#' @seealso [clusterProfiler::enrichGO()], [clusterProfiler::enrichKEGG()], [ReactomePA::enrichPathway()]
#' @export
gly_enrich_go <- function(exp, add_info = TRUE, ...) {
  rlang::check_installed("clusterProfiler")

  checkmate::assert_logical(add_info, len = 1)

  genes <- .extract_genes_from_exp(exp)
  dots <- rlang::list2(...)
  rlang::exec(gly_enrich_go_, genes, !!!dots)
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_go_ <- function(proteins, ...) {
  rlang::check_installed("clusterProfiler")

  checkmate::assert_character(proteins, min.len = 1)

  dots <- rlang::list2(...)
  call_args <- c(list(gene = proteins), dots)
  if (!"OrgDb" %in% names(call_args)) {
    rlang::check_installed("org.Hs.eg.db")
    call_args$OrgDb <- "org.Hs.eg.db"
  }
  if (!"keyType" %in% names(call_args)) {
    call_args$keyType <- "UNIPROT"
  }
  if (!"readable" %in% names(call_args)) {
    call_args$readable <- TRUE
  }
  if ("gene" %in% names(dots)) {
    cli::cli_abort("{.field gene} should not be supplied through `...`; use the function's first argument instead.")
  }

  suppressMessages(
    raw_result <- do.call(clusterProfiler::enrichGO, call_args)
  )

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
    class = c("glystats_go_ora_res", "glystats_res")
  )
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_kegg <- function(exp, add_info = TRUE, ...) {
  rlang::check_installed("clusterProfiler")

  checkmate::assert_logical(add_info, len = 1)

  genes <- .extract_genes_from_exp(exp)
  dots <- rlang::list2(...)
  rlang::exec(gly_enrich_kegg_, genes, !!!dots)
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_kegg_ <- function(proteins, ...) {
  rlang::check_installed("clusterProfiler")

  checkmate::assert_character(proteins, min.len = 1)

  dots <- rlang::list2(...)
  if ("gene" %in% names(dots)) {
    cli::cli_abort("{.field gene} should not be supplied through `...`; use the function's first argument instead.")
  }
  call_args <- c(list(gene = proteins), dots)
  if (!"keyType" %in% names(call_args)) {
    call_args$keyType <- "uniprot"
  }

  suppressMessages(
    raw_result <- do.call(clusterProfiler::enrichKEGG, call_args)
  )

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
    class = c("glystats_kegg_ora_res", "glystats_res")
  )
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_reactome <- function(exp, add_info = TRUE, ...) {
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("ReactomePA")

  checkmate::assert_logical(add_info, len = 1)

  dots <- rlang::list2(...)
  if (!"OrgDb" %in% names(dots)) {
    rlang::check_installed("org.Hs.eg.db")
  }

  uniprot_ids <- .extract_genes_from_exp(exp)
  rlang::exec(gly_enrich_reactome_, uniprot_ids, !!!dots)
}

#' @rdname gly_enrich_go
#' @export
gly_enrich_reactome_ <- function(proteins, ...) {
  rlang::check_installed("clusterProfiler")
  rlang::check_installed("ReactomePA")

  checkmate::assert_character(proteins, min.len = 1)

  dots <- rlang::list2(...)
  if ("gene" %in% names(dots)) {
    cli::cli_abort("{.field gene} should not be supplied through `...`; use the function's first argument instead.")
  }

  bitr_arg_names <- c("OrgDb", "fromType", "toType", "drop")
  bitr_args <- dots[intersect(names(dots), bitr_arg_names)]
  dots <- dots[setdiff(names(dots), bitr_arg_names)]

  orgdb <- bitr_args$OrgDb
  if (is.null(orgdb)) {
    rlang::check_installed("org.Hs.eg.db")
    orgdb <- "org.Hs.eg.db"
  }
  orgdb_obj <- .resolve_orgdb(orgdb)

  from_type <- if (!is.null(bitr_args$fromType)) bitr_args$fromType else "UNIPROT"
  to_type <- if (!is.null(bitr_args$toType)) bitr_args$toType else "ENTREZID"
  drop_arg <- bitr_args$drop

  # Convert UniProt to Entrez IDs
  bitr_call <- list(
    proteins,
    fromType = from_type,
    toType = to_type,
    OrgDb = orgdb_obj
  )
  if (!is.null(drop_arg)) {
    bitr_call$drop <- drop_arg
  }

  suppressWarnings(suppressMessages(
    entrez_ids <- do.call(clusterProfiler::bitr, bitr_call)$ENTREZID
  ))
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  n_failed <- length(proteins) - length(entrez_ids)
  if (n_failed > 0) {
    pct_failed <- round(n_failed / length(proteins) * 100, 1)
    cli::cli_alert_warning("{.val {n_failed}} of {.val {length(proteins)}} ({.val {pct_failed}}%) proteins failed to map to Entrez IDs.")
  }

  # Perform Reactome pathway analysis
  call_args <- c(list(gene = entrez_ids), dots)
  if (!"organism" %in% names(call_args)) {
    call_args$organism <- "human"
  }
  if (!"readable" %in% names(call_args)) {
    call_args$readable <- TRUE
  }
  suppressMessages(
    raw_result <- do.call(ReactomePA::enrichPathway, call_args)
  )

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

.resolve_orgdb <- function(orgdb) {
  if (is.character(orgdb) && length(orgdb) == 1) {
    rlang::check_installed(orgdb)
    get(orgdb, envir = asNamespace(orgdb))
  } else {
    orgdb
  }
}
