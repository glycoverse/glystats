#' GO, KEGG, and Reactome over-representation analysis (ORA)
#'
#' Perform GO, KEGG, and Reactome ORA for protein UniProt accessions.
#' For glycoproteomics experiments, the function extracts unique proteins from the variable information.
#'
#' @param exp (Only for [gly_enrich_go()], [gly_enrich_kegg()], [gly_enrich_reactome()]) A `glyexp::experiment()` object.
#' @param proteins (Only for [gly_enrich_go_()], [gly_enrich_kegg_()], [gly_enrich_reactome_()]) A character vector of UniProt accession IDs.
#' @param add_info A logical value. This parameter is included for API consistency but has no effect
#'  since enrichment results do not contain variable or sample columns.
#'  Only applicable to top-level APIs.
#' @param ... Additional arguments passed to `clusterProfiler::enrichGO()`, `clusterProfiler::enrichKEGG()`, or `ReactomePA::enrichPathway()`.
#'  See the "Additional arguments" section for more information.
#'
#' @section Required packages:
#' These functions require the following packages to be installed:
#' - `clusterProfiler` for enrichment analysis
#' - `ReactomePA` for Reactome pathway analysis
#' - `org.Hs.eg.db` for human gene annotation (GO analysis only)
#'
#' @details
#' # Details
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
#' # Additional arguments
#'
#' `universe` can be passed as a character vector of uniprot accession,
#' or a [glyexp::experiment()] object.
#' If latter, the background proteins will be extracted from the experiment.
#' This can be convenient when you first perform differential analysis and then perform enrichment
#' on the differentially expressed proteins:
#'
#' ```r
#' dea_res <- gly_limma(exp)
#' sig_exp <- filter_sig_vars(exp, dea_res)
#' enrich_res <- gly_enrich_go(sig_exp, universe = exp)
#' ```
#'
#' For `gly_enrich_reactome()`, an `OrgDb` can be passed through `...` to `clusterProfiler::bitr()`
#' to convert UniProt to Entrez IDs.
#' By default, `org.Hs.eg.db` is used.
#' This is useful when you want to use a different organism than human:
#'
#' ```r
#' enrich_res <- gly_enrich_reactome(exp, OrgDb = "org.Mm.eg.db")
#' ```
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

  genes <- .extract_uniprot_from_exp(exp)
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
  if ("universe" %in% names(dots)) {
    if (glyexp::is_experiment(dots$universe)) {
      call_args$universe <- .extract_uniprot_from_exp(dots$universe)
    }
  }
  if ("gene" %in% names(dots)) {
    cli::cli_abort(
      "{.field gene} should not be supplied through `...`; use the function's first argument instead."
    )
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

  genes <- .extract_uniprot_from_exp(exp)
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
    cli::cli_abort(
      "{.field gene} should not be supplied through `...`; use the function's first argument instead."
    )
  }
  call_args <- c(list(gene = proteins), dots)
  if (!"keyType" %in% names(call_args)) {
    call_args$keyType <- "uniprot"
  }
  if ("universe" %in% names(dots)) {
    if (glyexp::is_experiment(dots$universe)) {
      call_args$universe <- .extract_uniprot_from_exp(dots$universe)
    }
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

  uniprot_ids <- .extract_uniprot_from_exp(exp)
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
    cli::cli_abort(
      "{.field gene} should not be supplied through `...`; use the function's first argument instead."
    )
  }
  if ("OrgDb" %in% names(dots)) {
    orgdb <- dots$OrgDb
    dots$OrgDb <- NULL
  } else {
    orgdb <- "org.Hs.eg.db"
  }
  orgdb_obj <- .resolve_orgdb(orgdb)

  cli::cli_alert_info("Converting foreground proteins to Entrez IDs")
  fg_genes <- .uniprot_to_entrez(proteins, orgdb_obj)

  # Perform Reactome pathway analysis
  call_args <- c(list(gene = fg_genes), dots)
  if (!"organism" %in% names(call_args)) {
    call_args$organism <- "human"
  }
  if (!"readable" %in% names(call_args)) {
    call_args$readable <- TRUE
  }
  if ("universe" %in% names(dots)) {
    cli::cli_alert_info("Converting background proteins to Entrez IDs")
    if (glyexp::is_experiment(dots$universe)) {
      proteins <- .extract_uniprot_from_exp(dots$universe)
      bg_genes <- .uniprot_to_entrez(proteins, orgdb_obj)
      call_args$universe <- bg_genes
    } else {
      call_args$universe <- .uniprot_to_entrez(dots$universe, orgdb_obj)
    }
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

# Extract UniProt accessions from experiment object helper function
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

.resolve_orgdb <- function(orgdb) {
  if (is.character(orgdb) && length(orgdb) == 1) {
    rlang::check_installed(orgdb)
    get(orgdb, envir = asNamespace(orgdb))
  } else {
    orgdb
  }
}

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
