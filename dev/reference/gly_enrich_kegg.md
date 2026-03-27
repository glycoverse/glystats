# KEGG over-representation analysis (ORA)

Perform KEGG ORA for protein UniProt accessions using
[`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- `gly_enrich_kegg()` accepts a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  and extracts protein information from the "protein" column in the
  variable information tibble.

- `gly_enrich_kegg_()` accepts a character vector of UniProt IDs.

## Usage

``` r
gly_enrich_kegg(
  exp,
  add_info = TRUE,
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)

gly_enrich_kegg_(
  proteins,
  keytype = "uniprot",
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)
```

## Arguments

- exp:

  (Only for `gly_enrich_kegg()`) A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- add_info:

  A logical value. This parameter is included for API consistency but
  has no effect since enrichment results do not contain variable or
  sample columns. Only applicable to top-level APIs.

- organism:

  Passed to `organism` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).
  KEGG organism code (e.g., "hsa" for human, "mmu" for mouse). Defaults
  to "hsa".

- universe:

  Background genes. If a character vector, directly passed to `universe`
  of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).
  You can also provide a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with "glycoproteomics" type. In this case all detected proteins
  in this experiment will be extracted and passed to `universe` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- p_adj_method:

  Passed to `pAdjustMethod` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- p_cutoff:

  Passed to `pvalueCutoff` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- q_cutoff:

  Passed to `qvalueCutoff` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- proteins:

  (Only for `gly_enrich_kegg_()`) A character vector of UniProt
  accession IDs.

- keytype:

  Passed to `keyType` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).
  Defaults to "uniprot".

## Value

A list with two elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `id`: Term ID (e.g., hsa:XXXXX)

  - `description`: Term description

  - `gene_ratio`: Ratio of genes in the term to total genes in the input

  - `bg_ratio`: Ratio of genes in the term to total genes in the
    background

  - `p_val`: Raw p-value from hypergeometric test

  - `p_adj`: Adjusted p-value

  - `q_value`: Q-value (FDR)

  - `gene_id`: Gene IDs in the term (separated by "/")

  - `count`: Number of genes in the term

- `raw_result`: The raw clusterProfiler enrichResult object The list has
  classes `glystats_kegg_ora_res` and `glystats_res`.

## Required packages

These functions require the following packages to be installed:

- `clusterProfiler` for enrichment analysis

## See also

[`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html)
