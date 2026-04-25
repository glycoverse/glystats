# GO over-representation analysis (ORA)

**\[deprecated\]**

This function was deprecated because we decided to move all enrichment
analysis functions to the separate `glyfun` package, which has more
features and better API design.

Perform GO ORA for protein UniProt accessions using
[`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- `gly_enrich_go()` accepts a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  and extracts protein information from the "protein" column in the
  variable information tibble.

- `gly_enrich_go_()` accepts a character vector of UniProt IDs.

## Usage

``` r
gly_enrich_go(
  exp,
  add_info = TRUE,
  orgdb = "org.Hs.eg.db",
  ont = "MF",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)

gly_enrich_go_(
  proteins,
  orgdb = "org.Hs.eg.db",
  keytype = "UNIPROT",
  ont = "MF",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)
```

## Arguments

- exp:

  (Only for `gly_enrich_go()`) A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- add_info:

  A logical value. This parameter is included for API consistency but
  has no effect since enrichment results do not contain variable or
  sample columns. Only applicable to top-level APIs.

- orgdb:

  Passed to `OrgDb` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- ont:

  Passed to `ont` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).
  "BP", "MF", "CC", or "ALL". Defaults to "MF".

- universe:

  Background genes. If a character vector, directly passed to `universe`
  of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).
  You can also provide a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with "glycoproteomics" type. In this case all detected proteins
  in this experiment will be extracted and passed to `universe` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- p_adj_method:

  Passed to `pAdjustMethod` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- p_cutoff:

  Passed to `pvalueCutoff` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- q_cutoff:

  Passed to `qvalueCutoff` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- proteins:

  (Only for `gly_enrich_go_()`) A character vector of UniProt accession
  IDs.

- keytype:

  Passed to `keyType` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

## Value

A list with three elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `id`: Term ID

  - `description`: Term description

  - `gene_ratio`: Ratio of genes in the term to total genes in the input

  - `bg_ratio`: Ratio of genes in the term to total genes in the
    background

  - `p_val`: Raw p-value from hypergeometric test

  - `p_adj`: Adjusted p-value

  - `q_value`: Q-value (FDR)

  - `gene_id`: Gene IDs in the term (separated by "/")

  - `count`: Number of genes in the term

- `raw_result`: The raw clusterProfiler enrichResult object

- `meta_data` (only for `gly_enrich_go()`): A list containing metadata
  from the input experiment The list has classes `glystats_go_ora_res`
  and `glystats_res`.

## Required packages

These functions require the following packages to be installed:

- `clusterProfiler` for enrichment analysis

- `org.Hs.eg.db` for human gene annotation or other OrgDb packages

## See also

[`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)
