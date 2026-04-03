# WikiPathways over-representation analysis (ORA)

Perform WikiPathways ORA for protein UniProt accessions using
[`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html).

- `gly_enrich_wikipathways()` accepts a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  and extracts protein information from the "protein" column in the
  variable information tibble.

- `gly_enrich_wikipathways_()` accepts a character vector of UniProt
  IDs.

As
[`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html)
only accepts Entrez IDs, the UniProt IDs will be first transformed into
Entrez IDs with
[`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html).

## Usage

``` r
gly_enrich_wikipathways(
  exp,
  add_info = TRUE,
  organism = "Homo sapiens",
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)

gly_enrich_wikipathways_(
  proteins,
  organism = "Homo sapiens",
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)
```

## Arguments

- exp:

  (Only for `gly_enrich_wikipathways()`) A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- add_info:

  A logical value. This parameter is included for API consistency but
  has no effect since enrichment results do not contain variable or
  sample columns. Only applicable to top-level APIs.

- organism:

  Passed to `organism` of
  [`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html).
  Species name (e.g., "Homo sapiens", "Mus musculus", "Rattus
  norvegicus"). Defaults to "Homo sapiens".

- orgdb:

  Passed to `OrgDb` of
  [`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html).
  Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to
  "org.Hs.eg.db".

- universe:

  Background genes. If a character vector, it is expected to contain
  UniProt accession IDs; these will be converted to Entrez Gene IDs and
  then passed to `universe` of
  [`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html).
  You can also provide a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with "glycoproteomics" type. In this case all detected proteins
  in this experiment will be extracted as UniProt IDs, converted to
  Entrez IDs, and then passed to `universe` of
  [`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html).

- p_adj_method:

  Passed to `pAdjustMethod` of
  [`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html).

- p_cutoff:

  Passed to `pvalueCutoff` of
  [`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html).

- q_cutoff:

  Passed to `qvalueCutoff` of
  [`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html).

- proteins:

  (Only for `gly_enrich_wikipathways_()`) A character vector of UniProt
  accession IDs.

## Value

A list with three elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `id`: Term ID (e.g., WPXXXXX)

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

- `meta_data` (only for `gly_enrich_wikipathways()`): A list containing
  metadata from the input experiment The list has classes
  `glystats_wikipathways_ora_res` and `glystats_res`.

## Required packages

These functions require the following packages to be installed:

- `clusterProfiler` for enrichment analysis and ID conversion

- `org.Hs.eg.db` for human gene annotation or other OrgDb packages

## See also

[`clusterProfiler::enrichWP()`](https://rdrr.io/pkg/clusterProfiler/man/enrichWP.html)
