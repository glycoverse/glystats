# Disease Ontology over-representation analysis (ORA)

**\[deprecated\]**

This function was deprecated because we decided to move all enrichment
analysis functions to the separate `glyfun` package, which has more
features and better API design.

Perform Disease Ontology ORA for protein UniProt accessions using
[`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).

- `gly_enrich_do()` accepts a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  and extracts protein information from the "protein" column in the
  variable information tibble.

As [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html) only
accepts Entrez IDs, the UniProt IDs will be first transformed into
Entrez IDs with
[`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html).

## Usage

``` r
gly_enrich_do(
  exp,
  add_info = TRUE,
  ont = "HDO",
  organism = "hsa",
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- add_info:

  A logical value. This parameter is included for API consistency but
  has no effect since enrichment results do not contain variable or
  sample columns.

- ont:

  Passed to `ont` of
  [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).
  Disease Ontology type. Options: `"HDO"` (Human Disease Ontology),
  `"HPO"` (Human Phenotype Ontology), `"MPO"` (Mammalian Phenotype
  Ontology). Defaults to `"HDO"`.

- organism:

  Passed to `organism` of
  [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).
  Organism code. `"hsa"` for human (Homo sapiens), `"mmu"` for mouse
  (Mus musculus). Defaults to `"hsa"`.

- orgdb:

  Passed to `OrgDb` of
  [`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html).
  Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to
  "org.Hs.eg.db".

- universe:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  defining the background proteins, or `NULL` to use the default
  background.

- p_adj_method:

  Passed to `pAdjustMethod` of
  [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).

- p_cutoff:

  Passed to `pvalueCutoff` of
  [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).

- q_cutoff:

  Passed to `qvalueCutoff` of
  [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).

## Value

A list with three elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `id`: Term ID (e.g., DOID:XXXX)

  - `description`: Term description

  - `gene_ratio`: Ratio of genes in the term to total genes in the input

  - `bg_ratio`: Ratio of genes in the term to total genes in the
    background

  - `p_val`: Raw p-value from hypergeometric test

  - `p_adj`: Adjusted p-value

  - `q_value`: Q-value (FDR)

  - `gene_id`: Gene IDs in the term (separated by "/")

  - `count`: Number of genes in the term

- `raw_result`: The raw DOSE enrichResult object

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_do_ora_res` and `glystats_res`.

## Required packages

These functions require the following packages to be installed:

- `clusterProfiler` for ID conversion

- `DOSE` for enrichment analysis

- `org.Hs.eg.db` for human gene annotation or other OrgDb packages

## See also

[`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html)
