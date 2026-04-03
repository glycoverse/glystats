# Network of Cancer Genes (NCG) over-representation analysis (ORA)

Perform NCG ORA for protein UniProt accessions using
[`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html).

- `gly_enrich_ncg()` accepts a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  and extracts protein information from the "protein" column in the
  variable information tibble.

- `gly_enrich_ncg_()` accepts a character vector of UniProt IDs.

As [`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html)
only accepts Entrez IDs, the UniProt IDs will be first transformed into
Entrez IDs with
[`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html).

## Usage

``` r
gly_enrich_ncg(
  exp,
  add_info = TRUE,
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)

gly_enrich_ncg_(
  proteins,
  orgdb = "org.Hs.eg.db",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)
```

## Arguments

- exp:

  (Only for `gly_enrich_ncg()`) A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- add_info:

  A logical value. This parameter is included for API consistency but
  has no effect since enrichment results do not contain variable or
  sample columns. Only applicable to top-level APIs.

- orgdb:

  Passed to `OrgDb` of
  [`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html).
  Organism database name (e.g., "org.Hs.eg.db" for human). Defaults to
  "org.Hs.eg.db".

- universe:

  Background genes. If a character vector, it is expected to contain
  UniProt accession IDs; these will be converted to Entrez Gene IDs and
  then passed to `universe` of
  [`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html).
  You can also provide a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with "glycoproteomics" type. In this case all detected proteins
  in this experiment will be extracted as UniProt IDs, converted to
  Entrez IDs, and then passed to `universe` of
  [`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html).

- p_adj_method:

  Passed to `pAdjustMethod` of
  [`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html).

- p_cutoff:

  Passed to `pvalueCutoff` of
  [`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html).

- q_cutoff:

  Passed to `qvalueCutoff` of
  [`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html).

- proteins:

  (Only for `gly_enrich_ncg_()`) A character vector of UniProt accession
  IDs.

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

- `raw_result`: The raw DOSE enrichResult object

- `meta_data`: A list containing metadata from the input experiment The
  list has classes `glystats_ncg_ora_res` and `glystats_res`.

## Required packages

These functions require the following packages to be installed:

- `clusterProfiler` for ID conversion

- `DOSE` for enrichment analysis

- `org.Hs.eg.db` for human gene annotation or other OrgDb packages

## See also

[`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html)
