# GO, KEGG, and Reactome over-representation analysis (ORA)

Perform GO, KEGG, and Reactome ORA for protein UniProt accessions. For
glycoproteomics experiments, the function extracts unique proteins from
the variable information.

## Usage

``` r
gly_enrich_go(exp, add_info = TRUE, ...)

gly_enrich_go_(proteins, ...)

gly_enrich_kegg(exp, add_info = TRUE, ...)

gly_enrich_kegg_(proteins, ...)

gly_enrich_reactome(exp, add_info = TRUE, ...)

gly_enrich_reactome_(proteins, ...)
```

## Arguments

- exp:

  (Only for `gly_enrich_go()`, `gly_enrich_kegg()`,
  `gly_enrich_reactome()`) A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- add_info:

  A logical value. This parameter is included for API consistency but
  has no effect since enrichment results do not contain variable or
  sample columns. Only applicable to top-level APIs.

- ...:

  Additional arguments passed to
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html),
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html),
  or
  [`ReactomePA::enrichPathway()`](https://rdrr.io/pkg/ReactomePA/man/enrichPathway.html).
  See the "Additional arguments" section for more information.

- proteins:

  (Only for `gly_enrich_go_()`, `gly_enrich_kegg_()`,
  `gly_enrich_reactome_()`) A character vector of UniProt accession IDs.

## Value

A list with two elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `id`: Term ID (GO:XXXXXXX, hsa:XXXXX, or R-HSA-XXXXX)

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
  classes
  `glystats_go_ora_res`/`glystats_kegg_ora_res`/`glystats_reactome_ora_res`
  and `glystats_res`.

## Required packages

These functions require the following packages to be installed:

- `clusterProfiler` for enrichment analysis

- `ReactomePA` for Reactome pathway analysis

- `org.Hs.eg.db` for human gene annotation (GO analysis only)

## Details

These functions perform over-representation analysis using the specified
database.

`gly_enrich_go()`, `gly_enrich_kegg()`, and `gly_enrich_reactome()` are
the top-level APIs that work with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects and extract protein information automatically from the "protein"
column in the variable information tibble.

`gly_enrich_go_()`, `gly_enrich_kegg_()`, and `gly_enrich_reactome_()`
are the underlying APIs that work with protein vectors directly,
providing more flexibility for users who don't use the glyexp package.

**Gene Extraction (top-level APIs only):** Proteins are extracted from
the experiment's variable information. The function looks for columns
containing protein identifiers and uses them for enrichment analysis.
Protein identifiers should be UniProt accessions.

**GO Analysis:** Uses
[`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)
with UniProt IDs as input.

**KEGG Analysis:** Uses
[`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html)
with UniProt IDs as input.

**Reactome Analysis:** Converts UniProt IDs to Entrez IDs and uses
[`ReactomePA::enrichPathway()`](https://rdrr.io/pkg/ReactomePA/man/enrichPathway.html).

## Additional arguments

`universe` can be passed as a character vector of uniprot accession, or
a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object. If latter, the background proteins will be extracted from the
experiment. This can be convenient when you first perform differential
analysis and then perform enrichment on the differentially expressed
proteins:

    dea_res <- gly_limma(exp)
    sig_exp <- filter_sig_vars(exp, dea_res)
    enrich_res <- gly_enrich_go(sig_exp, universe = exp)

For `gly_enrich_reactome()`, an `OrgDb` can be passed through `...` to
[`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html)
to convert UniProt to Entrez IDs. By default, `org.Hs.eg.db` is used.
This is useful when you want to use a different organism than human:

    enrich_res <- gly_enrich_reactome(exp, OrgDb = "org.Mm.eg.db")

## See also

[`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html),
[`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html),
[`ReactomePA::enrichPathway()`](https://rdrr.io/pkg/ReactomePA/man/enrichPathway.html)
