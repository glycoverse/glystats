# Real-world glycoproteomics experiment for testing

# `test_gp_exp` is a real-world glycoproteomics experiment for testing.
# It has 12 samples, with 4 groups: "H", "M", "Y", and "C".
# Each group has 3 samples.
# It has 500 variables.
# It has the following columns in `sample_info`:
# - `sample`: sample name
# - `group`: group information
# It has the following columns in `var_info`:
# - `variable`: variable name
# - `proteins`: protein accessions, multiple proteins are separated by ";"
# - `genes`: gene names (symbols), multiple genes are separated by ";"
# - `glycan_composition`: `glyrepr::glycan_composition()` vector
# - `glycan_structure`: `glyrepr::glycan_structure()` vector
# - `protein_sites`: site of glycosylation, multiple sites are separated by ";"
# It has the following metadata fields:
# - `exp_type`: "glycoproteomics",
# - `glycan_type`: "N",
# - `quant_method`: "label-free"
# - `aggr_level`: "gfs"

library(tidyverse)
library(glyexp)
library(glyread)
library(glyclean)

exp <- read_pglyco3_pglycoquant(
  "data-raw/test_gp_res.list",
  glycan_type = "N"
)

set.seed(123)
test_gp_exp <- exp |>
  auto_clean() |>
  slice_sample_var(n = 500) |>
  mutate_obs(
    sample = str_split_i(sample, "-", -1),
    group = str_split_i(sample, "_", 1)
  )
usethis::use_data(test_gp_exp, internal = TRUE, overwrite = TRUE)
