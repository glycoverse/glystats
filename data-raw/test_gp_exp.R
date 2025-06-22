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
    group = str_split_i(sample, "_", -1)
  )
usethis::use_data(test_gp_exp, internal = TRUE, overwrite = TRUE)
