#!/usr/bin/env Rscript

source("utils.R")

# Effect size is σ_D / (σ_D + σ_P). Convert that to σ_D / σ_P as the "between
# variance" in the ANOVA power function.

results <- results_base %>%
  crossing(effect_size = seq(0, 1 - 1e-12, length.out = global_n_grid)) %>%
  mutate(
    between_var = 1 / ((1 / effect_size) - 1),
    estimate = pmap_dbl(
      list(n_donors, patients_per_donor, between_var),
      ~ power.anova.test(groups = ..1, n = ..2, within.var = 1, between.var = ..3)$power
    ),
    lci = NA,
    uci = NA
  )

results %>%
  select(n_donors, n_patients, effect_size, estimate, lci, uci) %>%
  write_tsv("results/anova.tsv")
