#!/usr/bin/env Rscript

source("utils.R")

# Effect size is σ_D / (σ_D + σ_P). Convert that to σ_D / σ_P as the "between
# variance" in the ANOVA power function.

anova_power <- function(n_donors, patients_per_donor, effect_size) {
  if (patients_per_donor == 1) {
    1.0
  } else {
    power.anova.test(
      groups = n_donors, n = patients_per_donor,
      within.var = 1, between.var = effect_size
    )$power
  }
}

results <- crossing(
  n_patients = global_n_patients,
  n_donors = global_n_donors,
  effect_size = seq(0, 1 - 1e-12, length.out = global_n_grid)
) %>%
  filter(n_patients > n_donors) %>%
  mutate(
    patients_per_donor = n_patients / n_donors,
    between_var = 1 / ((1 / effect_size) - 1),
    estimate = pmap_dbl(
      list(n_donors, patients_per_donor, between_var),
      anova_power
    )
  ) %>%
  select(n_donors, n_patients, effect_size, estimate)

write_tsv(results, "results/anova.tsv")
