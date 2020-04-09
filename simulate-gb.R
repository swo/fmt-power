#!/usr/bin/env Rscript

source("utils.R")

# Safe test function
chisq_p <- function(x) {
  p <- suppressWarnings(chisq.test(x)$p.value)
  ifelse(is.na(p), 1.0, p)
}

simulate_base <- function(n_donors, patients_per_donor, phi, p0, p1) {
  quality <- rbinom(n_donors, 1, phi)
  p <- p0 + quality * (p1 - p0)
  success <- rbinom(n_donors, patients_per_donor, p)
  chisq_p(cbind(success, patients_per_donor - success))
}

simulate_f <- function(n_donors, patients_per_donor, p1) {
  p0 <- 1 - p1
  simulate_base(n_donors, patients_per_donor, 0.5, p0, p1)
}


simulate_trials <- function(n_donors, patients_per_donor, effect_size) {
  map_dbl(
    1:global_n_trials,
    ~ simulate_f(n_donors, patients_per_donor, effect_size)
  )
}

f <- memoise(simulate_trials, cache = cache_filesystem("cache/gb"))

results <- crossing(
  n_patients = global_n_patients,
  n_donors = global_n_donors,
  effect_size = seq(0.5, 1.0, length.out = global_n_grid)
) %>%
  mutate(
    patients_per_donor = n_patients / n_donors,
    p_values = pmap(list(n_donors, patients_per_donor, effect_size), f),
    x = map_dbl(p_values, ~ sum(. <= 0.05)),
    n = map_dbl(p_values, length),
    estimate = x / n
  ) %>%
  select(n_donors, n_patients, effect_size, x, n, estimate)

write_tsv(results, "results/gb.tsv")
