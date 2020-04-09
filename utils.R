library(tidyverse)
library(memoise)

set.seed(5)
global_n_grid <- 11
global_n_trials <- 1e3
global_n_donors = c(2, 4, 6, 12),

# Safe test function
chisq_p <- function(x) {
  p <- suppressWarnings(chisq.test(x)$p.value)
  ifelse(is.na(p), 1.0, p)
}

# Simulation functions ------------------------------------------------
# simulate_f :: n_donors, patients_per_donor, effect_size -> p_value

results_base <- crossing(
  n_donors = c(2, 4, 6, 12),
  n_patients = c(12, 24, 48, 96, 192),
) %>%
  filter(n_patients > n_donors) %>%
  mutate(patients_per_donor = n_patients / n_donors)

results_f <- function(
  results_base, simulate_f, lower, upper, cache,
  n_grid = global_n_grid, n_trials = global_n_trials
) {
  simulate_trials <- function(n_donors, patients_per_donor, effect_size) {
    map_dbl(1:n_trials, ~ simulate_f(n_donors, patients_per_donor, effect_size))
  }

  f <- memoise(simulate_trials, cache = cache_filesystem(cache))

  results_base %>%
    crossing(effect_size = seq(lower, upper, length.out = n_grid)) %>%
    mutate(
      p_values = pmap(list(n_donors, patients_per_donor, effect_size), f),
      x = map_dbl(p_values, ~ sum(. <= 0.05)),
      n = map_dbl(p_values, length),
      estimate = x / n
    )
}
