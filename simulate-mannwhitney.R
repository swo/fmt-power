#!/usr/bin/env Rscript

source("utils.R")

invlogit <- function(lo) 1 / (1 + exp(-lo))

# Simulate one trial
simulate <- function(n, lor) {
  # Simulate donor biomarkers
  x <- rnorm(n)
  # Donor efficacies determined by their biomarkers
  p <- invlogit(x * lor)
  # Simulate patient outcomes
  outcome <- rbinom(n, 1, p)
  # Compare donor biomarkers based on patient outcomes
  wilcox.test(x[outcome], x[!outcome])$p.value
}

simulate_trials <- function(n, lor) {
  map_dbl(1:global_n_trials, ~ simulate(n, lor))
} %>%
  memoise(cache = cache_filesystem("cache/mannwhitney"))

results <- crossing(
  simulation = "mannwhitney",
  n_patients = global_n_patients,
  effect_size = seq(0, 3, length.out = global_n_grid)
) %>%
  mutate(
    n_donors = n_patients,
    p_values = pmap(list(n_patients, effect_size), simulate_trials),
    x = map_dbl(p_values, ~ sum(. <= 0.05)),
    n = map_dbl(p_values, length),
    estimate = x / n
  ) %>%
  select(simulation, n_donors, n_patients, effect_size, x, n, estimate)

write_tsv(results, "results/mannwhitney.tsv")
