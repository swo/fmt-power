library(tidyverse)
library(memoise)

set.seed(5)

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

results_f <- function(simulate_f, lower, upper, cache, n_grid = 10, n_trials = 1e2) {
  simulate_trials <- function(n_trials, n_donors, patients_per_donor, effect_size) {
    map_dbl(1:n_trials, ~ simulate_f(n_donors, patients_per_donor, effect_size))
  }

  f <- memoise(simulate_trials, cache = cache_filesystem(cache))

  results_base %>%
    crossing(effect_size = seq(lower, upper, length.out = n_grid)) %>%
    mutate(
      p_values = pmap(list(n_trials, n_donors, patients_per_donor, effect_size), f),
      x = map_dbl(p_values, ~ sum(. <= 0.05)),
      n = map_dbl(p_values, length),
      test = map2(x, n, binom.test),
      estimate = x / n,
      lci = map_dbl(test, ~ .$conf.int[1]),
      uci = map_dbl(test, ~ .$conf.int[2])
    )
}

plot_f <- function(results) {
  results %>%
    mutate_at(c("n_patients"), ~ fct_rev(factor(.))) %>%
    ggplot(aes(effect_size, estimate)) +
    facet_grid(n_patients ~ n_donors) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    geom_hline(yintercept = c(0, 1)) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = "gray") +
    geom_line(size = 1.5) +
    scale_y_continuous(
      name = "Statistical power",
      labels = scales::percent,
      breaks = c(0, 0.5, 0.8, 1)
    ) +
    cowplot::theme_half_open()
}
